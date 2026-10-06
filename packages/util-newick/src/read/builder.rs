use crate::grammar::Rule;
use crate::model::comment::{EdgeComment, EdgeField, LabelSide, NewickComment, NodeComment, ValueSide};
use crate::model::data::{NewickEdgeData, NewickHybrid, NewickNodeData, SupportSource};
use crate::model::graph::{NewickEdgeEntry, NewickGraph};
use crate::read::comments::{is_comment_rule, read_comment};
use crate::read::context::MapContext;
use crate::read::error::{Location, NewickError, NewickErrorKind};
use crate::read::labels::{LabelToken, read_label};
use crate::read::options::InternalLabel;
use pest::iterators::Pair;
use std::collections::BTreeMap;
use std::mem;

pub(crate) fn build_graph<'i>(
  tokens: impl Iterator<Item = Pair<'i, Rule>>,
  context: &mut MapContext<'_, 'i>,
  start: Location,
) -> Result<NewickGraph, NewickError> {
  let mut builder = Builder {
    context,
    nodes: Vec::new(),
    edges: Vec::new(),
    hybrids: BTreeMap::new(),
    open: Vec::new(),
    slot: Slot::default(),
    rooted: None,
    weight: None,
  };
  for token in tokens {
    builder.consume(token)?;
  }
  builder.finish(start)
}

struct Builder<'c, 'o, 'i> {
  context: &'c mut MapContext<'o, 'i>,
  nodes: Vec<NewickNodeData>,
  edges: Vec<NewickEdgeEntry>,
  hybrids: BTreeMap<(Option<String>, u32), HybridEntry>,
  open: Vec<OpenFrame<'i>>,
  slot: Slot,
  rooted: Option<bool>,
  weight: Option<f64>,
}

impl<'i> Builder<'_, '_, 'i> {
  fn consume(&mut self, token: Pair<'i, Rule>) -> Result<(), NewickError> {
    match token.as_rule() {
      Rule::rooting => self.set_rooting(&token),
      Rule::weight => self.set_weight(&token),
      Rule::open => {
        let comments = mem::take(&mut self.slot.before_label);
        self.open.push(OpenFrame {
          token,
          children: Vec::new(),
          comments,
        });
        Ok(())
      },
      Rule::comma => {
        let child = self.finish_slot()?;
        let Some(frame) = self.open.last_mut() else {
          return Err(self.syntax_error(&token, "unexpected ',' outside of parentheses"));
        };
        frame.children.push(child);
        Ok(())
      },
      Rule::close => {
        let child = self.finish_slot()?;
        let Some(mut frame) = self.open.pop() else {
          return Err(self.syntax_error(&token, "unexpected ')' without a matching '('"));
        };
        frame.children.push(child);
        self.slot.children = Some(frame.children);
        self.slot.before_label = frame.comments;
        Ok(())
      },
      Rule::label | Rule::network_label => {
        self.slot.label_location = Some(self.context.locate(&token));
        self.slot.label = Some(read_label(token, self.context)?);
        Ok(())
      },
      Rule::colon => {
        self.slot.fields.push(FieldSlot::default());
        Ok(())
      },
      Rule::number => self.set_field_value(&token),
      rule if is_comment_rule(rule) => {
        let comment = read_comment(token, self.context)?;
        self.slot.push_comment(comment);
        Ok(())
      },
      _ => Ok(()),
    }
  }

  fn set_rooting(&mut self, token: &Pair<'i, Rule>) -> Result<(), NewickError> {
    if self.rooted.is_some() {
      self.context.tolerate(
        NewickErrorKind::Structure,
        token,
        "The tree has more than one rooting comment",
      )?;
    }
    let value = token.clone().into_inner().as_str();
    self.rooted = Some(value.eq_ignore_ascii_case("r"));
    Ok(())
  }

  fn set_weight(&mut self, token: &Pair<'i, Rule>) -> Result<(), NewickError> {
    if self.weight.is_some() {
      self.context.tolerate(
        NewickErrorKind::Structure,
        token,
        "The tree has more than one weight comment",
      )?;
    }
    let text = token.clone().into_inner().as_str();
    self.weight = Some(self.parse_number(text, token)?);
    Ok(())
  }

  fn set_field_value(&mut self, token: &Pair<'i, Rule>) -> Result<(), NewickError> {
    let value = self.parse_number(token.as_str(), token)?;
    match self.slot.fields.last_mut() {
      Some(field) => field.value = Some(value),
      None => return Err(self.syntax_error(token, "a number outside of a ':' field")),
    }
    Ok(())
  }

  fn parse_number(&self, text: &str, token: &Pair<'i, Rule>) -> Result<f64, NewickError> {
    text.parse::<f64>().map_err(|error| {
      self.context.error(
        NewickErrorKind::Syntax,
        token,
        format!("{text:?} is not a number: {error}"),
      )
    })
  }

  fn syntax_error(&self, token: &Pair<'i, Rule>, message: &str) -> NewickError {
    self.context.error(NewickErrorKind::Syntax, token, message)
  }

  fn finish(mut self, start: Location) -> Result<NewickGraph, NewickError> {
    if let Some(frame) = self.open.last() {
      return Err(self.syntax_error(&frame.token, "the '(' is never closed"));
    }
    let (root, root_edge) = self.finish_slot()?;
    let mut graph = NewickGraph::from_parts(self.nodes, self.edges, root);
    graph.set_root_edge(root_edge);
    graph.set_rooted(self.rooted);
    graph.set_weight(self.weight);
    graph
      .validate()
      .map_err(|error| NewickError::new(NewickErrorKind::Structure, start, format!("{error:#}")))?;
    Ok(graph)
  }

  fn finish_slot(&mut self) -> Result<(usize, NewickEdgeData), NewickError> {
    let slot = mem::take(&mut self.slot);
    let is_internal = slot.children.is_some();
    let mut edge = NewickEdgeData::new();
    let field_support = self.read_fields(slot.fields, &mut edge);
    let mut node = NewickNodeData::new();
    let mut hybrid = None;
    if let Some(label) = slot.label {
      let location = slot.label_location.unwrap_or(Location::START);
      let decoded = self.decode_label(label, is_internal, field_support, location)?;
      if let Some(name) = decoded.name {
        node = node.with_name(name);
      }
      if let Some(support) = decoded.support {
        edge.set_support(support, SupportSource::Label);
      }
      if let Some(tag) = decoded.hybrid {
        node = node.with_hybrid(tag.clone());
        edge = edge.with_acceptor(decoded.is_acceptor);
        hybrid = Some((tag, location));
      }
    }
    let has_label = node.has_label();
    let side = |side: LabelSide| if has_label { side } else { LabelSide::AfterLabel };
    let comments = slot
      .before_label
      .into_iter()
      .map(|comment| NodeComment::new(side(LabelSide::BeforeLabel), comment))
      .chain(
        slot
          .after_label
          .into_iter()
          .map(|comment| NodeComment::new(LabelSide::AfterLabel, comment)),
      );
    node.comments_mut().extend(comments);
    let children = slot.children.unwrap_or_default();
    let idx = if let Some((tag, location)) = hybrid {
      self.merge_hybrid(tag, node, !children.is_empty(), location)?
    } else {
      self.nodes.push(node);
      self.nodes.len() - 1
    };
    for (child, child_edge) in children {
      self.edges.push(NewickEdgeEntry::new(idx, child, child_edge));
    }
    Ok((idx, edge))
  }

  fn read_fields(&self, fields: Vec<FieldSlot>, edge: &mut NewickEdgeData) -> bool {
    let mut field_support = false;
    for (field, slot) in EdgeField::ALL.into_iter().zip(fields) {
      match (field, slot.value) {
        (EdgeField::Length, value) => edge.set_branch_length(value),
        (EdgeField::Support, Some(value)) => {
          edge.set_support(vec![value], SupportSource::Field);
          field_support = true;
        },
        (EdgeField::Probability, value) => edge.set_probability(value),
        (EdgeField::Support, None) => {},
      }
      let comments = slot
        .before
        .into_iter()
        .map(|comment| EdgeComment::new(field, ValueSide::BeforeValue, comment))
        .chain(
          slot
            .after
            .into_iter()
            .map(|comment| EdgeComment::new(field, ValueSide::AfterValue, comment)),
        );
      edge.comments_mut().extend(comments);
    }
    field_support
  }

  fn decode_label(
    &mut self,
    label: LabelToken,
    is_internal: bool,
    field_support: bool,
    location: Location,
  ) -> Result<DecodedLabel, NewickError> {
    let mut decoded = DecodedLabel {
      name: label.name,
      support: None,
      hybrid: label.hybrid,
      is_acceptor: label.is_acceptor,
    };
    let reads_support = is_internal && !field_support && self.context.options.internal_label != InternalLabel::Name;
    match label.support {
      Some(support) if reads_support => {
        decoded.support = Some(support);
        decoded.name = None;
      },
      None if is_internal && self.context.options.internal_label == InternalLabel::Support => {
        let message = format!("The internal label {:?} is not a support value", label.text);
        self
          .context
          .tolerate_at(NewickErrorKind::Structure, location, message)?;
      },
      Some(_) | None => {},
    }
    Ok(decoded)
  }

  fn merge_hybrid(
    &mut self,
    hybrid: NewickHybrid,
    node: NewickNodeData,
    has_children: bool,
    location: Location,
  ) -> Result<usize, NewickError> {
    let tag = hybrid.tag(false);
    let NewickHybrid { kind, index } = hybrid;
    let key = (kind, index);
    let Some(entry) = self.hybrids.get_mut(&key) else {
      self.nodes.push(node);
      let idx = self.nodes.len() - 1;
      self.hybrids.insert(key, HybridEntry { idx, has_children });
      return Ok(idx);
    };
    if has_children && entry.has_children {
      return Err(NewickError::new(
        NewickErrorKind::Structure,
        location,
        format!("The hybrid node {tag} has children in more than one of its occurrences"),
      ));
    }
    entry.has_children |= has_children;
    let idx = entry.idx;
    let existing = &mut self.nodes[idx];
    match (existing.name(), node.name()) {
      (_, None) => {},
      (None, Some(name)) => existing.set_name(Some(name.to_owned())),
      (Some(previous), Some(name)) if previous == name => {},
      (Some(previous), Some(name)) => {
        return Err(NewickError::new(
          NewickErrorKind::Structure,
          location,
          format!("The occurrences of the hybrid node {tag} have different names: {previous:?} and {name:?}"),
        ));
      },
    }
    existing.comments_mut().extend(node.into_comments());
    Ok(idx)
  }
}

#[derive(Default)]
struct Slot {
  before_label: Vec<NewickComment>,
  after_label: Vec<NewickComment>,
  label: Option<LabelToken>,
  label_location: Option<Location>,
  fields: Vec<FieldSlot>,
  children: Option<Vec<(usize, NewickEdgeData)>>,
}

impl Slot {
  fn push_comment(&mut self, comment: NewickComment) {
    if let Some(field) = self.fields.last_mut() {
      match field.value {
        Some(_) => field.after.push(comment),
        None => field.before.push(comment),
      }
    } else if self.label.is_some() {
      self.after_label.push(comment);
    } else {
      self.before_label.push(comment);
    }
  }
}

#[derive(Default)]
struct FieldSlot {
  value: Option<f64>,
  before: Vec<NewickComment>,
  after: Vec<NewickComment>,
}

struct OpenFrame<'i> {
  token: Pair<'i, Rule>,
  children: Vec<(usize, NewickEdgeData)>,
  comments: Vec<NewickComment>,
}

struct HybridEntry {
  idx: usize,
  has_children: bool,
}

struct DecodedLabel {
  name: Option<String>,
  support: Option<Vec<f64>>,
  hybrid: Option<NewickHybrid>,
  is_acceptor: bool,
}
