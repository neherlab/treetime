use crate::dialect::NewickAnnotations;
use crate::grammar::{Rule, matches, parse};
use crate::model::comment::{EdgeComment, LabelSide, NewickComment, NodeComment, ValueSide};
use crate::model::data::{NewickEdgeData, SupportSource};
use crate::model::graph::NewickGraph;
use crate::model::validate::describe_node;
use crate::model::value::NewickValue;
use crate::write::comments::{encode_beast, encode_comment, encode_nhx};
use crate::write::conversions::{Conversion, DataKind, conversion};
use crate::write::options::{BranchAnnotations, NewickWriteOptions, Quoting, Spaces, SupportPlacement};
use eyre::{Report, WrapErr, eyre};
use std::borrow::Cow;
use std::collections::BTreeMap;
use std::{io, mem};

pub fn newick_to_writer(
  writer: &mut impl io::Write,
  graph: &NewickGraph,
  options: &NewickWriteOptions,
) -> Result<(), Report> {
  let text = newick_to_string(graph, options)?;
  writer.write_all(text.as_bytes())?;
  Ok(())
}

pub fn newick_to_string(graph: &NewickGraph, options: &NewickWriteOptions) -> Result<String, Report> {
  let mut text = String::new();
  write_tree(&mut text, graph, options, None).wrap_err("When writing Newick")?;
  Ok(text)
}

pub fn write_newick_trees<'g>(
  writer: &mut impl io::Write,
  trees: impl IntoIterator<Item = &'g NewickGraph>,
  options: &NewickWriteOptions,
) -> Result<(), Report> {
  for (idx, graph) in trees.into_iter().enumerate() {
    let mut text = newick_to_string(graph, options).wrap_err_with(|| format!("When writing tree {}", idx + 1))?;
    text.push('\n');
    writer.write_all(text.as_bytes())?;
  }
  Ok(())
}

pub(crate) fn write_tree(
  out: &mut String,
  graph: &NewickGraph,
  options: &NewickWriteOptions,
  translate: Option<&BTreeMap<String, String>>,
) -> Result<(), Report> {
  graph.validate()?;
  let writer = TreeWriter {
    graph,
    options,
    translate,
  };
  writer.check()?;
  writer.write_tree_comments(out)?;
  writer.write_nodes(out)?;
  out.push(';');
  Ok(())
}

struct TreeWriter<'g> {
  graph: &'g NewickGraph,
  options: &'g NewickWriteOptions,
  translate: Option<&'g BTreeMap<String, String>>,
}

impl TreeWriter<'_> {
  fn check(&self) -> Result<(), Report> {
    let dialect = self.options.dialect;
    if let Some((node, _)) = self.graph.nodes().find(|(_, data)| data.hybrid().is_some())
      && conversion(dialect, DataKind::HybridNodes) == Conversion::Fail
    {
      return Err(eyre!(
        "The {dialect} dialect cannot hold hybrid nodes such as {}; write networks with the enewick or rich structure",
        describe_node(self.graph, node)
      ));
    }
    match &self.options.support {
      SupportPlacement::Field if dialect.structure.field_count() < 2 => Err(eyre!(
        "Support in a colon field needs the rich structure, not {dialect}"
      )),
      SupportPlacement::Annotation(_)
        if !matches!(dialect.annotations, NewickAnnotations::Beast | NewickAnnotations::Nhx) =>
      {
        Err(eyre!(
          "Support in an annotation needs beast or nhx annotations, not {dialect}"
        ))
      },
      SupportPlacement::Source
      | SupportPlacement::Label
      | SupportPlacement::Field
      | SupportPlacement::Annotation(_) => Ok(()),
    }
  }

  fn write_tree_comments(&self, out: &mut String) -> Result<(), Report> {
    let dialect = self.options.dialect;
    if conversion(dialect, DataKind::Rooting) == Conversion::Keep {
      match self.graph.rooted() {
        Some(true) => out.push_str("[&R]"),
        Some(false) => out.push_str("[&U]"),
        None => {},
      }
    }
    if conversion(dialect, DataKind::Weight) == Conversion::Keep
      && let Some(weight) = self.graph.weight()
    {
      out.push_str("[&W ");
      out.push_str(&self.options.numbers.format(weight)?);
      out.push(']');
    }
    Ok(())
  }

  fn write_nodes(&self, out: &mut String) -> Result<(), Report> {
    let graph = self.graph;
    let mut defined = vec![false; graph.node_count()];
    defined[graph.root()] = true;
    let mut stack = vec![Frame {
      node: graph.root(),
      edge: None,
      next_child: 0,
      is_definition: true,
    }];
    while let Some(depth) = (!stack.is_empty()).then_some(stack.len()) {
      let Some(frame) = stack.last_mut() else {
        break;
      };
      let children: &[usize] = if frame.is_definition {
        graph.child_edges(frame.node)
      } else {
        &[]
      };
      if let Some(&edge) = children.get(frame.next_child) {
        out.push(if frame.next_child == 0 { '(' } else { ',' });
        self.indent(out, depth);
        frame.next_child += 1;
        let child = graph.edge(edge).child();
        let is_definition = graph.node(child).hybrid().is_none() || !mem::replace(&mut defined[child], true);
        stack.push(Frame {
          node: child,
          edge: Some(edge),
          next_child: 0,
          is_definition,
        });
        continue;
      }
      let finished = *frame;
      stack.pop();
      if finished.next_child > 0 {
        self.indent(out, depth - 1);
        out.push(')');
      }
      self.write_node(out, finished)?;
    }
    Ok(())
  }

  fn indent(&self, out: &mut String, depth: usize) {
    if let Some(width) = self.options.indent {
      out.push('\n');
      out.push_str(&" ".repeat(width * depth));
    }
  }

  fn write_node(&self, out: &mut String, frame: Frame) -> Result<(), Report> {
    let empty = NewickEdgeData::new();
    let edge = match frame.edge {
      Some(edge) => self.graph.edge(edge).data(),
      None if self.options.root_edge => self.graph.root_edge(),
      None => &empty,
    };
    let node_context = || format!("When writing {}", describe_node(self.graph, frame.node));
    let branch_context = || {
      format!(
        "When writing the branch above {}",
        describe_node(self.graph, frame.node)
      )
    };
    let support = self.support_target(frame, edge).wrap_err_with(branch_context)?;
    let comments = self.label_comments(frame);
    let node = self.graph.node(frame.node);
    let starts_tree = frame.node == self.graph.root() && frame.next_child == 0;
    let has_label = node.name().is_some() || node.hybrid().is_some();
    self
      .write_label_comments(out, comments, LabelSide::BeforeLabel, starts_tree)
      .wrap_err_with(node_context)?;
    self
      .write_label(out, frame, edge, &support)
      .wrap_err_with(node_context)?;
    self
      .write_label_comments(out, comments, LabelSide::AfterLabel, starts_tree && !has_label)
      .wrap_err_with(node_context)?;
    if let SupportTarget::Annotation(key, values) = &support {
      let value = match values.as_slice() {
        [single] => NewickValue::Number(*single),
        _ => NewickValue::Array(
          values
            .iter()
            .map(|&value| NewickValue::Number(value))
            .collect::<Vec<_>>()
            .into(),
        ),
      };
      let pairs = [(key.clone(), value)];
      let encoded = match self.options.dialect.annotations {
        NewickAnnotations::Nhx => encode_nhx(&pairs),
        NewickAnnotations::Beast | NewickAnnotations::Plain | NewickAnnotations::MrBayes => encode_beast(&pairs),
      };
      out.push_str(&encoded.wrap_err_with(branch_context)?);
    }
    let field_support = match support {
      SupportTarget::Field(value) => Some(value),
      SupportTarget::None | SupportTarget::Label(_) | SupportTarget::Annotation(..) => None,
    };
    self
      .write_fields(out, edge, field_support)
      .wrap_err_with(branch_context)
  }

  fn support_target<'e>(&self, frame: Frame, edge: &'e NewickEdgeData) -> Result<SupportTarget<'e>, Report> {
    let values = edge.support();
    if values.is_empty() {
      return Ok(SupportTarget::None);
    }
    let target = match (&self.options.support, edge.support_source()) {
      (SupportPlacement::Source, SupportSource::Label) | (SupportPlacement::Label, _) => SupportTarget::Label(values),
      (SupportPlacement::Source, SupportSource::Field) | (SupportPlacement::Field, _) => {
        if conversion(self.options.dialect, DataKind::FieldSupport) == Conversion::Drop {
          return Ok(SupportTarget::None);
        }
        match values {
          [single] => SupportTarget::Field(*single),
          _ => {
            return Err(eyre!(
              "A colon field holds one support value, but the branch has {}",
              values.len()
            ));
          },
        }
      },
      (SupportPlacement::Annotation(key), _) => SupportTarget::Annotation(key.clone(), values.to_vec()),
    };
    if let SupportTarget::Label(_) = target {
      let node = self.graph.node(frame.node);
      if self.graph.is_leaf(frame.node) {
        return Err(eyre!(
          "A leaf cannot carry support in its label, because a leaf label is read as a name"
        ));
      }
      if node.name().is_some() || node.hybrid().is_some() {
        return Err(eyre!(
          "The label of the node holds its name, so it cannot also hold the support of the branch above"
        ));
      }
    }
    Ok(target)
  }

  fn write_label(
    &self,
    out: &mut String,
    frame: Frame,
    edge: &NewickEdgeData,
    support: &SupportTarget<'_>,
  ) -> Result<(), Report> {
    let node = self.graph.node(frame.node);
    let tag = node.hybrid().map(|hybrid| hybrid.tag(edge.is_acceptor()));
    if let SupportTarget::Label(values) = support {
      let texts = values
        .iter()
        .map(|&value| self.options.numbers.format(value))
        .collect::<Result<Vec<_>, _>>()?;
      out.push_str(&texts.join("/"));
    }
    if let Some(name) = node.name() {
      let is_leaf = self.graph.is_leaf(frame.node);
      let translated = self.translate.filter(|_| is_leaf).and_then(|table| table.get(name));
      match translated {
        Some(key) => out.push_str(key),
        None => self.push_name(out, name, !is_leaf, tag.as_deref()),
      }
    }
    if let Some(tag) = tag {
      out.push_str(&tag);
    }
    Ok(())
  }

  fn push_name(&self, out: &mut String, name: &str, is_internal: bool, tag: Option<&str>) {
    let written: Cow<'_, str> = match self.options.spaces {
      Spaces::Quote => Cow::Borrowed(name),
      Spaces::Underscore => Cow::Owned(name.replace(' ', "_")),
    };
    let needs_quotes = self.options.quoting == Quoting::Always
      || !matches(Rule::safe_name_exact, &written)
      || is_internal && matches(Rule::support_label_exact, &written)
      || tag.is_some_and(|tag| !splits_back(&written, tag));
    if needs_quotes {
      out.push('\'');
      out.push_str(&name.replace('\'', "''"));
      out.push('\'');
    } else {
      out.push_str(&written);
    }
  }

  fn label_comments(&self, frame: Frame) -> &[NodeComment] {
    let node = self.graph.node(frame.node);
    if node.hybrid().is_none() {
      return node.comments();
    }
    match frame.edge {
      Some(edge) => self.graph.edge(edge).data().occurrence_comments(),
      None => self.graph.root_edge().occurrence_comments(),
    }
  }

  fn write_label_comments(
    &self,
    out: &mut String,
    comments: &[NodeComment],
    side: LabelSide,
    starts_tree: bool,
  ) -> Result<(), Report> {
    for comment in comments.iter().filter(|comment| comment.position == side) {
      if starts_tree {
        self.check_tree_start_comment(&comment.comment)?;
      }
      if let Some(text) = encode_comment(&comment.comment, self.options.dialect)? {
        out.push_str(&text);
      }
    }
    Ok(())
  }

  fn check_tree_start_comment(&self, comment: &NewickComment) -> Result<(), Report> {
    let NewickComment::Plain(text) = comment else {
      return Ok(());
    };
    let bracketed = format!("[{text}]");
    let reads_as = [
      (Rule::rooting_exact, DataKind::RootingPlainComments, "rooting"),
      (Rule::weight_exact, DataKind::WeightPlainComments, "tree weight"),
    ]
    .into_iter()
    .find(|&(rule, data, _)| conversion(self.options.dialect, data) == Conversion::Fail && matches(rule, &bracketed));
    match reads_as {
      Some((_, _, what)) => Err(eyre!(
        "The comment {bracketed} at the start of the tree cannot be written in the {} dialect, which would read it as the {what} comment",
        self.options.dialect
      )),
      None => Ok(()),
    }
  }

  fn write_fields(&self, out: &mut String, edge: &NewickEdgeData, field_support: Option<f64>) -> Result<(), Report> {
    let mut fields: [FieldText; 3] = Default::default();
    if let Some(length) = edge.branch_length() {
      fields[0].value = Some(self.options.numbers.format(length)?);
    }
    if let Some(support) = field_support {
      fields[1].value = Some(self.options.numbers.format(support)?);
    }
    if let Some(probability) = edge.probability()
      && conversion(self.options.dialect, DataKind::Probability) == Conversion::Keep
    {
      fields[2].value = Some(self.options.numbers.format(probability)?);
    }
    for comment in edge.comments() {
      let index = comment.field.index();
      if index > 0 && conversion(self.options.dialect, DataKind::FieldComments) == Conversion::Drop {
        continue;
      }
      if let Some(text) = encode_comment(&comment.comment, self.options.dialect)? {
        let field = &mut fields[index];
        match self.side(comment) {
          ValueSide::BeforeValue => field.before.push_str(&text),
          ValueSide::AfterValue => field.after.push_str(&text),
        }
      }
    }
    let Some(last) = fields.iter().rposition(|field| !field.is_empty()) else {
      return Ok(());
    };
    for field in fields.iter().take(last + 1) {
      out.push(':');
      out.push_str(&field.before);
      out.push_str(field.value.as_deref().unwrap_or_default());
      out.push_str(&field.after);
    }
    Ok(())
  }

  fn side(&self, comment: &EdgeComment) -> ValueSide {
    let is_annotation = !matches!(comment.comment, NewickComment::Plain(_));
    match self.options.branch_annotations {
      BranchAnnotations::BeforeLength if is_annotation && comment.field.index() == 0 => ValueSide::BeforeValue,
      BranchAnnotations::AfterLength if is_annotation && comment.field.index() == 0 => ValueSide::AfterValue,
      BranchAnnotations::Recorded | BranchAnnotations::BeforeLength | BranchAnnotations::AfterLength => comment.side,
    }
  }
}

fn splits_back(name: &str, tag: &str) -> bool {
  let text = format!("{name}{tag}");
  let Ok(pairs) = parse(Rule::network_label_exact, &text) else {
    return false;
  };
  let parts: Vec<(Rule, &str)> = pairs
    .flatten()
    .filter(|pair| matches!(pair.as_rule(), Rule::unquoted_label | Rule::hybrid_tag))
    .map(|pair| (pair.as_rule(), pair.as_str()))
    .collect();
  parts == [(Rule::unquoted_label, name), (Rule::hybrid_tag, tag)]
}

#[derive(Clone, Copy)]
struct Frame {
  node: usize,
  edge: Option<usize>,
  next_child: usize,
  is_definition: bool,
}

enum SupportTarget<'e> {
  None,
  Label(&'e [f64]),
  Field(f64),
  Annotation(String, Vec<f64>),
}

#[derive(Default)]
struct FieldText {
  value: Option<String>,
  before: String,
  after: String,
}

impl FieldText {
  fn is_empty(&self) -> bool {
    self.value.is_none() && self.before.is_empty() && self.after.is_empty()
  }
}
