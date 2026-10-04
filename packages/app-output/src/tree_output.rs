use crate::mutation_filter::UnknownBridge;
use chrono::Utc;
use eyre::{Report, WrapErr};
use maplit::{btreemap, btreeset};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::path::PathBuf;
use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
use treetime::progress::LogSink;
use treetime::progress_warn;
use treetime::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub, mutation_event_strings};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::auspice::auspice_write_file;
use treetime_io::auspice_types::{
  AuspiceColoring, AuspiceDisplayDefaults, AuspiceGenomeAnnotationCds, AuspiceGenomeAnnotationNuc,
  AuspiceGenomeAnnotations, AuspiceNumDate, AuspiceTree, AuspiceTreeBranchAttrs, AuspiceTreeBranchAttrsLabels,
  AuspiceTreeData, AuspiceTreeMeta, AuspiceTreeNode, AuspiceTreeNodeAttr, AuspiceTreeNodeAttrs,
};
use treetime_io::graph::TreeWriteKind;
use treetime_io::graphviz::graphviz_write_file;
use treetime_io::nex::{NexWriteOptions, nex_write_file_with};
use treetime_io::nwk::{CommentProviders, NwkWriteOptions, nwk_write_file_with, nwk_write_str};
use treetime_io::usher_mat::{
  UsherMatJsonOptions, UsherMetadata, UsherMutation, UsherMutationList, UsherTree, UsherTreeNode,
  usher_mat_json_write_file, usher_mat_pb_write_file,
};
use treetime_primitives::AsciiChar;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use treetime_utils::{make_error, make_internal_error, make_internal_report};

pub(crate) const COLORING_BAD_BRANCH: &str = "bad_branch";
const COLORING_GENOTYPE: &str = "gt";
pub(crate) const COLORING_NUM_DATE: &str = "num_date";
pub(crate) const NUC_TRACK: &str = "nuc";
pub(crate) const PANEL_ENTROPY: &str = "entropy";
pub(crate) const PANEL_TREE: &str = "tree";

pub(crate) fn write_tree_outputs<A, M>(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  graphviz_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
  command: &str,
  to_auspice: A,
  to_mat: M,
  log: &dyn LogSink,
) -> Result<(), Report>
where
  A: Fn() -> Result<AuspiceTree, Report>,
  M: Fn() -> Result<MatOutput, Report>,
{
  let mut mat = None;
  for (kind, path) in outputs {
    match kind {
      TreeWriteKind::Nwk(spec) => {
        graph
          .get_exactly_one_root()
          .wrap_err_with(|| format!("When converting {command} graph to Newick"))?;
        nwk_write_file_with(
          path,
          graph,
          names,
          nwk_weights,
          &NwkWriteOptions {
            style: spec.style,
            ..NwkWriteOptions::default()
          },
          providers,
        )?;
      },
      TreeWriteKind::Nexus(spec) => {
        graph
          .get_exactly_one_root()
          .wrap_err_with(|| format!("When converting {command} graph to Nexus"))?;
        nex_write_file_with(
          path,
          graph,
          names,
          nwk_weights,
          &NexWriteOptions {
            style: spec.style,
            ..NexWriteOptions::default()
          },
          providers,
        )?;
      },
      TreeWriteKind::Auspice => auspice_write_file(path, &to_auspice()?)?,
      TreeWriteKind::MatPb => usher_mat_pb_write_file(path, converted_mat(&mut mat, &to_mat, log)?)?,
      TreeWriteKind::MatJson => {
        usher_mat_json_write_file(
          path,
          converted_mat(&mut mat, &to_mat, log)?,
          &UsherMatJsonOptions::default(),
        )?;
      },
      TreeWriteKind::GraphJson => {
        json_write_file(path, graph, JsonPretty(true))?;
      },
      TreeWriteKind::Dot => graphviz_write_file(path, graph, names, graphviz_weights)?,
    }
  }
  Ok(())
}

fn converted_mat<'a, M>(mat: &'a mut Option<UsherTree>, to_mat: &M, log: &dyn LogSink) -> Result<&'a UsherTree, Report>
where
  M: Fn() -> Result<MatOutput, Report>,
{
  let tree = if let Some(tree) = mat.take() {
    tree
  } else {
    let MatOutput { tree, gaps } = to_mat()?;
    if let Some(warning) = gaps.warning() {
      progress_warn!(log, "{warning}");
    }
    tree
  };
  Ok(mat.insert(tree))
}

pub(crate) fn auspice_data(
  title: &str,
  updated: &str,
  mut colorings: Vec<AuspiceColoring>,
  filters: Vec<String>,
  color_by: Option<String>,
  cdses: BTreeMap<String, AuspiceGenomeAnnotationCds>,
  root_sequences: Option<BTreeMap<String, String>>,
  has_mutations: bool,
) -> Result<AuspiceTreeData, Report> {
  if has_mutations {
    colorings.push(coloring(COLORING_GENOTYPE, "Genotype", "categorical"));
  }
  let genome_annotations = genome_annotations(root_sequences.as_ref(), cdses)?;
  let mut panels = vec![PANEL_TREE.to_owned()];
  if genome_annotations.is_some() {
    panels.push(PANEL_ENTROPY.to_owned());
  }
  Ok(AuspiceTreeData {
    version: Some("v2".to_owned()),
    meta: AuspiceTreeMeta {
      title: Some(title.to_owned()),
      updated: Some(updated.to_owned()),
      panels,
      genome_annotations,
      colorings,
      filters,
      display_defaults: AuspiceDisplayDefaults {
        color_by,
        ..AuspiceDisplayDefaults::default()
      },
      ..AuspiceTreeMeta::default()
    },
    root_sequence: root_sequences.filter(|sequences| !sequences.is_empty()),
    other: Value::default(),
  })
}

fn genome_annotations(
  root_sequences: Option<&BTreeMap<String, String>>,
  cdses: BTreeMap<String, AuspiceGenomeAnnotationCds>,
) -> Result<Option<AuspiceGenomeAnnotations>, Report> {
  let nuc = root_sequences
    .and_then(|sequences| sequences.get(NUC_TRACK))
    .map(|sequence| -> Result<_, Report> {
      Ok(AuspiceGenomeAnnotationNuc {
        start: 1,
        end: isize::try_from(sequence.len()).wrap_err("Nucleotide sequence length does not fit Auspice coordinates")?,
        strand: Some("+".to_owned()),
        r#type: Some("source".to_owned()),
        other: Value::default(),
      })
    })
    .transpose()?;
  if nuc.is_none() && cdses.is_empty() {
    return Ok(None);
  }
  Ok(Some(AuspiceGenomeAnnotations {
    nuc,
    cdses,
    other: Value::default(),
  }))
}

pub(crate) fn sequence_auspice_node(
  name: &str,
  div: Option<f64>,
  branch_support: Option<f64>,
  mutations: Vec<Mutation>,
  command: &str,
) -> Result<AuspiceTreeNode, Report> {
  let node = auspice_node(
    name.to_owned(),
    finite_number(div, 6, command, name, "div")?,
    None,
    None,
    None,
    BTreeMap::new(),
    group_mutations(mutations)?,
    None,
  );
  with_branch_support(node, branch_support, command)
}

pub(crate) fn with_branch_support(
  mut node: AuspiceTreeNode,
  branch_support: Option<f64>,
  command: &str,
) -> Result<AuspiceTreeNode, Report> {
  let Some(branch_support) = branch_support else {
    return Ok(node);
  };
  ensure_finite(branch_support, command, &node.name, "input branch support")?;
  if let Value::Object(target) = &mut node.node_attrs.other {
    if target.contains_key("confidence") {
      return make_error!(
        "Node '{}' has a trait named 'confidence', which Auspice JSON also uses for the input branch support. \
         Rename the metadata column of the trait.",
        node.name
      );
    }
    target.insert("confidence".to_owned(), json!({ "value": branch_support }));
  }
  Ok(node)
}

pub(crate) fn auspice_node(
  name: String,
  div: Option<f64>,
  date: Option<f64>,
  date_confidence: Option<[f64; 2]>,
  bad_branch: Option<bool>,
  traits: BTreeMap<String, TraitValue>,
  mutations: BTreeMap<String, Vec<String>>,
  labels: Option<AuspiceTreeBranchAttrsLabels>,
) -> AuspiceTreeNode {
  AuspiceTreeNode {
    name,
    branch_attrs: AuspiceTreeBranchAttrs {
      mutations,
      labels,
      other: Value::default(),
    },
    node_attrs: AuspiceTreeNodeAttrs {
      div,
      num_date: date.map(|value| AuspiceNumDate {
        value,
        confidence: date_confidence,
      }),
      bad_branch: bad_branch.map(|bad| AuspiceTreeNodeAttr::new(if bad { "Yes" } else { "No" })),
      clade_membership: None,
      region: None,
      country: None,
      division: None,
      other: build_trait_attrs(traits),
    },
    children: vec![],
    other: Value::default(),
  }
}

pub(crate) fn auspice_from_graph<F>(graph: &Graph, data: AuspiceTreeData, mut convert: F) -> Result<AuspiceTree, Report>
where
  F: FnMut(&GraphNodeContext) -> Result<AuspiceTreeNode, Report>,
{
  let root_key = graph
    .get_exactly_one_root()
    .wrap_err("When converting graph to Auspice v2 JSON")?
    .key();
  let mut node_map = btreemap! {};
  let mut queue = VecDeque::from([(root_key, None)]);
  while let Some((node_key, edge_key)) = queue.pop_front() {
    let converted = convert(&GraphNodeContext { node_key, edge_key })?;
    if converted.node_attrs.div.is_none() && converted.node_attrs.num_date.is_none() {
      return make_error!(
        "Auspice v2 node '{}' requires divergence or numerical date data",
        converted.name
      );
    }
    node_map.insert(node_key, converted);
    let node = graph
      .get_node(node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} not found in graph"))?;
    for (child_key, child_edge_key) in graph.children_keys_of(node) {
      queue.push_back((child_key, Some(child_edge_key)));
    }
  }
  attach_auspice_children(graph, root_key, &mut node_map)?;
  let tree = node_map
    .remove(&root_key)
    .ok_or_else(|| make_internal_report!("Auspice root node {root_key} was not converted"))?;
  Ok(AuspiceTree { data, tree })
}

#[allow(
  clippy::field_scoped_visibility_modifiers,
  reason = "crate-internal fields are the record interface"
)]
pub(crate) struct GraphNodeContext {
  pub(crate) node_key: GraphNodeKey,
  pub(crate) edge_key: Option<GraphEdgeKey>,
}

fn attach_auspice_children(
  graph: &Graph,
  root_key: GraphNodeKey,
  node_map: &mut BTreeMap<GraphNodeKey, AuspiceTreeNode>,
) -> Result<(), Report> {
  let mut visited = btreeset! {};
  let mut stack = vec![root_key];
  while let Some(key) = stack.pop() {
    let node = graph
      .get_node(key)
      .ok_or_else(|| make_internal_report!("Node {key} not found in graph"))?;
    if visited.contains(&key) {
      let children = graph
        .children_keys_of(node)
        .map(|(child_key, _)| {
          node_map
            .remove(&child_key)
            .ok_or_else(|| make_internal_report!("Auspice child node {child_key} was not converted"))
        })
        .collect::<Result<Vec<_>, _>>()?;
      node_map
        .get_mut(&key)
        .ok_or_else(|| make_internal_report!("Auspice parent node {key} was not converted"))?
        .children = children;
    } else {
      visited.insert(key);
      stack.push(key);
      stack.extend(graph.children_keys_of(node).map(|(child_key, _)| child_key));
    }
  }
  Ok(())
}

#[derive(Debug)]
#[allow(
  clippy::field_scoped_visibility_modifiers,
  reason = "crate-internal fields are the record interface"
)]
pub(crate) struct MatOutput {
  pub(crate) tree: UsherTree,
  pub(crate) gaps: MatGapCounts,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[allow(
  clippy::field_scoped_visibility_modifiers,
  reason = "crate-internal fields are the record interface"
)]
pub(crate) struct MatGapCounts {
  pub(crate) deletions: usize,
  pub(crate) insertions: usize,
  pub(crate) substitutions: usize,
}

impl MatGapCounts {
  fn count(&mut self, mutations: &[Mutation], reference_gaps: &BTreeSet<usize>) {
    for mutation in mutations
      .iter()
      .filter(|mutation| mutation.track == MutationTrack::Nucleotide)
    {
      match &mutation.event {
        MutationEvent::Deletion(_) => self.deletions += 1,
        MutationEvent::Insertion(segment) => {
          if reference_gaps.range(segment.range.0..segment.range.1).next().is_some() {
            self.insertions += 1;
          }
        },
        MutationEvent::Substitution(substitution) => {
          if reference_gaps.contains(&substitution.pos()) {
            self.substitutions += 1;
          }
        },
      }
    }
  }

  pub(crate) fn warning(self) -> Option<String> {
    let mut clauses = vec![];
    if self.deletions > 0 {
      clauses.push(format!("wrote {} deletion(s) as missing data (N)", self.deletions));
    }
    if self.insertions > 0 || self.substitutions > 0 {
      clauses.push(format!(
        "left out {} insertion(s) and {} substitution(s) in alignment columns where the root sequence, which is the MAT reference, has a gap",
        self.insertions, self.substitutions
      ));
    }
    (!clauses.is_empty()).then(|| format!("UShER MAT has no gap state: {}", clauses.join("; ")))
  }
}

pub(crate) fn mutation_free_mat(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<MatOutput, Report> {
  mat_from_graph(graph, names, nwk_weights, None, |_node_key, _edge_key| Ok(vec![]))
}

pub(crate) fn mat_from_graph<F>(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  reference: Option<&str>,
  mut edge_mutations: F,
) -> Result<MatOutput, Report>
where
  F: FnMut(GraphNodeKey, GraphEdgeKey) -> Result<Vec<Mutation>, Report>,
{
  graph
    .get_exactly_one_root()
    .wrap_err("When converting graph to UShER MAT")?;
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;
  let reference_gaps: BTreeSet<usize> = reference.map_or_else(BTreeSet::new, |reference| {
    reference
      .bytes()
      .enumerate()
      .filter(|&(_, state)| state == u8::from(alphabet.gap()))
      .map(|(pos, _)| pos)
      .collect()
  });
  let mut missing_data = UnknownBridge::new(alphabet.unknown());
  let mut gaps = MatGapCounts::default();
  let mut node_mutations = vec![];
  let mut condensed_nodes = vec![];
  let mut metadata = vec![];
  graph.iter_depth_first_preorder_forward(|node| {
    let name = names[&node.key].clone().unwrap_or_default();
    let mutations = match node.parent_keys.as_slice() {
      [(parent_key, edge_key)] => {
        let mutations = edge_mutations(node.key, *edge_key)?;
        gaps.count(&mutations, &reference_gaps);
        let mutations = gaps_as_missing_data(mutations, alphabet.unknown())?;
        missing_data.bridge_edge(*parent_key, node.key, node.child_edge_keys.len(), mutations)?
      },
      _ => vec![],
    };
    let mut mutations = mutations
      .iter()
      .filter(|mutation| !is_in_reference_gap(mutation, &reference_gaps))
      .map(|mutation| mat_mutation(mutation, reference, &alphabet, &name))
      .collect::<Result<Vec<_>, _>>()?;
    mutations.sort_by_key(|mutation| mutation.position);
    node_mutations.push(UsherMutationList { mutation: mutations });
    condensed_nodes.push(UsherTreeNode {
      node_name: name,
      condensed_leaves: vec![],
    });
    metadata.push(UsherMetadata {
      clade_annotations: vec![],
    });
    Ok(())
  })?;
  Ok(MatOutput {
    tree: UsherTree {
      newick: nwk_write_str(graph, names, nwk_weights, &NwkWriteOptions::default())?,
      node_mutations,
      condensed_nodes,
      metadata,
    },
    gaps,
  })
}

fn gaps_as_missing_data(mutations: Vec<Mutation>, unknown: AsciiChar) -> Result<Vec<Mutation>, Report> {
  let mut converted = Vec::with_capacity(mutations.len());
  for mutation in mutations {
    let (segment, is_deletion) = match (&mutation.track, &mutation.event) {
      (MutationTrack::Nucleotide, MutationEvent::Deletion(segment)) => (segment, true),
      (MutationTrack::Nucleotide, MutationEvent::Insertion(segment)) => (segment, false),
      _ => {
        converted.push(mutation);
        continue;
      },
    };
    for (pos, &state) in (segment.range.0..segment.range.1).zip(segment.sequence.iter()) {
      if state == unknown {
        continue;
      }
      let substitution = if is_deletion {
        Sub::new(state, pos, unknown)?
      } else {
        Sub::new(unknown, pos, state)?
      };
      converted.push(Mutation::substitution(MutationTrack::Nucleotide, substitution));
    }
  }
  Ok(converted)
}

fn is_in_reference_gap(mutation: &Mutation, reference_gaps: &BTreeSet<usize>) -> bool {
  matches!(
    (&mutation.track, &mutation.event),
    (MutationTrack::Nucleotide, MutationEvent::Substitution(substitution)) if reference_gaps.contains(&substitution.pos())
  )
}

pub(crate) fn mat_mutation(
  mutation: &Mutation,
  reference: Option<&str>,
  alphabet: &Alphabet,
  node_name: &str,
) -> Result<UsherMutation, Report> {
  if mutation.track != MutationTrack::Nucleotide {
    return make_internal_error!(
      "Node '{node_name}' has an amino-acid mutation, but UShER MAT stores nucleotide mutations only"
    );
  }
  let MutationEvent::Substitution(substitution) = &mutation.event else {
    return make_internal_error!(
      "Node '{node_name}' has an insertion or deletion that was not converted to missing data"
    );
  };
  let reference = reference.ok_or_else(|| {
    eyre::eyre!("Node '{node_name}' has nucleotide mutations, but UShER MAT requires a root nucleotide reference")
  })?;
  let position = substitution
    .pos()
    .checked_add(1)
    .ok_or_else(|| eyre::eyre!("Node '{node_name}' mutation coordinate overflow"))?;
  let position = i32::try_from(position).wrap_err_with(|| {
    format!("Node '{node_name}' mutation position {position} exceeds the UShER MAT i32 coordinate range")
  })?;
  let reference_state = reference.as_bytes().get(substitution.pos()).copied().ok_or_else(|| {
    eyre::eyre!(
      "Node '{node_name}' mutation position {} is outside the root nucleotide reference of length {}",
      substitution.pos() + 1,
      reference.len()
    )
  })?;
  let reference_state = AsciiChar::try_new(reference_state)?;
  Ok(UsherMutation {
    position,
    ref_nuc: mat_nucleotide(reference_state, node_name, "root reference")?,
    par_nuc: mat_nucleotide(substitution.reff(), node_name, "parent")?,
    mut_nuc: mat_nucleotide_states(substitution.qry(), alphabet, node_name)?,
    chromosome: String::new(),
  })
}

fn mat_nucleotide_states(nucleotide: AsciiChar, alphabet: &Alphabet, node_name: &str) -> Result<Vec<i32>, Report> {
  let rejected = || {
    format!(
      "Node '{node_name}' has child nucleotide '{nucleotide}', but UShER MAT accepts only A, C, G, T, IUPAC ambiguity codes, or N"
    )
  };
  let states = alphabet.canonical_states(nucleotide).wrap_err_with(rejected)?;
  if states.is_empty() {
    return make_internal_error!(
      "Node '{node_name}' has child nucleotide '{nucleotide}', which stands for no nucleotide"
    );
  }
  states
    .iter()
    .map(|state| mat_nucleotide(state, node_name, "child"))
    .collect()
}

fn mat_nucleotide(nucleotide: AsciiChar, node_name: &str, role: &str) -> Result<i32, Report> {
  match char::from(nucleotide).to_ascii_uppercase() {
    'A' => Ok(0),
    'C' => Ok(1),
    'G' => Ok(2),
    'T' => Ok(3),
    state => {
      make_error!("Node '{node_name}' has {role} nucleotide '{state}', but UShER MAT accepts only A, C, G, or T")
    },
  }
}

pub(crate) fn group_mutations(mutations: Vec<Mutation>) -> Result<BTreeMap<String, Vec<String>>, Report> {
  let mut grouped = BTreeMap::new();
  for mutation in mutations {
    let track = match &mutation.track {
      MutationTrack::Nucleotide => NUC_TRACK.to_owned(),
      MutationTrack::AminoAcid(track) => {
        if track.is_empty()
          || !track
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'*' | b'_' | b'.' | b'(' | b')' | b'-'))
        {
          return make_error!("Auspice v2 cannot represent amino-acid mutation track '{track}'");
        }
        track.clone()
      },
    };
    if matches!(mutation.track, MutationTrack::Nucleotide) && !matches!(mutation.event, MutationEvent::Substitution(_))
    {
      continue;
    }
    grouped
      .entry(track)
      .or_insert_with(Vec::new)
      .extend(mutation_event_strings(&mutation.event)?);
  }
  Ok(grouped)
}

fn build_trait_attrs(traits: BTreeMap<String, TraitValue>) -> Value {
  Value::Object(
    traits
      .into_iter()
      .map(|(attribute, value)| {
        let mut fields = serde_json::Map::new();
        fields.insert("value".to_owned(), json!(value.value));
        if !value.confidence.is_empty() {
          fields.insert(
            "confidence".to_owned(),
            json!(
              value
                .confidence
                .iter()
                .map(|(state, probability)| (state.clone(), format_number(*probability, 3)))
                .collect::<BTreeMap<_, _>>()
            ),
          );
        }
        if let Some(entropy) = value.entropy {
          fields.insert("entropy".to_owned(), json!(format_number(entropy, 3)));
        }
        (attribute, Value::Object(fields))
      })
      .collect(),
  )
}

#[derive(Clone, Debug, Default, PartialEq)]
#[allow(
  clippy::field_scoped_visibility_modifiers,
  reason = "crate-internal fields are the record interface"
)]
pub(crate) struct TraitValue {
  pub(crate) value: String,
  pub(crate) confidence: BTreeMap<String, f64>,
  pub(crate) entropy: Option<f64>,
}

pub(crate) fn cumulative_branch_length_from(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  mut key: GraphNodeKey,
) -> Result<Option<f64>, Report> {
  let mut total = 0.0;
  while let Some((parent, edge_key)) = graph.node_parent(key)? {
    let Some(length) = branch_lengths[&edge_key] else {
      return Ok(None);
    };
    total += length;
    key = parent;
  }
  Ok(Some(total))
}

pub(crate) fn node_name_value(key: GraphNodeKey, name: Option<&str>) -> String {
  name.map_or_else(|| format!("node_{}", key.as_usize()), str::to_owned)
}

pub(crate) fn finite_number(
  value: Option<f64>,
  precision: i32,
  command: &str,
  node_name: &str,
  field: &str,
) -> Result<Option<f64>, Report> {
  ensure_optional_finite(value, command, node_name, field)?;
  Ok(value.map(|value| format_number(value, precision)))
}

fn ensure_optional_finite(value: Option<f64>, command: &str, node_name: &str, field: &str) -> Result<(), Report> {
  if let Some(value) = value {
    ensure_finite(value, command, node_name, field)?;
  }
  Ok(())
}

pub(crate) fn ensure_finite(value: f64, command: &str, node_name: &str, field: &str) -> Result<(), Report> {
  if value.is_finite() {
    Ok(())
  } else {
    make_error!("{command} node '{node_name}' has non-finite {field}={value}")
  }
}

#[allow(
  clippy::as_conversions,
  clippy::expect_used,
  reason = "count/index numeric cast is exact for the domain range; expect on a value an upstream invariant guarantees is present"
)]
pub(crate) fn format_number(number: f64, precision: i32) -> f64 {
  if number == 0.0 || !number.is_finite() {
    return number;
  }
  let integral = number.abs().trunc();
  let significand = if integral >= 1.0 {
    integral.log10().floor() as i32 + 1
  } else {
    0
  };
  let significant_figures = (significand + precision).max(1) as usize;
  format!("{number:.*e}", significant_figures - 1)
    .parse()
    .expect("a float formatted in scientific notation must parse back")
}

pub(crate) fn coloring(key: &str, title: &str, type_: &str) -> AuspiceColoring {
  AuspiceColoring {
    key: key.to_owned(),
    title: title.to_owned(),
    type_: type_.to_owned(),
    ..AuspiceColoring::default()
  }
}

pub(crate) fn generation_date() -> String {
  Utc::now().format("%Y-%m-%d").to_string()
}
