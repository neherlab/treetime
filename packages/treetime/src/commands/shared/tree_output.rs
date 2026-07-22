use crate::ancestral::pipeline::AncestralPartition;
use crate::clock::clock_graph::{EdgeClock, NodeClock};
use crate::commands::ancestral::result::AncestralGraphData;
use crate::commands::clock::run::ClockGraphData;
use crate::commands::mugration::augur_node_data::{build_confidence_map, compute_entropy};
use crate::commands::optimize::result::OptimizeGraphData;
use crate::commands::prune::result::PruneGraphData;
use crate::commands::timetree::result::TimetreeGraphData;
use crate::mugration::result::MugrationGraphData;
use crate::partition::augur::AugurNodeDataJsonAncestralPartition;
use crate::partition::traits::{BranchTopology, PartitionBranchOps};
use crate::payload::ancestral::{EdgeAncestral, NodeAncestral};
use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
use crate::seq::mutation::Sub;
use eyre::Report;
use log::warn;
use serde_json::{Value, json};
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_io::auspice::{AuspiceGraphContext, AuspiceWrite};
use treetime_io::auspice_types::{
  AuspiceColoring, AuspiceDisplayDefaults, AuspiceNumDate, AuspiceTreeBranchAttrs, AuspiceTreeData, AuspiceTreeMeta,
  AuspiceTreeNode, AuspiceTreeNodeAttr, AuspiceTreeNodeAttrs,
};
use treetime_io::phyloxml::{
  Phyloxml, PhyloxmlClade, PhyloxmlDate, PhyloxmlFromGraph, PhyloxmlGraphContext, PhyloxmlPhylogeny, PhyloxmlProperty,
};
use treetime_io::usher_mat::{
  UsherGraphContext, UsherMetadata, UsherMutation, UsherMutationList, UsherTreeNode, UsherWrite,
};
use treetime_primitives::AsciiChar;
use treetime_utils::make_error;

const APPLIES_BRANCH: &str = "parent_branch";
const APPLIES_NODE: &str = "node";
const BRANCH_LENGTH_UNIT: &str = "subs/site";
const COLORING_BAD_BRANCH: &str = "bad_branch";
const COLORING_GENOTYPE: &str = "gt";
const COLORING_NUM_DATE: &str = "num_date";
const DT_BOOLEAN: &str = "xsd:boolean";
const DT_DOUBLE: &str = "xsd:double";
const DT_STRING: &str = "xsd:string";
const NUC_GENE: &str = "nuc";
const REF_BAD_BRANCH: &str = "treetime:bad_branch";
const REF_DATE_INFERRED: &str = "treetime:date_inferred";
const REF_DIV: &str = "treetime:divergence";
const REF_GAMMA: &str = "treetime:gamma";
const REF_MUTATION: &str = "treetime:mutation";
const REF_TRAIT_PREFIX: &str = "treetime:trait:";

pub struct TreeOutputAdapter {
  has_bad_branch: bool,
  warned_ambiguous: bool,
  warned_non_nuc: bool,
}

pub trait TreeOutputNode: GraphNode + Named {
  fn output_divergence(&self) -> Option<f64> {
    None
  }

  fn output_date(&self) -> Option<f64> {
    None
  }

  fn output_bad_branch(&self) -> bool {
    false
  }
}

pub trait TreeOutputEdge: GraphEdge + HasBranchLength {
  fn output_gamma(&self) -> Option<f64> {
    None
  }
}

pub trait TreeOutputData<N, E>: Send + Sync
where
  N: TreeOutputNode,
  E: TreeOutputEdge,
{
  fn title(&self) -> Option<&str> {
    None
  }

  fn description(&self) -> Option<&str> {
    None
  }

  fn trait_attributes(&self) -> Vec<String> {
    vec![]
  }

  fn has_dates(&self) -> bool {
    false
  }

  fn has_bad_branch(&self) -> bool {
    false
  }

  fn has_mutations(&self) -> bool {
    false
  }

  fn root_sequences(&self, _graph: &Graph<N, E, Self>) -> Result<BTreeMap<String, String>, Report>
  where
    Self: Sized,
  {
    Ok(BTreeMap::new())
  }

  fn divergence(&self, _graph: &Graph<N, E, Self>, _key: GraphNodeKey, node: &N) -> Result<Option<f64>, Report>
  where
    Self: Sized,
  {
    Ok(node.output_divergence())
  }

  fn date_confidence(&self, _key: GraphNodeKey) -> Option<[f64; 2]> {
    None
  }

  fn date_is_inferred(&self, _key: GraphNodeKey) -> bool {
    false
  }

  fn traits(&self, _key: GraphNodeKey) -> BTreeMap<String, TreeOutputTrait> {
    BTreeMap::new()
  }

  fn mutations(&self, _graph: &Graph<N, E, Self>, _key: GraphEdgeKey) -> Result<Vec<TreeOutputMutation>, Report>
  where
    Self: Sized,
  {
    Ok(vec![])
  }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TreeOutputMutation {
  pub gene: String,
  pub position: usize,
  pub parent: AsciiChar,
  pub child: AsciiChar,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct TreeOutputTrait {
  pub value: String,
  pub confidence: BTreeMap<String, f64>,
  pub entropy: Option<f64>,
}

impl TreeOutputNode for NodeAncestral {}

impl TreeOutputNode for NodeClock {
  fn output_divergence(&self) -> Option<f64> {
    Some(self.div)
  }

  fn output_date(&self) -> Option<f64> {
    self.time
  }

  fn output_bad_branch(&self) -> bool {
    self.bad_branch || self.is_outlier
  }
}

impl TreeOutputNode for NodeTimetree {
  fn output_divergence(&self) -> Option<f64> {
    Some(self.div)
  }

  fn output_date(&self) -> Option<f64> {
    self.time
  }

  fn output_bad_branch(&self) -> bool {
    self.bad_branch || self.is_outlier
  }
}

impl TreeOutputEdge for EdgeAncestral {}
impl TreeOutputEdge for EdgeClock {}

impl TreeOutputEdge for EdgeTimetree {
  fn output_gamma(&self) -> Option<f64> {
    Some(self.gamma)
  }
}

impl<N, E> TreeOutputData<N, E> for ()
where
  N: TreeOutputNode,
  E: TreeOutputEdge,
{
}

impl TreeOutputData<NodeAncestral, EdgeAncestral> for AncestralGraphData {
  fn has_mutations(&self) -> bool {
    self.partition.is_some()
  }

  fn root_sequences(
    &self,
    graph: &Graph<NodeAncestral, EdgeAncestral, Self>,
  ) -> Result<BTreeMap<String, String>, Report> {
    self.partition.as_ref().map_or_else(
      || Ok(BTreeMap::new()),
      |partition| {
        Ok(BTreeMap::from([(
          NUC_GENE.to_owned(),
          ancestral_root_sequence(partition, graph)?,
        )]))
      },
    )
  }

  fn mutations(
    &self,
    graph: &Graph<NodeAncestral, EdgeAncestral, Self>,
    key: GraphEdgeKey,
  ) -> Result<Vec<TreeOutputMutation>, Report> {
    self
      .partition
      .as_ref()
      .map_or_else(|| Ok(vec![]), |partition| ancestral_mutations(partition, graph, key))
  }
}

impl TreeOutputData<NodeAncestral, EdgeAncestral> for OptimizeGraphData {
  fn has_mutations(&self) -> bool {
    !self.dense_partitions.is_empty() || !self.sparse_partitions.is_empty()
  }

  fn root_sequences(
    &self,
    graph: &Graph<NodeAncestral, EdgeAncestral, Self>,
  ) -> Result<BTreeMap<String, String>, Report> {
    let sequence = if let Some(partition) = self.dense_partitions.first() {
      Some(partition.read_arc().root_sequence(graph)?)
    } else if let Some(partition) = self.sparse_partitions.first() {
      Some(partition.read_arc().root_sequence(graph)?)
    } else {
      None
    };
    Ok(
      sequence
        .map(|sequence| BTreeMap::from([(NUC_GENE.to_owned(), sequence.to_string())]))
        .unwrap_or_default(),
    )
  }

  fn mutations(
    &self,
    graph: &Graph<NodeAncestral, EdgeAncestral, Self>,
    key: GraphEdgeKey,
  ) -> Result<Vec<TreeOutputMutation>, Report> {
    if let Some(partition) = self.dense_partitions.first() {
      subs_to_output(PartitionBranchOps::edge_subs(&*partition.read_arc(), graph, key)?)
    } else if let Some(partition) = self.sparse_partitions.first() {
      subs_to_output(PartitionBranchOps::edge_subs(&*partition.read_arc(), graph, key)?)
    } else {
      Ok(vec![])
    }
  }
}

impl TreeOutputData<NodeAncestral, EdgeAncestral> for PruneGraphData {
  fn has_mutations(&self) -> bool {
    !self.partitions.is_empty()
  }

  fn mutations(
    &self,
    graph: &Graph<NodeAncestral, EdgeAncestral, Self>,
    key: GraphEdgeKey,
  ) -> Result<Vec<TreeOutputMutation>, Report> {
    self.partitions.first().map_or_else(
      || Ok(vec![]),
      |partition| subs_to_output(PartitionBranchOps::edge_subs(&*partition.read_arc(), graph, key)?),
    )
  }
}

impl TreeOutputData<NodeClock, EdgeClock> for ClockGraphData {
  fn title(&self) -> Option<&str> {
    Some("TreeTime clock analysis")
  }

  fn has_dates(&self) -> bool {
    true
  }

  fn has_bad_branch(&self) -> bool {
    true
  }
}

impl TreeOutputData<NodeAncestral, EdgeAncestral> for MugrationGraphData {
  fn trait_attributes(&self) -> Vec<String> {
    vec![self.traits.attribute.clone()]
  }

  fn traits(&self, key: GraphNodeKey) -> BTreeMap<String, TreeOutputTrait> {
    let Some(value) = self.partition.get_reconstructed_trait(key) else {
      return BTreeMap::new();
    };
    let profile = self.partition.get_confidence(key);
    let confidence = profile
      .as_ref()
      .map(|profile| build_confidence_map(&self.partition.states, profile))
      .unwrap_or_default();
    let entropy = profile.as_ref().map(compute_entropy);
    BTreeMap::from([(
      self.traits.attribute.clone(),
      TreeOutputTrait {
        value,
        confidence,
        entropy,
      },
    )])
  }
}

impl TreeOutputData<NodeTimetree, EdgeTimetree> for TimetreeGraphData {
  fn title(&self) -> Option<&str> {
    Some("TreeTime timetree analysis")
  }

  fn has_dates(&self) -> bool {
    true
  }

  fn has_bad_branch(&self) -> bool {
    true
  }

  fn has_mutations(&self) -> bool {
    !self.partitions.is_empty()
  }

  fn divergence(
    &self,
    graph: &Graph<NodeTimetree, EdgeTimetree, Self>,
    key: GraphNodeKey,
    node: &NodeTimetree,
  ) -> Result<Option<f64>, Report> {
    self.mutation_counts.as_ref().map_or(Ok(Some(node.div)), |counts| {
      cumulative_mutation_count(graph, key, counts).map(Some)
    })
  }

  fn date_confidence(&self, key: GraphNodeKey) -> Option<[f64; 2]> {
    self
      .confidence_intervals
      .as_ref()?
      .iter()
      .find_map(|interval| (interval.key == key).then_some([interval.lower, interval.upper]))
  }

  fn mutations(
    &self,
    graph: &Graph<NodeTimetree, EdgeTimetree, Self>,
    key: GraphEdgeKey,
  ) -> Result<Vec<TreeOutputMutation>, Report> {
    self.partitions.first().map_or_else(
      || Ok(vec![]),
      |partition| subs_to_output(partition.read_arc().edge_subs(graph, key)?),
    )
  }
}

impl<N, E, D> AuspiceWrite<N, E, D> for TreeOutputAdapter
where
  N: TreeOutputNode,
  E: TreeOutputEdge,
  D: TreeOutputData<N, E>,
{
  fn new(graph: &Graph<N, E, D>) -> Result<Self, Report> {
    Ok(Self {
      has_bad_branch: graph.data().has_bad_branch(),
      warned_ambiguous: false,
      warned_non_nuc: false,
    })
  }

  fn auspice_data_from_graph_data(&self, graph: &Graph<N, E, D>) -> Result<AuspiceTreeData, Report> {
    let data = graph.data();
    let trait_attributes = data.trait_attributes();
    let mut colorings = vec![];
    if data.has_dates() {
      colorings.push(coloring(COLORING_NUM_DATE, "Date", "continuous"));
    }
    if data.has_bad_branch() {
      colorings.push(coloring(COLORING_BAD_BRANCH, "Excluded", "categorical"));
    }
    for attribute in &trait_attributes {
      colorings.push(coloring(attribute, attribute, "categorical"));
    }
    if data.has_mutations() {
      colorings.push(coloring(COLORING_GENOTYPE, "Genotype", "categorical"));
    }

    let color_by = trait_attributes
      .first()
      .cloned()
      .or_else(|| data.has_bad_branch().then(|| COLORING_BAD_BRANCH.to_owned()))
      .or_else(|| data.has_dates().then(|| COLORING_NUM_DATE.to_owned()));
    let mut filters = trait_attributes;
    if data.has_bad_branch() {
      filters.push(COLORING_BAD_BRANCH.to_owned());
    }
    let root_sequences = data.root_sequences(graph)?;

    Ok(AuspiceTreeData {
      version: Some("v2".to_owned()),
      meta: AuspiceTreeMeta {
        title: data
          .title()
          .map(str::to_owned)
          .or_else(|| Some("TreeTime analysis".to_owned())),
        description: data.description().map(str::to_owned),
        panels: vec!["tree".to_owned()],
        colorings,
        filters,
        display_defaults: AuspiceDisplayDefaults {
          color_by,
          ..AuspiceDisplayDefaults::default()
        },
        ..AuspiceTreeMeta::default()
      },
      root_sequence: (!root_sequences.is_empty()).then_some(root_sequences),
      other: Value::default(),
    })
  }

  fn auspice_node_from_graph_components(
    &mut self,
    context: &AuspiceGraphContext<N, E, D>,
  ) -> Result<AuspiceTreeNode, Report> {
    let name = context.node.name().map_or_else(
      || format!("node_{}", context.node_key.as_usize()),
      |name| name.as_ref().to_owned(),
    );
    let div = context
      .graph
      .data()
      .divergence(context.graph, context.node_key, context.node)?;
    let div = finite_number(div, 6, &name, "div")?;
    let date = finite_number(context.node.output_date(), 3, &name, "date")?;
    let confidence = context.graph.data().date_confidence(context.node_key);
    let confidence = match confidence {
      Some([lower, upper]) if lower.is_finite() && upper.is_finite() => {
        Some([format_number(lower, 3), format_number(upper, 3)])
      },
      Some([lower, upper]) => return make_error!("Node '{name}' has non-finite date confidence [{lower}, {upper}]"),
      None => None,
    };
    let num_date = date.map(|value| AuspiceNumDate { value, confidence });
    let bad_branch = self
      .has_bad_branch
      .then(|| AuspiceTreeNodeAttr::new(if context.node.output_bad_branch() { "Yes" } else { "No" }));
    let traits = context.graph.data().traits(context.node_key);
    let mutations = match context.edge_key {
      Some(edge_key) => group_mutations(&context.graph.data().mutations(context.graph, edge_key)?),
      None => BTreeMap::new(),
    };

    Ok(AuspiceTreeNode {
      name,
      branch_attrs: AuspiceTreeBranchAttrs {
        mutations,
        labels: None,
        other: Value::default(),
      },
      node_attrs: AuspiceTreeNodeAttrs {
        div,
        num_date,
        bad_branch,
        clade_membership: None,
        region: None,
        country: None,
        division: None,
        other: build_trait_attrs(&traits),
      },
      children: vec![],
      other: Value::default(),
    })
  }
}

impl<N, E, D> PhyloxmlFromGraph<N, E, D> for TreeOutputAdapter
where
  N: TreeOutputNode,
  E: TreeOutputEdge,
  D: TreeOutputData<N, E>,
{
  fn phyloxml_data_from_graph_data(graph: &Graph<N, E, D>) -> Result<Phyloxml, Report> {
    let data = graph.data();
    Ok(Phyloxml {
      phylogeny: vec![PhyloxmlPhylogeny {
        rooted: true,
        rerootable: None,
        branch_length_unit: Some(BRANCH_LENGTH_UNIT.to_owned()),
        phylogeny_type: None,
        name: data.title().map(str::to_owned),
        id: None,
        description: data.description().map(str::to_owned),
        date: None,
        confidence: vec![],
        clade: None,
        clade_relation: vec![],
        sequence_relation: vec![],
        property: vec![],
        other: BTreeMap::new(),
      }],
      other: BTreeMap::new(),
    })
  }

  fn phyloxml_node_from_graph_components(context: &PhyloxmlGraphContext<N, E, D>) -> Result<PhyloxmlClade, Report> {
    let mut property = vec![];
    if let Some(div) = context
      .graph
      .data()
      .divergence(context.graph, context.node_key, context.node)?
    {
      property.push(make_property(REF_DIV, DT_DOUBLE, APPLIES_NODE, &div.to_string()));
    }
    if context.node.output_bad_branch() {
      property.push(make_property(REF_BAD_BRANCH, DT_BOOLEAN, APPLIES_NODE, "true"));
    }
    if context.graph.data().date_is_inferred(context.node_key) {
      property.push(make_property(REF_DATE_INFERRED, DT_BOOLEAN, APPLIES_NODE, "true"));
    }
    for (attribute, value) in context.graph.data().traits(context.node_key) {
      property.push(make_property(
        &format!("{REF_TRAIT_PREFIX}{attribute}"),
        DT_STRING,
        APPLIES_NODE,
        &value.value,
      ));
    }
    if let Some(edge) = context.edge {
      if let Some(gamma) = edge.output_gamma() {
        property.push(make_property(REF_GAMMA, DT_DOUBLE, APPLIES_BRANCH, &gamma.to_string()));
      }
    }
    if let Some(edge_key) = context.edge_key {
      for mutation in context.graph.data().mutations(context.graph, edge_key)? {
        property.push(make_property(
          REF_MUTATION,
          DT_STRING,
          APPLIES_BRANCH,
          &format!("{}:{}", mutation.gene, mutation_string(&mutation)),
        ));
      }
    }
    let date = context.node.output_date().map(|value| {
      let (minimum, maximum) = context
        .graph
        .data()
        .date_confidence(context.node_key)
        .map_or((None, None), |[lower, upper]| (Some(lower), Some(upper)));
      PhyloxmlDate {
        desc: None,
        value: Some(value),
        minimum,
        maximum,
        unit: Some("year".to_owned()),
      }
    });

    Ok(PhyloxmlClade {
      name: context.node.name().map(|name| name.as_ref().to_owned()),
      branch_length_elem: context.edge.and_then(HasBranchLength::branch_length),
      branch_length_attr: None,
      confidence: vec![],
      width: None,
      color: None,
      node_id: None,
      taxonomy: vec![],
      sequence: vec![],
      events: None,
      binary_characters: None,
      distribution: vec![],
      date,
      reference: vec![],
      property,
      clade: vec![],
      other: BTreeMap::new(),
    })
  }
}

impl<N, E, D> UsherWrite<N, E, D> for TreeOutputAdapter
where
  N: TreeOutputNode + treetime_io::nwk::NodeToNwk,
  E: TreeOutputEdge + treetime_io::nwk::EdgeToNwk,
  D: TreeOutputData<N, E>,
{
  fn new(graph: &Graph<N, E, D>) -> Result<Self, Report> {
    Ok(Self {
      has_bad_branch: graph.data().has_bad_branch(),
      warned_ambiguous: false,
      warned_non_nuc: false,
    })
  }

  fn usher_node_from_graph_components(
    &mut self,
    context: &UsherGraphContext<N, E, D>,
  ) -> Result<(UsherTreeNode, UsherMutationList, UsherMetadata), Report> {
    let reference = context
      .graph
      .data()
      .root_sequences(context.graph)?
      .remove(NUC_GENE)
      .unwrap_or_default();
    let mutations = match context.edge_key {
      Some(edge_key) => context.graph.data().mutations(context.graph, edge_key)?,
      None => vec![],
    };
    let mut encoded = vec![];
    for mutation in mutations {
      if mutation.gene != NUC_GENE {
        if !self.warned_non_nuc {
          warn!("UShER MAT format carries only nucleotide substitutions; dropping amino-acid mutations");
          self.warned_non_nuc = true;
        }
        continue;
      }
      let (Ok(par_nuc), Ok(mut_nuc)) = (nuc_to_int(mutation.parent), nuc_to_int(mutation.child)) else {
        if !self.warned_ambiguous {
          warn!("UShER MAT format cannot encode non-ACGT nucleotides; skipping ambiguous substitutions");
          self.warned_ambiguous = true;
        }
        continue;
      };
      let ref_nuc = reference
        .as_bytes()
        .get(mutation.position)
        .copied()
        .and_then(|nucleotide| AsciiChar::try_new(nucleotide).ok())
        .and_then(|nucleotide| nuc_to_int(nucleotide).ok())
        .unwrap_or(par_nuc);
      encoded.push(UsherMutation {
        position: i32::try_from(mutation.position + 1)?,
        ref_nuc,
        par_nuc,
        mut_nuc: vec![mut_nuc],
        chromosome: String::new(),
      });
    }

    Ok((
      UsherTreeNode {
        node_name: context
          .node
          .name()
          .map(|name| name.as_ref().to_owned())
          .unwrap_or_default(),
        condensed_leaves: vec![],
      },
      UsherMutationList { mutation: encoded },
      UsherMetadata {
        clade_annotations: vec![],
      },
    ))
  }
}

fn ancestral_root_sequence(partition: &AncestralPartition, graph: &dyn BranchTopology) -> Result<String, Report> {
  let sequence = match partition {
    AncestralPartition::Fitch(partition) => partition.read_arc().root_sequence(graph)?,
    AncestralPartition::Sparse(partition) => partition.read_arc().root_sequence(graph)?,
    AncestralPartition::Dense(partition) => partition.read_arc().root_sequence(graph)?,
  };
  Ok(sequence.to_string())
}

fn ancestral_mutations(
  partition: &AncestralPartition,
  graph: &dyn BranchTopology,
  key: GraphEdgeKey,
) -> Result<Vec<TreeOutputMutation>, Report> {
  let substitutions = match partition {
    AncestralPartition::Fitch(partition) => partition.read_arc().edge_subs(graph, key)?,
    AncestralPartition::Sparse(partition) => PartitionBranchOps::edge_subs(&*partition.read_arc(), graph, key)?,
    AncestralPartition::Dense(partition) => PartitionBranchOps::edge_subs(&*partition.read_arc(), graph, key)?,
  };
  subs_to_output(substitutions)
}

fn subs_to_output(substitutions: Vec<Sub>) -> Result<Vec<TreeOutputMutation>, Report> {
  Ok(
    substitutions
      .into_iter()
      .map(|substitution| TreeOutputMutation {
        gene: NUC_GENE.to_owned(),
        position: substitution.pos(),
        parent: substitution.reff(),
        child: substitution.qry(),
      })
      .collect(),
  )
}

fn cumulative_mutation_count<N, E, D>(
  graph: &Graph<N, E, D>,
  mut key: GraphNodeKey,
  counts: &BTreeMap<GraphEdgeKey, usize>,
) -> Result<f64, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let mut count = 0;
  while let Some((parent, edge)) = graph.node_parent(key)? {
    count += counts.get(&edge).copied().unwrap_or_default();
    key = parent;
  }
  Ok(count as f64)
}

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

fn finite_number(value: Option<f64>, precision: i32, node_name: &str, field: &str) -> Result<Option<f64>, Report> {
  match value {
    Some(value) if value.is_finite() => Ok(Some(format_number(value, precision))),
    Some(value) => make_error!("Node '{node_name}' has non-finite {field}={value}"),
    None => Ok(None),
  }
}

fn coloring(key: &str, title: &str, type_: &str) -> AuspiceColoring {
  AuspiceColoring {
    key: key.to_owned(),
    title: title.to_owned(),
    type_: type_.to_owned(),
    ..AuspiceColoring::default()
  }
}

fn build_trait_attrs(traits: &BTreeMap<String, TreeOutputTrait>) -> Value {
  let map = traits
    .iter()
    .map(|(attribute, value)| {
      let mut fields = serde_json::Map::new();
      fields.insert("value".to_owned(), json!(value.value));
      if !value.confidence.is_empty() {
        let confidence = value
          .confidence
          .iter()
          .map(|(state, probability)| (state.clone(), format_number(*probability, 3)))
          .collect::<BTreeMap<_, _>>();
        fields.insert("confidence".to_owned(), json!(confidence));
      }
      if let Some(entropy) = value.entropy {
        fields.insert("entropy".to_owned(), json!(format_number(entropy, 3)));
      }
      (attribute.clone(), Value::Object(fields))
    })
    .collect();
  Value::Object(map)
}

fn group_mutations(mutations: &[TreeOutputMutation]) -> BTreeMap<String, Vec<String>> {
  let mut grouped: BTreeMap<String, Vec<String>> = BTreeMap::new();
  for mutation in mutations {
    grouped
      .entry(mutation.gene.clone())
      .or_default()
      .push(mutation_string(mutation));
  }
  grouped
}

fn mutation_string(mutation: &TreeOutputMutation) -> String {
  format!("{}{}{}", mutation.parent, mutation.position + 1, mutation.child)
}

fn make_property(ref_: &str, datatype: &str, applies_to: &str, value: &str) -> PhyloxmlProperty {
  PhyloxmlProperty {
    value: value.to_owned(),
    ref_: ref_.to_owned(),
    unit: None,
    datatype: datatype.to_owned(),
    applies_to: applies_to.to_owned(),
    id_ref: None,
  }
}

fn nuc_to_int(nucleotide: AsciiChar) -> Result<i32, Report> {
  match char::from(nucleotide).to_ascii_uppercase() {
    'A' => Ok(0),
    'C' => Ok(1),
    'G' => Ok(2),
    'T' => Ok(3),
    other => make_error!("UShER MAT cannot encode nucleotide '{other}': expected A, C, G, or T"),
  }
}
