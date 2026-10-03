use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::attach::complete_alignment_for_leaves;
use crate::ancestral::plan::{ReconstructionOptions, ReconstructionPlan, reconstruct_partition};
use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::partition::create::Representation;
use crate::partition::marginal::sample::SampleMode;
use crate::progress::{LogSink, NoopProgress};
use crate::seq::alignment::{get_common_length, node_seq_inputs};
use crate::seq::mutation::{MutationEvent, MutationTrack, SequenceMutations, Sub};
use crate::seq::sink::SeqSink;
use crate::{make_error, make_internal_report};
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AlignmentRecord, AsciiChar, Seq};
use treetime_utils::sync::random::get_random_number_generator;
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

pub fn reconstruct_aa(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  params: &AaParams,
  cdses: Vec<CdsInput>,
  mut seq_sink: Option<&mut dyn SeqSink>,
  cancel: &dyn Cancel,
  log: &dyn LogSink,
) -> Result<AaNodeData, Report> {
  let mut rng = get_random_number_generator(params.seed);
  let options = ReconstructionOptions::new(params.impute_missing_data, params.sample_from_profile);
  let representation = Representation::resolve(params.dense);
  let branch_lengths = branch_lengths_or_zero(branch_lengths);
  let mut aa_node_data = AaNodeData::default();
  if let Some(sink) = seq_sink.as_deref_mut() {
    sink.on_topology(graph)?;
  }
  for (index, cds) in cdses.into_iter().enumerate() {
    let CdsInput {
      name,
      alphabet,
      gtr_model,
      sequences,
      annotation,
      reference_override,
    } = cds;
    let sequences = complete_alignment_for_leaves(graph, sequences, &alphabet, params.ignore_missing_alns, names, log)?;
    validate_cds_length(&name, annotation.as_ref(), get_common_length(&sequences)?)?;
    let unknown = alphabet.unknown();
    let node_inputs = node_seq_inputs(graph, names, sequences);
    let plan = ReconstructionPlan::Marginal {
      representation,
      model: gtr_model,
      gtr_refinement: None,
    };
    let partition = reconstruct_partition(
      graph,
      &plan,
      index,
      alphabet,
      &node_inputs,
      &branch_lengths,
      &options,
      &mut rng,
      cancel,
      &NoopProgress,
      log,
    )?;
    let mutations = partition
      .stream_sequences(
        graph,
        &MutationTrack::AminoAcid(name.clone()),
        params.include_leaves,
        seq_sink.as_deref_mut(),
      )
      .map_err(OperationError::into_report)?;
    let cds_data = collect_aa_cds_node_data(graph, mutations, unknown, &name, reference_override.as_ref())?;
    aa_node_data.add_cds(&name, cds_data);
  }

  Ok(aa_node_data)
}

pub struct CdsInput {
  pub name: String,
  pub alphabet: Alphabet,
  pub gtr_model: GtrModelName,
  pub sequences: Vec<AlignmentRecord>,
  pub annotation: Option<AugurNodeDataJsonAnnotationEntry>,
  pub reference_override: Option<Seq>,
}

pub struct AaParams {
  pub dense: Option<bool>,
  pub include_leaves: bool,
  pub impute_missing_data: bool,
  pub sample_from_profile: SampleMode,
  pub seed: Option<u64>,
  pub ignore_missing_alns: bool,
}

fn validate_cds_length(
  name: &str,
  annotation: Option<&AugurNodeDataJsonAnnotationEntry>,
  alignment_length: usize,
) -> Result<(), Report> {
  if let Some(cds_len) = annotation.and_then(annotation_cds_nuc_length) {
    let aa_len = i64::try_from(alignment_length)?;
    if 3 * aa_len != cds_len {
      return make_error!(
        "Translated alignment for CDS '{name}' has {aa_len} amino acids ({} nucleotides), which does not match \
         the annotated CDS length of {cds_len} nucleotides. Check that the annotation matches the translations.",
        3 * aa_len
      );
    }
  }
  Ok(())
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct AaNodeData {
  pub reference: BTreeMap<String, String>,
  pub node_aa_mutations: BTreeMap<GraphNodeKey, BTreeMap<String, Vec<MutationEvent>>>,
  pub root_aa_sequences: BTreeMap<String, String>,
}

impl AaNodeData {
  pub fn add_cds(&mut self, cds: &str, cds_data: AaCdsNodeData) {
    self.reference.insert(cds.to_owned(), cds_data.reference);
    self.root_aa_sequences.insert(cds.to_owned(), cds_data.root_sequence);
    for (node_key, mutations) in cds_data.node_mutations {
      self
        .node_aa_mutations
        .entry(node_key)
        .or_default()
        .insert(cds.to_owned(), mutations);
    }
  }
}

pub(crate) fn collect_aa_cds_node_data(
  graph: &Graph,
  mutations: SequenceMutations,
  unknown: AsciiChar,
  cds: &str,
  reference_override: Option<&Seq>,
) -> Result<AaCdsNodeData, Report> {
  let SequenceMutations {
    root_sequence: inferred_root,
    mut edge_mutations,
  } = mutations;
  let reference = reference_override.cloned().unwrap_or_else(|| inferred_root.clone());

  if reference.len() != inferred_root.len() {
    return make_error!(
      "AA root/reference sequence for CDS '{cds}' has length {}, but inferred root has length {}",
      reference.len(),
      inferred_root.len()
    );
  }

  let mut node_mutations = BTreeMap::new();
  for node in graph.get_nodes() {
    let node_key = node.key();
    let mutations = match graph.node_parent(node_key)? {
      None => diff_sequences(&reference, &inferred_root, unknown)?
        .into_iter()
        .map(MutationEvent::Substitution)
        .collect(),
      Some((_parent_key, edge_key)) => edge_mutations
        .remove(&edge_key)
        .ok_or_else(|| {
          make_internal_report!("No mutations were derived for the edge above node {node_key} in CDS '{cds}'")
        })?
        .into_iter()
        .map(|mutation| mutation.event)
        .collect(),
    };
    node_mutations.insert(node_key, mutations);
  }

  Ok(AaCdsNodeData {
    reference: reference.as_str().to_owned(),
    root_sequence: inferred_root.as_str().to_owned(),
    node_mutations,
  })
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct AaCdsNodeData {
  pub reference: String,
  pub root_sequence: String,
  pub node_mutations: BTreeMap<GraphNodeKey, Vec<MutationEvent>>,
}

pub(super) fn diff_sequences(reference: &Seq, query: &Seq, unknown: AsciiChar) -> Result<Vec<Sub>, Report> {
  if reference.len() != query.len() {
    return make_error!(
      "Cannot diff sequences with lengths {} and {}",
      reference.len(),
      query.len()
    );
  }

  reference
    .iter()
    .zip(query.iter())
    .enumerate()
    .filter(|(_pos, (reff, qry))| reff != qry && is_reportable_sub(**reff, **qry, unknown))
    .map(|(pos, (reff, qry))| Sub::new(*reff, pos, *qry))
    .collect()
}

fn is_reportable_sub(reff: AsciiChar, qry: AsciiChar, unknown: AsciiChar) -> bool {
  let gap = AsciiChar::from_byte_unchecked(b'-');
  reff != gap && qry != gap && reff != unknown && qry != unknown
}

pub(crate) fn annotation_cds_nuc_length(entry: &AugurNodeDataJsonAnnotationEntry) -> Option<i64> {
  if let Some(segments) = &entry.segments {
    Some(segments.iter().map(|segment| segment.end - segment.start + 1).sum())
  } else if let (Some(start), Some(end)) = (entry.start, entry.end) {
    Some(end - start + 1)
  } else {
    None
  }
}
