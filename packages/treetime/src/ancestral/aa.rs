#[cfg(test)]
mod __tests__;

use crate::ancestral::multi::{MarginalPartitionParams, PartitionPlan, reconstruct_marginal_partition};
use crate::ancestral::pipeline::AncestralPartition;
use crate::make_error;
use crate::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub};
use crate::seq::sink::{SeqItem, SeqSink, SeqTrack};
use eyre::Report;
use itertools::Itertools;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AsciiChar, Seq};
use treetime_utils::sync::random::get_random_number_generator;
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

pub fn reconstruct_aa(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  params: &MarginalPartitionParams,
  plans: Vec<PartitionPlan>,
  mut seq_sink: Option<Box<dyn SeqSink>>,
) -> Result<AaNodeData, Report> {
  let mut rng = get_random_number_generator(params.seed);
  let mut aa_node_data = AaNodeData::default();
  if let Some(sink) = seq_sink.as_mut() {
    sink.on_topology(graph)?;
  }
  for (index, plan) in plans.into_iter().enumerate() {
    let reconstructed = reconstruct_marginal_partition(graph, index, plan, params, names, branch_lengths, &mut rng)?;
    let guard = &reconstructed.partition;

    if let Some(annotation) = &reconstructed.annotation
      && let Some(cds_len) = annotation_cds_nuc_length(annotation)
    {
      let aa_len = i64::try_from(guard.sequence_length())?;
      if 3 * aa_len != cds_len {
        return make_error!(
          "Translated alignment for CDS '{}' has {aa_len} amino acids ({} nucleotides), which does not match \
           the annotated CDS length of {cds_len} nucleotides. Check that the annotation matches the translations.",
          reconstructed.name,
          3 * aa_len
        );
      }
    }

    let cds_data = collect_aa_cds_node_data(
      graph,
      guard,
      &reconstructed.name,
      names,
      reconstructed.reference_override.as_ref(),
    )?;
    aa_node_data.add_cds(&reconstructed.name, cds_data);

    if let Some(sink) = seq_sink.as_mut() {
      for node in graph.get_nodes() {
        let node_key = node.key();
        let seq = guard.augur_node_sequence(node_key);
        sink.emit(SeqItem {
          key: node_key,
          track: SeqTrack::Aa(reconstructed.name.as_str()),
          seq: &seq,
        })?;
      }
    }
  }

  Ok(aa_node_data)
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

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct AaCdsNodeData {
  pub reference: String,
  pub root_sequence: String,
  pub node_mutations: BTreeMap<GraphNodeKey, Vec<MutationEvent>>,
}

pub(crate) fn collect_aa_cds_node_data(
  graph: &Graph,
  partition: &AncestralPartition,
  cds: &str,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  reference_override: Option<&Seq>,
) -> Result<AaCdsNodeData, Report> {
  let root_key = graph.root_key()?;
  let inferred_root = partition.augur_root_sequence(graph)?;
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
    let node_guard = node;
    let node_key = node_guard.key();
    let node_name = names[&node_key]
      .as_deref()
      .map_or_else(|| format!("node_{}", node_key.0), str::to_owned);

    let mutations = if node_key == root_key {
      diff_sequences(&reference, &inferred_root, partition.ambiguous_char())?
        .into_iter()
        .map(MutationEvent::Substitution)
        .collect()
    } else {
      let (_parent_key, edge_key) = graph
        .node_parent(node_key)?
        .ok_or_else(|| eyre::eyre!("Non-root node '{node_name}' has no parent while collecting AA node data"))?;
      let substitutions = partition
        .edge_subs(graph, edge_key)?
        .into_iter()
        .sorted_by_key(Sub::pos)
        .map(MutationEvent::Substitution)
        .map(Ok);
      let indels = partition
        .edge_indels(edge_key)
        .into_iter()
        .map(|indel| Mutation::indel(MutationTrack::AminoAcid(cds.to_owned()), &indel).map(|mutation| mutation.event));
      substitutions.chain(indels).collect::<Result<Vec<_>, Report>>()?
    };

    node_mutations.insert(node_key, mutations);
  }

  Ok(AaCdsNodeData {
    reference: reference.as_str().to_owned(),
    root_sequence: inferred_root.as_str().to_owned(),
    node_mutations,
  })
}

fn diff_sequences(reference: &Seq, query: &Seq, unknown: AsciiChar) -> Result<Vec<Sub>, Report> {
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
