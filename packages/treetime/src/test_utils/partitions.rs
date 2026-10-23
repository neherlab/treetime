use crate::gtr::get_gtr::GtrModel;
use crate::partition::create::Representation;
use crate::partition::resolve::{PartitionEntry, resolve_partitions};
use crate::partition::set::{PartitionSet, build_partition_set};
use crate::partition::sites::SiteRange;
use crate::partition::split::{SequenceInput, SplitInput};
use crate::progress::NoopProgress;
use crate::seq::alignment::SequenceSource;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::str::FromStr;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_utils::make_internal_error;

pub(crate) fn codon_positions_input(
  alignments: Vec<SequenceSource>,
  model: &GtrModel,
) -> Result<SequenceInput, Report> {
  let [alignment] = alignments.as_slice() else {
    return make_internal_error!(
      "Codon positions split one alignment, but the input has {}",
      alignments.len()
    );
  };
  let length = alignment.length;
  if length % 3 != 0 {
    return make_internal_error!("Codon positions need a length divisible by 3, but the alignment has {length} sites");
  }
  let entries: Vec<PartitionEntry> = (1..=3)
    .map(|position| {
      Ok(PartitionEntry {
        name: None,
        alignment: Some(1),
        sites: Some(vec![SiteRange::from_str(&format!("{position}-{length}\\3"))?]),
        model: model.clone(),
      })
    })
    .try_collect::<_, _, Report>()?;
  let lengths = [length];
  let partitions = resolve_partitions(&entries, &lengths, model, Some(alignment.alphabet.name()))?.plans;
  Ok(SequenceInput { alignments, partitions })
}

pub(crate) fn codon_partition_set(
  representation: Representation,
  model: &GtrModel,
  graph: &Graph,
  sources: Vec<SequenceSource>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<PartitionSet, Report> {
  let SplitInput { partitions, .. } = codon_positions_input(sources, model)?.split()?;
  build_partition_set(representation, graph, partitions, branch_lengths, &NoopProgress)
}
