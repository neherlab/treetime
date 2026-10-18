use crate::branch_lengths::branch_lengths_or_zero;
use crate::optimize::indel::EdgeIndels;
use crate::partition::set::PartitionSet;
use crate::progress::LogSink;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_utils::cast::ToF64Exact;

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub enum PartitionRates {
  #[default]
  Estimate,
  Equal,
}

impl PartitionRates {
  fn estimates(self, partitions: &PartitionSet) -> bool {
    self == Self::Estimate && partitions.len() > 1
  }

  fn fixes_equal(self, partitions: &PartitionSet) -> bool {
    self == Self::Equal && partitions.len() > 1
  }
}

pub(crate) fn update_partition_rates(
  graph: &Graph,
  partitions: PartitionSet,
  mode: PartitionRates,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  indels: Option<&mut EdgeIndels>,
  log: &dyn LogSink,
) -> Result<PartitionSet, Report> {
  if !mode.estimates(&partitions) {
    return Ok(partitions);
  }
  let mut partitions = partitions.estimate_rates(graph, &branch_lengths_or_zero(branch_lengths), log)?;
  normalize_partition_rates(&mut partitions, branch_lengths, indels)?;
  debug!(
    "Relative rates of the alignments: {}",
    partitions.rates().iter().map(|rate| format!("{rate:.6}")).join(", ")
  );
  Ok(partitions)
}

pub(crate) fn fix_equal_rates(
  graph: &Graph,
  partitions: PartitionSet,
  mode: PartitionRates,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<PartitionSet, Report> {
  let equal_rates = vec![1.0; partitions.len()];
  if !mode.fixes_equal(&partitions) || partitions.rates() == equal_rates {
    return Ok(partitions);
  }
  let mut partitions = partitions;
  partitions.set_rates(&equal_rates)?;
  Ok(
    partitions
      .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?
      .0,
  )
}

pub(crate) fn normalize_partition_rates(
  partitions: &mut PartitionSet,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  indels: Option<&mut EdgeIndels>,
) -> Result<(), Report> {
  let mut partitions = partitions.lengths_and_gtrs_mut();
  let total_length: usize = partitions.iter().map(|(length, _)| *length).sum();

  if total_length == 0 {
    return Ok(());
  }

  let weighted_rate: f64 = partitions
    .iter()
    .map(|(length, gtr)| length.to_f64_exact() * gtr.mu())
    .sum();

  let total_average = weighted_rate / total_length.to_f64_exact();

  if total_average == 0.0 {
    return Ok(());
  }

  for (_, gtr) in &mut partitions {
    gtr.set_mu(gtr.mu() / total_average)?;
  }

  for value in branch_lengths.values_mut().flatten() {
    *value *= total_average;
  }

  if let Some(indels) = indels {
    indels.divide_rates(total_average);
  }
  Ok(())
}
