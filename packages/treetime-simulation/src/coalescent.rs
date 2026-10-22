use crate::tc_schedule::TcSchedule;
use eyre::Report;
use itertools::Itertools;
use rand::{Rng, RngExt};
use rand_distr::{Distribution, Exp1};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_internal_error};

pub const PRESENT: f64 = 2025.0;

pub fn simulate_coalescent(
  schedule: &TcSchedule,
  sampling: &[SampleBatch],
  rng: &mut impl Rng,
) -> Result<SimulatedTree, Report> {
  let batches = sorted_batches(sampling)?;
  let mut pending = batches.iter().peekable();
  let mut graph = Graph::new();
  let mut tau = BTreeMap::new();
  let mut tips = vec![];
  let mut mergers = vec![];
  let mut active: Vec<GraphNodeKey> = vec![];
  let mut time = 0.0;

  loop {
    if active.len() < 2 {
      let Some(batch) = pending.next() else {
        break;
      };
      time = batch.tau;
      add_batch(batch, &mut graph, &mut tau, &mut tips, &mut active);
      continue;
    }

    let pair_count = active.len().to_f64_exact() * (active.len() - 1).to_f64_exact() * 0.5;
    let exp1: f64 = Exp1.sample(rng);
    let merger_time = schedule.time_after_hazard(time, exp1 / pair_count)?;

    if let Some(batch) = pending.next_if(|batch| batch.tau <= merger_time) {
      time = batch.tau;
      add_batch(batch, &mut graph, &mut tau, &mut tips, &mut active);
      continue;
    }

    if !merger_time.is_finite() {
      return make_internal_error!(
        "Coalescent merger time is not finite after time {time} with {} lineages",
        active.len()
      );
    }
    if merger_time <= time {
      return make_error!(
        "Coalescent merger time {merger_time} does not exceed the current time {time} under the Tc schedule {schedule:?}: the wait of the next merger is below the rounding of the time"
      );
    }

    let first = rng.random_range(0..active.len());
    let second = rng.random_range(0..active.len() - 1);
    let second = if second >= first { second + 1 } else { second };
    let parent = graph.add_node();
    for child in [active[first], active[second]] {
      graph.add_edge(parent, child)?;
    }
    active.swap_remove(first.max(second));
    active.swap_remove(first.min(second));
    active.push(parent);
    tau.insert(parent, merger_time);
    mergers.push(parent);
    time = merger_time;
  }

  graph.build()?;
  let names = NodeNames::from_input(node_names(&tips, &mergers), &graph)?.names;
  Ok(SimulatedTree { graph, names, tau })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SampleBatch {
  pub tau: f64,
  pub count: usize,
}

#[derive(Debug)]
pub struct SimulatedTree {
  pub graph: Graph,
  pub names: NodeNames,
  pub tau: BTreeMap<GraphNodeKey, f64>,
}

impl SimulatedTree {
  pub fn date(&self, node: GraphNodeKey) -> f64 {
    PRESENT - self.tau[&node]
  }

  pub fn edge_durations(&self) -> BTreeMap<GraphEdgeKey, f64> {
    self
      .graph
      .get_edges()
      .map(|edge| (edge.key(), self.tau[&edge.source()] - self.tau[&edge.target()]))
      .collect()
  }
}

fn sorted_batches(sampling: &[SampleBatch]) -> Result<Vec<SampleBatch>, Report> {
  if sampling.is_empty() {
    return make_error!("Sampling must have at least one batch");
  }
  for (index, batch) in sampling.iter().enumerate() {
    if !(batch.tau.is_finite() && batch.tau >= 0.0) {
      return make_error!(
        "Sample batch {index} has time {}, expected a finite non-negative number of years before the most recent sample",
        batch.tau
      );
    }
    if batch.count == 0 {
      return make_error!("Sample batch {index} at time {} has no tips", batch.tau);
    }
  }
  let tip_count: usize = sampling.iter().map(|batch| batch.count).sum();
  if tip_count < 2 {
    return make_error!("Sampling must have at least 2 tips in total, got {tip_count}");
  }
  Ok(
    sampling
      .iter()
      .copied()
      .sorted_by(|left, right| left.tau.total_cmp(&right.tau))
      .collect(),
  )
}

fn add_batch(
  batch: &SampleBatch,
  graph: &mut Graph,
  tau: &mut BTreeMap<GraphNodeKey, f64>,
  tips: &mut Vec<GraphNodeKey>,
  active: &mut Vec<GraphNodeKey>,
) {
  for _ in 0..batch.count {
    let tip = graph.add_node();
    tau.insert(tip, batch.tau);
    tips.push(tip);
    active.push(tip);
  }
}

fn node_names(tips: &[GraphNodeKey], mergers: &[GraphNodeKey]) -> BTreeMap<GraphNodeKey, Option<String>> {
  let tip_names = tips.iter().enumerate().map(|(index, &tip)| {
    let number = index + 1;
    (tip, format!("tip_{number:04}"))
  });
  let merger_names = mergers.iter().enumerate().map(|(index, &merger)| {
    let number = index + 1;
    let name = if number == mergers.len() {
      "node_root".to_owned()
    } else {
      format!("node_{number:04}")
    };
    (merger, name)
  });
  tip_names
    .chain(merger_names)
    .map(|(node, name)| (node, Some(name)))
    .collect()
}
