use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::node_time::{CoalescentNodeTime, CoalescentNodeTimes};
use crate::timetree::inference::bad_branches::{derive_bad_branches, undated_leaves};
use crate::timetree::inference::time_inference::{
  BranchLikelihood, NodePosterior, TimeBackward, TimeInference, likely_times,
};
use eyre::Report;
use treetime_graph::graph::Graph;

pub(crate) fn constraint_coalescent_node_times(
  graph: &Graph,
  constraints: &DateConstraints,
) -> Result<CoalescentNodeTimes, Report> {
  let bad_branches = derive_bad_branches(graph, constraints, &undated_leaves(graph, constraints))?;
  let times = likely_times(graph, constraints, None)?
    .into_iter()
    .map(|(key, time_dist_likely)| {
      let entry = CoalescentNodeTime {
        time: None,
        time_dist_likely,
        bad_branch: bad_branches[&key],
      };
      (key, entry)
    })
    .collect();
  Ok(times)
}

pub(crate) fn empty_time_inference(graph: &Graph) -> TimeInference {
  TimeInference {
    bad_branches: graph.get_nodes().map(|node| (node.key(), false)).collect(),
    branches: graph
      .get_edges()
      .map(|edge| {
        let branch = BranchLikelihood {
          distribution: None,
          time_length: None,
        };
        (edge.key(), branch)
      })
      .collect(),
    backward: TimeBackward {
      subtree: graph.get_nodes().map(|node| (node.key(), None)).collect(),
      messages: graph.get_edges().map(|edge| (edge.key(), None)).collect(),
    },
    posterior: graph
      .get_nodes()
      .map(|node| (node.key(), NodePosterior::default()))
      .collect(),
  }
}
