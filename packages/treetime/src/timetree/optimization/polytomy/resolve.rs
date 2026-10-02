use crate::optimize::topology::polytomy_nodes::find_polytomy_nodes;
use crate::timetree::branch_model::BranchModel;
use crate::timetree::optimization::polytomy::apply::{ChildRef, apply_plan};
use crate::timetree::optimization::polytomy::sweep::{Lineage, simulate_subtree};
use eyre::{Report, WrapErr};
use log::debug;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::{record_merge, remove_node_if_trivial, trivial_node_branch_lengths};
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_utils::make_error;

pub(crate) fn resolve_polytomies(
  mut graph: Graph,
  mut branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  branch_model: &BranchModel,
  mutation_rate: f64,
  total_length: usize,
  merger_rate: &PiecewiseConstantFn,
  rng: &mut dyn rand::RngCore,
  node_times: &BTreeMap<GraphNodeKey, Option<f64>>,
) -> Result<PolytomyResolution, Report> {
  let polytomy_keys = find_polytomy_nodes(&graph);
  if polytomy_keys.is_empty() {
    debug!("No polytomies to resolve");
    return Ok(PolytomyResolution {
      graph,
      branch_lengths,
      merger_times: BTreeMap::new(),
    });
  }

  let mut merger_times = BTreeMap::new();
  let mut topology_validated = false;

  for node_key in polytomy_keys {
    let created = resolve_single_polytomy(
      &mut graph,
      branch_model,
      node_key,
      mutation_rate,
      total_length,
      merger_rate,
      rng,
      &mut topology_validated,
      &mut branch_lengths,
      node_times,
    )?;
    merger_times.extend(created);
  }

  let obsolete_count = remove_single_child_nodes(&mut graph, &mut branch_lengths)?;
  if obsolete_count > 0 {
    debug!("Removed {obsolete_count} obsolete single-child nodes");
  }

  graph.build()?;

  if merger_times.is_empty() {
    debug!("Polytomies found but the sampled histories resolved none of them");
  } else {
    debug!("Polytomy resolution introduced {} new nodes", merger_times.len());
  }

  Ok(PolytomyResolution {
    graph,
    branch_lengths,
    merger_times,
  })
}

pub(crate) struct PolytomyResolution {
  pub graph: Graph,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub merger_times: BTreeMap<GraphNodeKey, f64>,
}

#[allow(
  clippy::too_many_arguments,
  reason = "one call site; splitting would only shuffle the arguments"
)]
fn resolve_single_polytomy(
  graph: &mut Graph,
  branch_model: &BranchModel,
  node_key: GraphNodeKey,
  mutation_rate: f64,
  total_length: usize,
  merger_rate: &PiecewiseConstantFn,
  rng: &mut dyn rand::RngCore,
  topology_validated: &mut bool,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  node_times: &BTreeMap<GraphNodeKey, Option<f64>>,
) -> Result<BTreeMap<GraphNodeKey, f64>, Report> {
  let parent_time = inferred_time(node_times, node_key)?;

  let children = collect_children(graph, branch_model, node_key, total_length, branch_lengths, node_times)?;
  if children.len() < 3 {
    return Ok(BTreeMap::new());
  }

  let lineages: Vec<Lineage> = children
    .iter()
    .map(|child| Lineage {
      time: child.time,
      mutations: child.mutations,
    })
    .collect();

  let plan = simulate_subtree(&lineages, parent_time, mutation_rate, merger_rate, rng)?;
  if plan.mergers.is_empty() {
    debug!(
      "Polytomy at node {node_key}: {} children, sampled history merged none",
      children.len()
    );
    return Ok(BTreeMap::new());
  }

  if !*topology_validated {
    require_internal_node_times(graph, node_times)?;
    *topology_validated = true;
  }

  let child_refs: Vec<ChildRef> = children
    .iter()
    .map(|child| ChildRef {
      edge_key: child.edge_key,
      time: child.time,
    })
    .collect();

  let created = apply_plan(graph, node_key, parent_time, &child_refs, &plan, branch_lengths)?;

  debug!(
    "Polytomy at node {node_key}: {} children -> {} children, created {} nodes",
    children.len(),
    plan.roots.len(),
    created.len()
  );

  Ok(created)
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
fn collect_children(
  graph: &Graph,
  branch_model: &BranchModel,
  node_key: GraphNodeKey,
  total_length: usize,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  node_times: &BTreeMap<GraphNodeKey, Option<f64>>,
) -> Result<Vec<ChildInfo>, Report> {
  let edge_keys = {
    let node = graph.get_node(node_key).expect("Node must exist");
    node.outbound().to_vec()
  };

  edge_keys
    .into_iter()
    .map(|edge_key| {
      let edge = graph.get_edge(edge_key).expect("Edge must exist");
      let child_key = edge.target();

      let time = inferred_time(node_times, child_key)?;
      let mutation_length = branch_lengths.get(&edge_key).copied().flatten();

      Ok(ChildInfo {
        edge_key,
        time,
        mutations: edge_mutation_count(graph, branch_model, edge_key, mutation_length, total_length)?,
      })
    })
    .collect()
}

struct ChildInfo {
  edge_key: GraphEdgeKey,
  time: f64,
  mutations: u32,
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn edge_mutation_count(
  graph: &Graph,
  branch_model: &BranchModel,
  edge_key: GraphEdgeKey,
  mutation_length: Option<f64>,
  total_length: usize,
) -> Result<u32, Report> {
  let exact = match branch_model {
    BranchModel::Input => None,
    BranchModel::Marginal(partition) => partition
      .edge_sub_count(graph, edge_key)
      .wrap_err_with(|| format!("When counting the substitutions on edge {edge_key}"))?,
  };

  let count = exact.unwrap_or_else(|| {
    let estimate = mutation_length.unwrap_or(0.0) * total_length as f64;
    if estimate.is_finite() && estimate > 0.0 {
      estimate.round() as usize
    } else {
      0
    }
  });

  Ok(u32::try_from(count).unwrap_or(u32::MAX))
}

fn inferred_time(node_times: &BTreeMap<GraphNodeKey, Option<f64>>, node_key: GraphNodeKey) -> Result<f64, Report> {
  let Some(time) = node_times[&node_key] else {
    return make_error!("Polytomy resolution requires an inferred time for node {node_key}, but it has none");
  };
  if !time.is_finite() {
    return make_error!("Polytomy resolution requires a finite inferred time for node {node_key}, but it has {time}");
  }
  Ok(time)
}

fn remove_single_child_nodes(
  graph: &mut Graph,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<usize, Report> {
  let mut removed_count = 0;

  loop {
    let obsolete_key = graph.get_nodes().into_iter().find_map(|node| {
      let is_trivial = node.inbound().len() == 1 && node.outbound().len() == 1;
      is_trivial.then_some(node.key())
    });

    let Some(node_key) = obsolete_key else {
      break;
    };

    let (parent_branch, child_branch) = trivial_node_branch_lengths(graph, node_key, branch_lengths);
    if let Some(merge) = remove_node_if_trivial(graph, node_key, parent_branch, child_branch)? {
      record_merge(branch_lengths, &merge);
      removed_count += 1;
    }
  }

  Ok(removed_count)
}

pub(crate) fn require_internal_node_times(
  graph: &Graph,
  node_times: &BTreeMap<GraphNodeKey, Option<f64>>,
) -> Result<(), Report> {
  for node in graph.get_nodes() {
    if node.is_leaf() {
      continue;
    }
    let Some(time) = node_times[&node.key()] else {
      return make_error!(
        "Topology rebuild requires an inferred time for every internal node, but node {:?} has none",
        node.key()
      );
    };
    if !time.is_finite() {
      return make_error!(
        "Topology rebuild requires a finite inferred time for every internal node, but node {:?} has {time}",
        node.key()
      );
    }
  }
  Ok(())
}
