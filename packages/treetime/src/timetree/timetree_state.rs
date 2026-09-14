use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::node_time::{CoalescentNodeTime, CoalescentNodeTimes};
use eyre::Report;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{
  GraphMapOutputs, GraphPass, GraphPassBackwardContext, GraphPassForwardContext, GraphPassNodeOutput,
};

/// Per-node date-inference fields, held as a value keyed by [`GraphNodeKey`].
///
/// `time_distribution` is the node's posterior over its date, refined in place across the two date
/// passes; `time` the committed point estimate; `bad_branch` the exclusion flag a parent reads when
/// gathering child messages; `contradicted` a per-pass flag the forward pass raises when the rest of
/// the tree gives the given date no probability, folded into a diagnostic count.
///
/// The fixed input date lives outside this working state, in the durable
/// [`DateConstraints`](crate::clock::date_constraints::DateConstraints) the passes borrow: the
/// backward pass lifts it back into the posterior on every pass, so keeping it out of the produced
/// results is what keeps the input recoverable and separately owned.
#[derive(Debug, Clone, Default)]
pub struct DateNodeState {
  pub time_distribution: Option<Arc<Distribution<NegLog>>>,
  pub time: Option<f64>,
  pub bad_branch: bool,
  pub contradicted: bool,
}

/// Per-edge date-inference fields, held as a value keyed by [`GraphEdgeKey`].
///
/// `branch_length_distribution` is the branch's time-duration law the passes convolve across;
/// `msg_to_parent` the backward message the child sends up, divided back out as the cavity on the
/// forward pass; `time_length` the branch's committed duration (the Newick weight); `gamma` the
/// per-branch relaxed-clock rate multiplier, `1.0` for a strict clock.
#[derive(Debug, Clone, SmartDefault)]
pub struct DateEdgeState {
  pub branch_length_distribution: Option<Arc<Distribution<NegLog>>>,
  pub msg_to_parent: Option<Arc<Distribution<NegLog>>>,
  pub time_length: Option<f64>,
  #[default = 1.0]
  pub gamma: f64,
}

/// The date-inference state for a whole tree, routed through the timetree date passes as the
/// per-node [`DateNodeState`] and per-edge [`DateEdgeState`] fields.
///
/// Keyed by stable node and edge ids, so the maps stay valid across a reroot and polytomy resolution,
/// which add and drop nodes and edges while leaving ids stable (gaps, never renumbered).
#[derive(Debug, Clone, Default)]
pub struct TimetreeState {
  pub nodes: BTreeMap<GraphNodeKey, DateNodeState>,
  pub edges: BTreeMap<GraphEdgeKey, DateEdgeState>,
}

impl TimetreeState {
  /// Empty per-node/per-edge state for every node and edge of `graph`, all fields default.
  pub fn new(graph: &Graph) -> Self {
    let nodes = graph
      .get_nodes()
      .iter()
      .map(|node| (node.read_arc().key(), DateNodeState::default()))
      .collect();
    let edges = graph
      .get_edges()
      .iter()
      .map(|edge| (edge.read_arc().key(), DateEdgeState::default()))
      .collect();
    Self { nodes, edges }
  }

  /// Seed date state from the value maps [`load_date_constraints`] returns.
  ///
  /// At the pipeline seed point each node's date constraint, initial time distribution (equal to the
  /// constraint), and bad-branch flag come from `constraints`; the committed time starts `None`,
  /// `contradicted` starts false, and every edge starts default (no branch-length distribution or
  /// backward message, `None` time length, strict-clock `gamma`).
  ///
  /// [`load_date_constraints`]: crate::clock::date_constraints::load_date_constraints
  pub fn seed_from_values(graph: &Graph, constraints: &DateConstraints) -> Self {
    let nodes = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let key = node.read_arc().key();
        let state = DateNodeState {
          time_distribution: constraints.time_distributions.get(&key).cloned().flatten(),
          time: None,
          bad_branch: constraints.bad_branches.get(&key).copied().unwrap_or(false),
          contradicted: false,
        };
        (key, state)
      })
      .collect();
    let edges = graph
      .get_edges()
      .iter()
      .map(|edge| (edge.read_arc().key(), DateEdgeState::default()))
      .collect();
    Self { nodes, edges }
  }

  /// Rebuild the per-node and per-edge maps to match the current graph, carrying every value-resident
  /// field forward.
  ///
  /// Stays valid across a reroot or polytomy resolution that added or dropped nodes and edges: each
  /// surviving node and edge keeps its bad-branch flag, committed time, time distribution, committed
  /// time length, branch-length distribution, backward message, and relaxed-clock rate multiplier
  /// from the previous state, and a node or edge a topology change introduced starts default. The
  /// fixed date constraint lives in the durable
  /// [`DateConstraints`](crate::clock::date_constraints::DateConstraints), not here. `contradicted` is
  /// a per-pass flag and always resets to false.
  ///
  /// The bad-branch flag and committed time length are written straight into this state by their
  /// producers -- the clock filter
  /// ([`apply_outlier_bad_branches`](crate::timetree::optimization::clock_filter::apply_outlier_bad_branches)),
  /// the topology rebuild
  /// ([`propagate_bad_branches`](crate::timetree::optimization::clock_filter::propagate_bad_branches)),
  /// and the branch-distribution builders and polytomy application for the time length.
  pub fn reseed_from_values(&mut self, graph: &Graph) {
    let nodes = graph
      .get_nodes()
      .iter()
      .map(|node| {
        let key = node.read_arc().key();
        let state = self
          .nodes
          .get(&key)
          .map_or_else(DateNodeState::default, |node| DateNodeState {
            time_distribution: node.time_distribution.clone(),
            time: node.time,
            bad_branch: node.bad_branch,
            contradicted: false,
          });
        (key, state)
      })
      .collect();
    let edges = graph
      .get_edges()
      .iter()
      .map(|edge| {
        let key = edge.read_arc().key();
        let state = self
          .edges
          .get(&key)
          .map_or_else(DateEdgeState::default, |edge| DateEdgeState {
            branch_length_distribution: edge.branch_length_distribution.clone(),
            msg_to_parent: edge.msg_to_parent.clone(),
            time_length: edge.time_length,
            gamma: edge.gamma,
          });
        (key, state)
      })
      .collect();
    self.nodes = nodes;
    self.edges = edges;
  }

  /// Clear the value-resident edge fields after a topology change, so the next pass rebuilds them.
  ///
  /// Re-parenting invalidates the branch-length distributions, backward messages, and relaxed-clock
  /// rate multiplier: they describe a parent-child pair that no longer exists. This blanks the
  /// distribution and message and resets `gamma` to the strict-clock `1.0` on every current edge (and
  /// adds default entries for edges and nodes the topology change introduced), the counterpart of the
  /// reset
  /// [`prepare_tree_after_topology_change`](crate::timetree::optimization::polytomy::prepare_tree_after_topology_change)
  /// does for the transitional fields. The following [`reseed_from_values`](Self::reseed_from_values) preserves
  /// these blanked values, so the branch-distribution builders start each surviving edge from `None`.
  pub fn reset_date_edges_for_topology_change(&mut self, graph: &Graph) {
    for edge_ref in graph.get_edges() {
      let key = edge_ref.read_arc().key();
      let entry = self.edges.entry(key).or_default();
      entry.branch_length_distribution = None;
      entry.msg_to_parent = None;
      entry.gamma = 1.0;
    }
    for node_ref in graph.get_nodes() {
      let key = node_ref.read_arc().key();
      self.nodes.entry(key).or_default();
    }
  }

  /// Per-node date the clock regression reads, keyed by node, computed from this state and the
  /// durable `constraints`.
  ///
  /// The likely time: the input date constraint where there is one, the refined time distribution's
  /// peak otherwise. Used to seed the clock inputs in the refinement loop.
  #[must_use]
  pub fn likely_times(&self, constraints: &DateConstraints) -> BTreeMap<GraphNodeKey, Option<f64>> {
    self
      .nodes
      .iter()
      .map(|(key, node)| {
        let time = constraints
          .date_constraints
          .get(key)
          .cloned()
          .flatten()
          .as_ref()
          .or(node.time_distribution.as_ref())
          .and_then(|dist| dist.likely_time());
        (*key, time)
      })
      .collect()
  }

  /// Build the coalescent node-time map from this state, so the coalescent collectors read node
  /// times as a value. Each entry carries both the committed point
  /// estimate and the distribution peak.
  #[must_use]
  pub fn coalescent_node_times(&self) -> CoalescentNodeTimes {
    self
      .nodes
      .iter()
      .map(|(key, node)| {
        let entry = CoalescentNodeTime {
          time: node.time,
          time_dist_likely: node.time_distribution.as_ref().and_then(|dist| dist.likely_time()),
          bad_branch: node.bad_branch,
        };
        (*key, entry)
      })
      .collect()
  }

  #[must_use]
  pub fn node(&self, key: GraphNodeKey) -> &DateNodeState {
    self
      .nodes
      .get(&key)
      .unwrap_or_else(|| panic!("Timetree state is missing node {key}"))
  }

  #[must_use]
  pub fn node_mut(&mut self, key: GraphNodeKey) -> &mut DateNodeState {
    self
      .nodes
      .get_mut(&key)
      .unwrap_or_else(|| panic!("Timetree state is missing node {key}"))
  }

  #[must_use]
  pub fn edge(&self, key: GraphEdgeKey) -> &DateEdgeState {
    self
      .edges
      .get(&key)
      .unwrap_or_else(|| panic!("Timetree state is missing edge {key}"))
  }

  #[must_use]
  pub fn edge_mut(&mut self, key: GraphEdgeKey) -> &mut DateEdgeState {
    self
      .edges
      .get_mut(&key)
      .unwrap_or_else(|| panic!("Timetree state is missing edge {key}"))
  }

  /// Run a value-returning backward pass over the date state through the graph's dependency engine,
  /// replacing the per-node and per-edge maps with the visitor's outputs. Uses a
  /// thread-independent, deterministic child fold order.
  pub fn map_backward<F>(&mut self, graph: &Graph, visit: F) -> Result<(), Report>
  where
    F: Fn(
        GraphPassBackwardContext<'_, DateNodeState, DateEdgeState, DateNodeState, DateEdgeState>,
      ) -> Result<GraphPassNodeOutput<DateNodeState, DateEdgeState>, Report>
      + Sync
      + Send,
  {
    let pass = GraphPass::new(graph)?;
    let GraphMapOutputs { nodes, edges } =
      pass.map_backward(&self.nodes, &self.edges, |_| Ok(DateNodeState::default()), visit)?;
    self.nodes = nodes;
    self.edges = edges;
    Ok(())
  }

  /// Run a value-returning forward pass over the date state through the graph's dependency engine,
  /// replacing the per-node and per-edge maps with the visitor's outputs. Each node reads its single
  /// parent's already-published output.
  pub fn map_forward<F>(&mut self, graph: &Graph, visit: F) -> Result<(), Report>
  where
    F: Fn(
        GraphPassForwardContext<'_, DateNodeState, DateEdgeState, DateNodeState>,
      ) -> Result<GraphPassNodeOutput<DateNodeState, DateEdgeState>, Report>
      + Sync
      + Send,
  {
    let pass = GraphPass::new(graph)?;
    let GraphMapOutputs { nodes, edges } =
      pass.map_forward(&self.nodes, &self.edges, |_| Ok(DateNodeState::default()), visit)?;
    self.nodes = nodes;
    self.edges = edges;
    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use treetime_distribution::Distribution;
  use treetime_io::nwk::nwk_read_str;

  /// A topology change clears the value-resident branch-length distribution and backward message on
  /// every edge, so the next branch-distribution build starts each surviving edge from scratch.
  #[test]
  fn test_timetree_state_reset_date_edges_clears_distribution_and_message() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:1.0,B:1.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut state = TimetreeState::new(&graph);
    for edge_ref in graph.get_edges() {
      let key = edge_ref.read_arc().key();
      let entry = state.edge_mut(key);
      entry.branch_length_distribution = Some(Arc::new(Distribution::point(1.0, 0.0)));
      entry.msg_to_parent = Some(Arc::new(Distribution::point(2.0, 0.0)));
    }

    state.reset_date_edges_for_topology_change(&graph);

    for edge_ref in graph.get_edges() {
      let key = edge_ref.read_arc().key();
      let entry = state.edge(key);
      assert_eq!(None, entry.branch_length_distribution);
      assert_eq!(None, entry.msg_to_parent);
    }

    Ok(())
  }

  /// Rebuilding the maps keeps the value-resident branch-length distribution and backward message for
  /// edges already in the state, so they carry across passes.
  #[test]
  fn test_timetree_state_reseed_preserves_distribution_and_message() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:1.0,B:1.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut state = TimetreeState::new(&graph);
    let key = graph
      .get_edges()
      .first()
      .expect("tree has at least one edge")
      .read_arc()
      .key();
    let dist = Arc::new(Distribution::point(3.0, 0.0));
    let msg = Arc::new(Distribution::point(4.0, 0.0));
    {
      let entry = state.edge_mut(key);
      entry.branch_length_distribution = Some(Arc::clone(&dist));
      entry.msg_to_parent = Some(Arc::clone(&msg));
    }

    state.reseed_from_values(&graph);

    let entry = state.edge(key);
    assert_eq!(Some(dist), entry.branch_length_distribution);
    assert_eq!(Some(msg), entry.msg_to_parent);

    Ok(())
  }
}
