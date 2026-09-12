#[cfg(test)]
mod __tests__;

use crate::dependency_queue::{run_dependency_queue, validate_dependency_graph};
use crate::edge::GraphEdgeKey;
use crate::graph::Graph;
use crate::node::GraphNodeKey;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;
use treetime_utils::make_internal_report;

/// Frozen, immutable topology view of a graph, prepared once in the caller before any worker runs.
///
/// The view carries every topology fact a parallel worker needs -- node keys, each node's parent edge,
/// and each node's children in the graph's own `children_of` (outbound) order -- so a worker never
/// reads the graph (and never takes a graph lock) during a pass. Build it once with [`GraphPass::new`],
/// then run [`GraphPass::map_backward`] or [`GraphPass::map_forward`] against borrowed input maps as
/// many times as needed. After a structural change to the graph, rebuild the view.
pub struct GraphPass {
  /// Nodes in `graph.get_nodes()` traversal order.
  nodes: Vec<GraphPassNode>,
  /// Node key -> index into `nodes`. A map (not a key-indexed vector) so deleted-key gaps in the
  /// graph's node storage cost nothing.
  node_index: BTreeMap<GraphNodeKey, usize>,
  /// Every edge key reachable as some node's parent edge, for input validation.
  edge_keys: BTreeSet<GraphEdgeKey>,
  /// Parent node index per node (`None` at a root).
  parents: Vec<Option<usize>>,
  /// Child node indices per node, in `children_of` (outbound edge) order. This is the canonical
  /// numerical child order the reductions fold in, so it is thread-count independent.
  children: Vec<Vec<usize>>,
}

impl GraphPass {
  /// Freeze the topology of `graph` into a reusable pass view and validate that it forms an acyclic
  /// dependency graph (each node has at most one parent, no cycles, no duplicate readiness).
  pub fn new(graph: &Graph) -> Result<Self, Report> {
    let safe_nodes = graph.get_nodes();
    let mut nodes = Vec::with_capacity(safe_nodes.len());
    for safe in &safe_nodes {
      let node = safe.read_arc();
      let key = node.key();
      let parent_edge = if let Some(edge_key) = graph.parent_inbound_edge(key)? {
        Some((graph.get_source_node_key(edge_key)?, edge_key))
      } else {
        None
      };
      nodes.push(GraphPassNode { key, parent_edge });
    }

    let node_index = nodes
      .iter()
      .enumerate()
      .map(|(index, node)| (node.key, index))
      .collect::<BTreeMap<_, _>>();

    let parents = nodes
      .iter()
      .map(|node| node.parent_edge.map(|(parent_key, _)| node_index[&parent_key]))
      .collect::<Vec<_>>();

    // Children in `children_of` (outbound) order, so the value the backward pass hands each visitor
    // already folds in the same canonical order the graph exposes, without any per-node graph read.
    let mut children = vec![Vec::new(); nodes.len()];
    for safe in &safe_nodes {
      let node = safe.read_arc();
      let parent_index = node_index[&node.key()];
      for (child, _edge) in graph.children_of(&node) {
        let child_key = child.read_arc().key();
        children[parent_index].push(node_index[&child_key]);
      }
    }

    let edge_keys = nodes
      .iter()
      .filter_map(|node| node.parent_edge.map(|(_, edge_key)| edge_key))
      .collect::<BTreeSet<_>>();

    // Validate the backward schedule: a valid rooted forest is acyclic in both traversal directions,
    // so this single check covers both maps.
    let prerequisites = children.iter().map(Vec::len).collect::<Vec<_>>();
    let successors = parents
      .iter()
      .map(|parent| parent.iter().copied().collect::<Vec<_>>())
      .collect::<Vec<_>>();
    validate_dependency_graph(&prerequisites, &successors)?;

    Ok(Self {
      nodes,
      node_index,
      edge_keys,
      parents,
      children,
    })
  }

  /// Run a value-returning backward map (children before parent) over borrowed input maps.
  ///
  /// Every node is visited exactly once, after all of its children have published their outputs, in a
  /// thread-count-independent child fold order. The input maps are read immutably and never mutated,
  /// so a failed visit leaves them intact and the caller may retry from the same inputs. Node keys the
  /// graph has but `nodes` lacks are filled by `missing_node`; edges must match the graph exactly.
  pub fn map_backward<N, E, NodeOut, EdgeOut>(
    &self,
    nodes: &BTreeMap<GraphNodeKey, N>,
    edges: &BTreeMap<GraphEdgeKey, E>,
    missing_node: impl FnMut(GraphNodeKey) -> Result<N, Report>,
    visit: impl Fn(
      GraphPassBackwardContext<'_, N, E, NodeOut, EdgeOut>,
    ) -> Result<GraphPassNodeOutput<NodeOut, EdgeOut>, Report>
    + Sync
    + Send,
  ) -> Result<GraphMapOutputs<NodeOut, EdgeOut>, Report>
  where
    N: Sync,
    E: Sync,
    NodeOut: Send + Sync,
    EdgeOut: Send + Sync,
  {
    let created = self.validate_and_create_missing(nodes, edges, missing_node)?;

    // Backward schedule: a node becomes ready once every child has completed, then unblocks its parent.
    let prerequisites = self.children.iter().map(Vec::len).collect::<Vec<_>>();
    let successors = self
      .parents
      .iter()
      .map(|parent| parent.iter().copied().collect::<Vec<_>>())
      .collect::<Vec<_>>();

    let completed = std::iter::repeat_with(OnceLock::<GraphPassNodeOutput<NodeOut, EdgeOut>>::new)
      .take(self.nodes.len())
      .collect::<Vec<_>>();

    run_dependency_queue(&prerequisites, &successors, |index| {
      let node = &self.nodes[index];
      let input = self.resolve_node(nodes, &created, index);
      let parent_edge = node.parent_edge.map(|(_, edge_key)| (edge_key, &edges[&edge_key]));

      // Every child has completed before this node is scheduled, so its output is published. Collect
      // them in the fixed `children_of` order recorded in `self.children[index]` so the reduction the
      // visitor folds does not depend on thread count.
      let children = self.children[index]
        .iter()
        .map(|&child_index| {
          let child = &self.nodes[child_index];
          let (_, edge_key) = child.parent_edge.expect("Backward child must have a parent edge");
          let output = completed[child_index]
            .get()
            .expect("Backward child must complete before its parent");
          GraphPassChildBackward {
            node_key: child.key,
            edge_key,
            node: &output.node,
            edge: output.parent_message.as_ref(),
          }
        })
        .collect::<Vec<_>>();

      let context = GraphPassBackwardContext {
        key: node.key,
        is_leaf: self.children[index].is_empty(),
        is_root: self.parents[index].is_none(),
        input,
        parent_edge,
        children: &children,
      };
      let output = visit(context)?;
      assert!(
        completed[index].set(output).is_ok(),
        "Dependency queue must publish each indexed slot once"
      );
      Ok(())
    })?;

    self.collect_map_outputs(completed)
  }

  /// Run a value-returning forward map (parent before children) over borrowed input maps.
  ///
  /// Every node is visited exactly once, after its single parent has published its output (roots run
  /// first). The input maps are read immutably and never mutated, so a failed visit leaves them intact
  /// and the caller may retry from the same inputs.
  pub fn map_forward<N, E, NodeOut, EdgeOut>(
    &self,
    nodes: &BTreeMap<GraphNodeKey, N>,
    edges: &BTreeMap<GraphEdgeKey, E>,
    missing_node: impl FnMut(GraphNodeKey) -> Result<N, Report>,
    visit: impl Fn(GraphPassForwardContext<'_, N, E, NodeOut>) -> Result<GraphPassNodeOutput<NodeOut, EdgeOut>, Report>
    + Sync
    + Send,
  ) -> Result<GraphMapOutputs<NodeOut, EdgeOut>, Report>
  where
    N: Sync,
    E: Sync,
    NodeOut: Send + Sync,
    EdgeOut: Send + Sync,
  {
    let created = self.validate_and_create_missing(nodes, edges, missing_node)?;

    // Forward schedule: a node becomes ready once its parent has completed (roots are ready
    // immediately), then unblocks its children.
    let prerequisites = self
      .parents
      .iter()
      .map(|parent| usize::from(parent.is_some()))
      .collect::<Vec<_>>();
    let successors = self.children.clone();

    let completed = std::iter::repeat_with(OnceLock::<GraphPassNodeOutput<NodeOut, EdgeOut>>::new)
      .take(self.nodes.len())
      .collect::<Vec<_>>();

    run_dependency_queue(&prerequisites, &successors, |index| {
      let node = &self.nodes[index];
      let input = self.resolve_node(nodes, &created, index);
      let parent_edge = node.parent_edge.map(|(_, edge_key)| (edge_key, &edges[&edge_key]));

      // The parent has completed before this node is scheduled, so its output is published.
      let parent = self.parents[index].map(|parent_index| {
        &completed[parent_index]
          .get()
          .expect("Forward parent must complete before its child")
          .node
      });

      let context = GraphPassForwardContext {
        key: node.key,
        is_leaf: self.children[index].is_empty(),
        is_root: self.parents[index].is_none(),
        input,
        parent_edge,
        parent,
      };
      let output = visit(context)?;
      assert!(
        completed[index].set(output).is_ok(),
        "Dependency queue must publish each indexed slot once"
      );
      Ok(())
    })?;

    self.collect_map_outputs(completed)
  }

  /// Validate the input maps against the frozen topology and create inputs for node keys the graph has
  /// but `nodes` lacks. Stale node/edge keys, or an edge set that does not match the graph exactly, are
  /// internal errors. Returns the created (missing) node inputs, owned so they can be borrowed by the
  /// workers alongside `nodes`.
  fn validate_and_create_missing<N, E>(
    &self,
    nodes: &BTreeMap<GraphNodeKey, N>,
    edges: &BTreeMap<GraphEdgeKey, E>,
    mut missing_node: impl FnMut(GraphNodeKey) -> Result<N, Report>,
  ) -> Result<BTreeMap<GraphNodeKey, N>, Report> {
    if nodes.keys().any(|key| !self.node_index.contains_key(key))
      || edges.keys().any(|key| !self.edge_keys.contains(key))
    {
      return Err(make_internal_report!(
        "Partition contains stale topology entries while indexing a pass"
      ));
    }
    if edges.len() != self.edge_keys.len() || self.edge_keys.iter().any(|key| !edges.contains_key(key)) {
      return Err(make_internal_report!(
        "Partition edge set does not match the graph while indexing a pass"
      ));
    }
    let mut created = BTreeMap::new();
    for node in &self.nodes {
      if !nodes.contains_key(&node.key) {
        created.insert(node.key, missing_node(node.key)?);
      }
    }
    Ok(created)
  }

  /// Borrow the input for node `index`, from the caller's map when present or from the created
  /// missing-node inputs otherwise.
  fn resolve_node<'a, N>(
    &self,
    nodes: &'a BTreeMap<GraphNodeKey, N>,
    created: &'a BTreeMap<GraphNodeKey, N>,
    index: usize,
  ) -> &'a N {
    let key = self.nodes[index].key;
    nodes
      .get(&key)
      .or_else(|| created.get(&key))
      .expect("Indexed node must have an input")
  }

  /// Drain the per-index published outputs into key-addressed maps: each node output keyed by its node,
  /// and each node's optional parent-edge message keyed by that parent edge. Duplicate keys signal a
  /// scheduling bug, so they are reported as internal errors.
  fn collect_map_outputs<NodeOut, EdgeOut>(
    &self,
    completed: Vec<OnceLock<GraphPassNodeOutput<NodeOut, EdgeOut>>>,
  ) -> Result<GraphMapOutputs<NodeOut, EdgeOut>, Report> {
    let mut nodes = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for (index, slot) in completed.into_iter().enumerate() {
      let output = slot
        .into_inner()
        .expect("Every indexed slot must publish an output after a successful map");
      let key = self.nodes[index].key;
      if nodes.insert(key, output.node).is_some() {
        return Err(make_internal_report!(
          "Duplicate node {key} while collecting a graph map"
        ));
      }
      if let Some(edge) = output.parent_message {
        let (_, edge_key) = self.nodes[index]
          .parent_edge
          .expect("Parent message must belong to a parent edge");
        if edges.insert(edge_key, edge).is_some() {
          return Err(make_internal_report!(
            "Duplicate edge {edge_key} while collecting a graph map"
          ));
        }
      }
    }
    Ok(GraphMapOutputs { nodes, edges })
  }
}

/// One completed child seen by a backward-mapping visitor: the child's returned node output and the
/// optional message it produced for the edge connecting it to the current (parent) node. Children are
/// presented in the graph's `children_of` (outbound) order.
pub struct GraphPassChildBackward<'a, NodeOut, EdgeOut> {
  pub node_key: GraphNodeKey,
  pub edge_key: GraphEdgeKey,
  pub node: &'a NodeOut,
  pub edge: Option<&'a EdgeOut>,
}

/// Input handed to a backward-mapping visitor for one node: the node's own borrowed input, its borrowed
/// parent-edge input, and the already-completed outputs of its children (in `children_of` order).
pub struct GraphPassBackwardContext<'a, N, E, NodeOut, EdgeOut> {
  pub key: GraphNodeKey,
  pub is_leaf: bool,
  pub is_root: bool,
  pub input: &'a N,
  pub parent_edge: Option<(GraphEdgeKey, &'a E)>,
  pub children: &'a [GraphPassChildBackward<'a, NodeOut, EdgeOut>],
}

/// Input handed to a forward-mapping visitor for one node: the node's own borrowed input, its borrowed
/// parent-edge input, and the already-completed forward output of its single parent.
pub struct GraphPassForwardContext<'a, N, E, NodeOut> {
  pub key: GraphNodeKey,
  pub is_leaf: bool,
  pub is_root: bool,
  pub input: &'a N,
  pub parent_edge: Option<(GraphEdgeKey, &'a E)>,
  pub parent: Option<&'a NodeOut>,
}

/// Output returned by a mapping visitor for one node: the node's output and the optional message it
/// sends along its own parent edge (`None` at the root, which has no parent edge, or whenever the node
/// produces no message for that edge).
pub struct GraphPassNodeOutput<NodeOut, EdgeOut> {
  pub node: NodeOut,
  pub parent_message: Option<EdgeOut>,
}

/// Collected outputs of a graph map: node outputs keyed by node, and per-edge messages keyed by the
/// edge each message travelled along.
pub struct GraphMapOutputs<NodeOut, EdgeOut> {
  pub nodes: BTreeMap<GraphNodeKey, NodeOut>,
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
}

/// A single node in the frozen pass topology: its key and, for a non-root, the parent node and the edge
/// connecting them.
struct GraphPassNode {
  key: GraphNodeKey,
  parent_edge: Option<(GraphNodeKey, GraphEdgeKey)>,
}
