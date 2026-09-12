use crate::gtr::gtr::GTR;
use crate::make_internal_error;
use crate::make_internal_report;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::seq::indel::InDel;
use crate::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub, mutation_event_strings};
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::NodeCommentProvider;
use treetime_primitives::Seq;

/// Access to a partition's substitution model and sequence length, shared by the concrete
/// representations so rate normalization can stay generic over them.
pub trait HasGtr {
  fn gtr(&self) -> &GTR;
  fn gtr_mut(&mut self) -> &mut GTR;
  fn sequence_length(&self) -> usize;

  fn weighted_rate(&self) -> f64 {
    self.sequence_length() as f64 * self.gtr().mu
  }

  fn normalize_rate(&mut self, scale: f64) {
    self.gtr_mut().mu /= scale;
  }
}

/// Minimal graph-structure abstraction used by per-branch partition operations.
///
/// Exists so that read accessors can serve any phylogenetic graph without binding to a concrete
/// `Graph`. Only the operations actually needed by branch-level computations are exposed: resolving an
/// edge to its endpoints and walking one step toward the root. Trait-object safe so `&dyn
/// BranchTopology` can flow through dynamic dispatch while callers keep their concrete `&Graph`.
pub trait BranchTopology: Send + Sync {
  /// Return `(parent_node_key, child_node_key)` for one edge.
  fn edge_endpoints(&self, edge_key: GraphEdgeKey) -> Result<(GraphNodeKey, GraphNodeKey), Report>;

  /// Return `Some((parent_node_key, parent_edge_key))` for a non-root node,
  /// or `None` when the node is the root. Errors when the node has more than
  /// one parent (the algorithm only supports trees).
  fn node_parent(&self, node_key: GraphNodeKey) -> Result<Option<(GraphNodeKey, GraphEdgeKey)>, Report>;

  /// Return the key of the single root node. Errors when the graph has zero or
  /// more than one root.
  fn root_key(&self) -> Result<GraphNodeKey, Report>;
}

impl BranchTopology for Graph {
  fn edge_endpoints(&self, edge_key: GraphEdgeKey) -> Result<(GraphNodeKey, GraphNodeKey), Report> {
    let edge = self
      .get_edge(edge_key)
      .ok_or_else(|| make_internal_report!("Edge {edge_key} not found"))?;
    let edge = edge.read_arc();
    Ok((edge.source(), edge.target()))
  }

  fn node_parent(&self, node_key: GraphNodeKey) -> Result<Option<(GraphNodeKey, GraphEdgeKey)>, Report> {
    let node = self
      .get_node(node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} not found"))?;
    let node = node.read_arc();
    let inbound = node.inbound();
    match inbound.len() {
      0 => Ok(None),
      1 => {
        let parent_edge_key = inbound[0];
        let parent_node_key = self.get_source_node_key(parent_edge_key)?;
        Ok(Some((parent_node_key, parent_edge_key)))
      },
      n => make_internal_error!("Node {node_key} has {n} parents; only trees are supported"),
    }
  }

  fn root_key(&self) -> Result<GraphNodeKey, Report> {
    Ok(self.get_exactly_one_root()?.read_arc().key())
  }
}

/// Read accessors over a completed marginal reconstruction, shared by dense and sparse representations.
///
/// Implemented by short-lived per-representation read views that borrow the durable partition inputs
/// together with the node states and edge messages/estimates the passes returned. The view is a
/// transient read projection assembled at a consumer boundary, never a stored stage-filled object.
///
/// Requires marginal inference to have run. Pre-marginal consumers (GTR inference, prune/merge) read
/// Fitch data directly.
pub trait PartitionBranchOps: Send + Sync {
  /// Return the sequence length represented by this partition.
  fn sequence_length(&self) -> usize;

  /// Return MAP-derived nucleotide substitutions for one edge.
  fn edge_subs(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report>;

  /// Return grouped aligned insertions and deletions for one edge.
  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel>;

  /// Return the reconstructed root sequence represented by this partition.
  fn root_sequence(&self, graph: &dyn BranchTopology) -> Result<Seq, Report>;

  /// Return the reconstructed sequence for one node.
  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq;

  fn edge_mutations(
    &self,
    graph: &dyn BranchTopology,
    edge_key: GraphEdgeKey,
    track: MutationTrack,
  ) -> Result<Vec<Mutation>, Report> {
    self
      .edge_subs(graph, edge_key)?
      .into_iter()
      .map(|substitution| Ok(Mutation::substitution(track.clone(), substitution)))
      .chain(
        self
          .edge_indels(edge_key)
          .iter()
          .map(|indel| Mutation::indel(track.clone(), indel)),
      )
      .collect()
  }

  /// Return the number of alignment positions where both parent and child
  /// have canonical (non-gap, non-ambiguous) states for one edge.
  fn edge_effective_length(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<usize, Report>;
}

pub struct MutationCommentProvider<'a> {
  partition: &'a dyn PartitionBranchOps,
  graph: &'a dyn BranchTopology,
}

impl<'a> MutationCommentProvider<'a> {
  pub fn new(partition: &'a dyn PartitionBranchOps, graph: &'a dyn BranchTopology) -> Self {
    Self { partition, graph }
  }
}

impl NodeCommentProvider for MutationCommentProvider<'_> {
  fn node_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report> {
    let Some((_parent_key, edge_key)) = self.graph.node_parent(key)? else {
      return Ok(BTreeMap::new());
    };
    let mut mutations = self
      .partition
      .edge_mutations(self.graph, edge_key, MutationTrack::Nucleotide)?;
    if mutations.is_empty() {
      return Ok(BTreeMap::new());
    }
    mutations.sort_by_key(|mutation| match &mutation.event {
      MutationEvent::Substitution(substitution) => substitution.pos(),
      MutationEvent::Insertion(segment) | MutationEvent::Deletion(segment) => segment.range.0,
    });
    let mutations = mutations
      .iter()
      .map(|mutation| mutation_event_strings(&mutation.event))
      .collect::<Result<Vec<_>, _>>()?
      .into_iter()
      .flatten()
      .join(",");
    Ok(btreemap! {
      "mutations".to_owned() => mutations,
    })
  }
}

/// Newick/Nexus node-comment provider that reads a gathered per-edge nucleotide mutation map.
///
/// Mirrors [`MutationCommentProvider`], but reads mutations from a value map instead of the partition,
/// so the tree writers no longer touch the partition during serialization.
pub struct EdgeMutationCommentProvider<'a> {
  edge_mutations: &'a BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  graph: &'a dyn BranchTopology,
}

impl<'a> EdgeMutationCommentProvider<'a> {
  pub fn new(edge_mutations: &'a BTreeMap<GraphEdgeKey, Vec<Mutation>>, graph: &'a dyn BranchTopology) -> Self {
    Self { edge_mutations, graph }
  }
}

impl NodeCommentProvider for EdgeMutationCommentProvider<'_> {
  fn node_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report> {
    let Some((_parent_key, edge_key)) = self.graph.node_parent(key)? else {
      return Ok(BTreeMap::new());
    };
    let mut mutations = self.edge_mutations[&edge_key].clone();
    if mutations.is_empty() {
      return Ok(BTreeMap::new());
    }
    mutations.sort_by_key(|mutation| match &mutation.event {
      MutationEvent::Substitution(substitution) => substitution.pos(),
      MutationEvent::Insertion(segment) | MutationEvent::Deletion(segment) => segment.range.0,
    });
    let mutations = mutations
      .iter()
      .map(|mutation| mutation_event_strings(&mutation.event))
      .collect::<Result<Vec<_>, _>>()?
      .into_iter()
      .flatten()
      .join(",");
    Ok(btreemap! {
      "mutations".to_owned() => mutations,
    })
  }
}

/// Optimize-specific read accessors, extending [`PartitionBranchOps`] with the per-edge likelihood
/// contribution and indel count the branch-length optimizer needs. Implemented by the same short-lived
/// per-representation read views.
pub trait PartitionOptimizeOps: PartitionBranchOps {
  /// Return the precomputed likelihood contribution for one edge.
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report>;

  /// Return the number of indel events on one edge.
  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize;
}
