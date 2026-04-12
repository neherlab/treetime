use crate::alphabet::alphabet::Alphabet;
use crate::commands::clock::reroot::RerootChanges;
use crate::commands::optimize::partition_ops::PartitionOptimizeOps;
use crate::commands::timetree::partition_ops::PartitionRerootOps;
use crate::gtr::gtr::GTR;
use crate::make_internal_report;
use crate::representation::partition::marginal_passes;
use crate::representation::partition::traits::HasLogLh;
use crate::representation::partition::traits::PartitionBranchOps;
use crate::representation::partition::traits::PartitionCompressed;
use crate::representation::partition::traits::{
  CurrentStateCache, ExactStateCache, GraphNodePathLookup, PartitionMarginal, PartitionMarginalOps,
  SequenceReconstructionCache,
};
use crate::representation::payload::ancestral::GraphAncestral;
use crate::representation::payload::sparse::{MarginalSparseSeqDistribution, SparseEdgePartition, SparseNodePartition};
use crate::seq::mutation::{Sub, compose_substitutions};
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::mem;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdge, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{AsciiChar, Seq, seq};
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::interval::range::range_contains;
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug)]
pub struct PartitionMarginalSparse {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, SparseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, SparseEdgePartition>,
}

impl PartitionCompressed for PartitionMarginalSparse {
  fn index(&self) -> usize {
    self.index
  }

  fn alphabet(&self) -> &Alphabet {
    &self.alphabet
  }

  fn length(&self) -> usize {
    self.length
  }

  fn nodes(&self) -> &BTreeMap<GraphNodeKey, SparseNodePartition> {
    &self.nodes
  }

  fn edges(&self) -> &BTreeMap<GraphEdgeKey, SparseEdgePartition> {
    &self.edges
  }

  fn nodes_mut(&mut self) -> &mut BTreeMap<GraphNodeKey, SparseNodePartition> {
    &mut self.nodes
  }

  fn edges_mut(&mut self) -> &mut BTreeMap<GraphEdgeKey, SparseEdgePartition> {
    &mut self.edges
  }
}

impl PartitionMarginalSparse {
  #[allow(clippy::same_name_method)]
  pub fn get_sequence_length(&self) -> usize {
    self.length
  }

  /// Return current nucleotide changes for one edge, generic over graph type.
  ///
  /// This is the generic implementation used by both `PartitionBranchOps::edge_subs()`
  /// (which constrains to `GraphAncestral`) and direct callers that have other graph
  /// types (e.g. GTR inference in the timetree context with `Graph<NodeTimetree, ...>`).
  pub fn edge_subs_from_graph<N: GraphNode, E: GraphEdge, D: Send + Sync>(
    &self,
    graph: &Graph<N, E, D>,
    edge_key: GraphEdgeKey,
  ) -> Result<Vec<Sub>, Report> {
    let parent_key = graph.get_source_node_key(edge_key)?;
    let child_key = graph.get_target_node_key(edge_key)?;
    let mut exact_cache = ExactStateCache::new();
    let mut current_cache = CurrentStateCache::new();
    let mut subs = Vec::new();

    for pos in self.edge_candidate_positions(graph, edge_key)? {
      let parent_state = self.node_state_at(graph, parent_key, pos, &mut exact_cache, &mut current_cache)?;
      let child_state = self.node_state_at(graph, child_key, pos, &mut exact_cache, &mut current_cache)?;

      if parent_state == child_state {
        continue;
      }

      // Gaps and unknowns are not nucleotide substitutions.
      if !self.alphabet.is_canonical(parent_state) || !self.alphabet.is_canonical(child_state) {
        continue;
      }

      subs.push(Sub::new(parent_state, pos, child_state)?);
    }

    Ok(subs)
  }

  /// Return positions that can change the reconstructed mutation set for one edge.
  ///
  /// In sparse mode, only a small set of sites can change the branch result:
  /// sites changed on this edge, or sites where the parent or child has a
  /// non-default marginal state. Everything else stays the same on both ends.
  fn edge_candidate_positions<N: GraphNode, E: GraphEdge, D: Send + Sync>(
    &self,
    graph: &Graph<N, E, D>,
    edge_key: GraphEdgeKey,
  ) -> Result<Vec<usize>, Report> {
    let Some(edge) = self.edges.get(&edge_key) else {
      return Ok(vec![]);
    };
    let parent_key = graph.get_source_node_key(edge_key)?;
    let child_key = graph.get_target_node_key(edge_key)?;

    // Nodes may be absent when called during topology changes (e.g. after
    // merge_sibling_pair creates a new node before marginal reconstruction
    // populates its profile). Missing nodes contribute no variable positions.
    let empty_iter = || -> Box<dyn Iterator<Item = usize>> { Box::new(std::iter::empty()) };
    let parent_vars: Box<dyn Iterator<Item = usize>> = self
      .nodes
      .get(&parent_key)
      .map_or_else(empty_iter, |n| Box::new(n.profile.variable.keys().copied()));
    let child_vars: Box<dyn Iterator<Item = usize>> = self
      .nodes
      .get(&child_key)
      .map_or_else(empty_iter, |n| Box::new(n.profile.variable.keys().copied()));

    Ok(
      edge
        .subs
        .iter()
        .map(Sub::pos)
        .chain(parent_vars)
        .chain(child_vars)
        // Stable ordering keeps `Vec<Sub>` equality meaningful in tests and
        // callers that compare branch mutation lists directly.
        .sorted()
        .dedup()
        .collect_vec(),
    )
  }

  /// Reconstruct the exact Fitch-resolved state at one node and site.
  ///
  /// This follows the canonical root sequence and the edge deltas only. It does
  /// not consult the marginal MAP profile.
  fn edge_state_from_parent(&self, edge_key: GraphEdgeKey, pos: usize, parent_state: AsciiChar) -> AsciiChar {
    self.edges.get(&edge_key).map_or(parent_state, |edge_data| {
      match edge_data
        .indels
        .iter()
        .find(|indel| indel.range.0 <= pos && pos < indel.range.1)
      {
        Some(indel) if indel.deletion => self.alphabet.gap(),
        Some(indel) => indel.seq[pos - indel.range.0],
        None => edge_data
          .subs
          .iter()
          .find(|sub| sub.pos() == pos)
          .map_or(parent_state, Sub::qry),
      }
    })
  }

  pub(crate) fn node_canonical_state_known_parent<G: GraphNodePathLookup + ?Sized>(
    &self,
    graph: &G,
    node_key: GraphNodeKey,
    parent: Option<(GraphNodeKey, GraphEdgeKey)>,
    pos: usize,
    cache: &mut ExactStateCache,
  ) -> Result<AsciiChar, Report> {
    if let Some(state) = cache.get(&(node_key, pos)).copied() {
      return Ok(state);
    }

    let node_data = self
      .nodes
      .get(&node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} has no partition data"))?;
    let state = if range_contains(&node_data.seq.gaps, pos) {
      self.alphabet.gap()
    } else if range_contains(&node_data.seq.unknown, pos) {
      self.alphabet.unknown()
    } else if let Some((parent_key, edge_key)) = parent {
      let parent_state = self.node_canonical_state(graph, parent_key, pos, cache)?;
      self.edge_state_from_parent(edge_key, pos, parent_state)
    } else {
      node_data
        .seq
        .sequence
        .get(pos)
        .copied()
        .unwrap_or_else(|| self.alphabet.char(0))
    };

    cache.insert((node_key, pos), state);
    Ok(state)
  }

  pub(crate) fn child_canonical_state_from_parent_state(
    &self,
    child_key: GraphNodeKey,
    edge_key: GraphEdgeKey,
    parent_state: AsciiChar,
    pos: usize,
    cache: &mut ExactStateCache,
  ) -> Result<AsciiChar, Report> {
    if let Some(state) = cache.get(&(child_key, pos)).copied() {
      return Ok(state);
    }

    let child_data = self
      .nodes
      .get(&child_key)
      .ok_or_else(|| make_internal_report!("Node {child_key} has no partition data"))?;
    let state = if range_contains(&child_data.seq.gaps, pos) {
      self.alphabet.gap()
    } else if range_contains(&child_data.seq.unknown, pos) {
      self.alphabet.unknown()
    } else {
      self.edge_state_from_parent(edge_key, pos, parent_state)
    };

    cache.insert((child_key, pos), state);
    Ok(state)
  }

  pub(crate) fn node_canonical_state<G: GraphNodePathLookup + ?Sized>(
    &self,
    graph: &G,
    node_key: GraphNodeKey,
    pos: usize,
    cache: &mut ExactStateCache,
  ) -> Result<AsciiChar, Report> {
    let parent = graph.parent_of(node_key)?;
    self.node_canonical_state_known_parent(graph, node_key, parent, pos, cache)
  }

  fn leaf_parsimony_state<G: GraphNodePathLookup + ?Sized>(
    &self,
    graph: &G,
    node_key: GraphNodeKey,
    pos: usize,
    exact_cache: &mut ExactStateCache,
  ) -> Result<AsciiChar, Report> {
    let node_data = self
      .nodes
      .get(&node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} has no partition data"))?;
    if range_contains(&node_data.seq.gaps, pos) {
      return Ok(self.alphabet.gap());
    }
    if range_contains(&node_data.seq.unknown, pos) {
      return Ok(self.alphabet.unknown());
    }

    // Sparse leaf semantics follow the exact parsimony-resolved state carried by
    // the compressed tree, not the raw observed ambiguity code from the input
    // FASTA. The ambiguity itself still shapes the leaf likelihood through the
    // backward `dis` vector; downstream state-based consumers need the
    // parsimony-resolved nucleotide that anchors those messages.
    self.node_canonical_state(graph, node_key, pos, exact_cache)
  }

  fn reconstruct_leaf_parsimony_sequence<G: GraphNodePathLookup + ?Sized>(
    &self,
    graph: &G,
    node_key: GraphNodeKey,
    exact_cache: &mut ExactStateCache,
  ) -> Result<Seq, Report> {
    (0..self.length)
      .map(|pos| self.leaf_parsimony_state(graph, node_key, pos, exact_cache))
      .collect()
  }

  /// Reconstruct one node state at one site from sparse data.
  ///
  /// This starts from the exact Fitch-resolved state and then applies the
  /// current marginal MAP profile when the position is explicit in the sparse
  /// posterior.
  fn node_state_at<G: GraphNodePathLookup + ?Sized>(
    &self,
    graph: &G,
    node_key: GraphNodeKey,
    pos: usize,
    exact_cache: &mut ExactStateCache,
    current_cache: &mut CurrentStateCache,
  ) -> Result<AsciiChar, Report> {
    if let Some(state) = current_cache.get(&(node_key, pos)).copied() {
      return Ok(state);
    }

    let node_data = self
      .nodes
      .get(&node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} has no partition data"))?;
    if graph.is_leaf(node_key)? {
      let state = self.leaf_parsimony_state(graph, node_key, pos, exact_cache)?;
      current_cache.insert((node_key, pos), state);
      return Ok(state);
    }
    let base_state = if let Some((parent_key, edge_key)) = graph.parent_of(node_key)? {
      let parent_state = self.node_canonical_state(graph, parent_key, pos, exact_cache)?;
      self.child_canonical_state_from_parent_state(node_key, edge_key, parent_state, pos, exact_cache)?
    } else {
      self.node_canonical_state_known_parent(graph, node_key, None, pos, exact_cache)?
    };
    let state = node_data.profile.variable.get(&pos).map_or(base_state, |states| {
      self.alphabet.char(argmax_first(&states.dis.view()).unwrap_or(0))
    });
    current_cache.insert((node_key, pos), state);
    Ok(state)
  }

  fn reconstruct_sequence_cached<G: GraphNodePathLookup + ?Sized>(
    &self,
    graph: &G,
    node_key: GraphNodeKey,
    cache: &mut SequenceReconstructionCache,
  ) -> Result<Seq, Report> {
    if let Some(seq) = cache.get(&node_key) {
      return Ok(seq.clone());
    }

    let node_data = self
      .nodes
      .get(&node_key)
      .ok_or_else(|| make_internal_report!("Node {node_key} has no partition data"))?;
    if graph.is_leaf(node_key)? {
      let mut exact_cache = ExactStateCache::default();
      let seq = self.reconstruct_leaf_parsimony_sequence(graph, node_key, &mut exact_cache)?;
      cache.insert(node_key, seq.clone());
      return Ok(seq);
    }
    let mut seq = if let Some((parent_key, edge_key)) = graph.parent_of(node_key)? {
      let mut seq = self.reconstruct_sequence_cached(graph, parent_key, cache)?;
      if let Some(edge_data) = self.edges.get(&edge_key) {
        for sub in &edge_data.subs {
          seq[sub.pos()] = sub.qry();
        }
        for indel in &edge_data.indels {
          if indel.deletion {
            seq[indel.range.0..indel.range.1].fill(self.alphabet.gap());
          } else {
            seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
          }
        }
      }
      seq
    } else {
      node_data.seq.sequence.clone()
    };

    for gap in &node_data.seq.gaps {
      seq[gap.0..gap.1].fill(self.alphabet.gap());
    }
    for unknown in &node_data.seq.unknown {
      seq[unknown.0..unknown.1].fill(self.alphabet.unknown());
    }
    for (&pos, states) in &node_data.profile.variable {
      seq[pos] = self.alphabet.char(argmax_first(&states.dis.view()).unwrap_or(0));
    }
    cache.insert(node_key, seq.clone());
    Ok(seq)
  }
}

impl HasLogLh for PartitionMarginalSparse {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.nodes.get(&node_key).map_or(0.0, |node| node.profile.log_lh)
  }
}

impl PartitionMarginal for PartitionMarginalSparse {}

impl PartitionRerootOps for PartitionMarginalSparse {
  fn apply_reroot(&mut self, changes: &RerootChanges) -> Result<(), Report> {
    // Handle edge split: create new node entry, move mutations to child-side edge
    if let Some(info) = &changes.edge_split {
      self
        .nodes
        .insert(info.new_node_key, SparseNodePartition::empty(&self.alphabet));

      let old_edge_data = self
        .edges
        .remove(&info.old_edge_key)
        .ok_or_else(|| make_internal_report!("Old edge {:?} must exist for split", info.old_edge_key))?;

      // Child-side edge gets all mutations (they describe parent->child relationship)
      self.edges.insert(info.child_side_edge_key, old_edge_data);

      // Parent-side edge is empty (no mutations between parent and new split node)
      self
        .edges
        .insert(info.parent_side_edge_key, SparseEdgePartition::default());
    }

    // Handle edge merge: compose mutations from two edges
    if let Some(info) = &changes.edge_merge {
      self.nodes.remove(&info.removed_node_key);

      let parent_edge = self
        .edges
        .remove(&info.parent_edge_key)
        .ok_or_else(|| make_internal_report!("Parent edge {:?} must exist for merge", info.parent_edge_key))?;
      let child_edge = self
        .edges
        .remove(&info.child_edge_key)
        .ok_or_else(|| make_internal_report!("Child edge {:?} must exist for merge", info.child_edge_key))?;

      // Compose substitutions: parent then child
      let merged_subs = compose_substitutions(&parent_edge.subs, &child_edge.subs)?;

      // Compose indels: concatenate (parent indels first, then child indels)
      let mut merged_indels = parent_edge.indels;
      merged_indels.extend(child_edge.indels);

      let merged_edge = SparseEdgePartition {
        subs: merged_subs,
        indels: merged_indels,
        ..SparseEdgePartition::default()
      };

      self.edges.insert(info.merged_edge_key, merged_edge);
    }

    // Invert edge data for each edge on the reroot path
    for edge_key in &changes.inverted_edge_keys {
      let edge_data = self
        .edges
        .get_mut(edge_key)
        .ok_or_else(|| make_internal_report!("Edge {edge_key:?} must exist on reroot path"))?;

      // Invert all substitutions
      for sub in &mut edge_data.subs {
        sub.invert();
      }

      // Invert all indels
      for indel in &mut edge_data.indels {
        indel.invert();
      }

      // Swap directional messages
      mem::swap(&mut edge_data.msg_to_parent, &mut edge_data.msg_to_child);

      // Reset msg_from_child (stale propagated message cache)
      edge_data.msg_from_child = MarginalSparseSeqDistribution::default();
    }

    Ok(())
  }
}

impl PartitionBranchOps for PartitionMarginalSparse {
  fn sequence_length(&self) -> usize {
    self.length
  }

  fn edge_subs(&self, graph: &GraphAncestral, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    self.edge_subs_from_graph(graph, edge_key)
  }

  fn edge_effective_length(&self, graph: &GraphAncestral, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    let parent_key = graph.get_source_node_key(edge_key)?;
    let child_key = graph.get_target_node_key(edge_key)?;
    let parent_non_char = &self.nodes[&parent_key].seq.non_char;
    let child_non_char = &self.nodes[&child_key].seq.non_char;

    // non_char covers both gaps and unknowns (positions that do not evolve
    // under the substitution model). For internal nodes, non_char is the
    // intersection of children's non_char (Fitch backward pass), so a
    // position is excluded only when all descendants lack data there.
    let non_char_positions: usize = range_union(&[parent_non_char.clone(), child_non_char.clone()])
      .iter()
      .map(|(start, end)| end - start)
      .sum();

    Ok(self.length.saturating_sub(non_char_positions))
  }
}

impl PartitionOptimizeOps for PartitionMarginalSparse {
  fn create_edge_contribution(
    &self,
    graph: &dyn GraphNodePathLookup,
    edge_key: GraphEdgeKey,
    cache: &mut ExactStateCache,
  ) -> Result<crate::commands::optimize::optimize_unified::OptimizationContribution, Report> {
    crate::commands::optimize::optimize_unified::OptimizationContribution::from_sparse(graph, edge_key, self, cache)
  }

  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.edges[&edge_key].indels.len()
  }
}

impl<N, E> crate::commands::timetree::partition_ops::PartitionTimetreeOps<N, E> for PartitionMarginalSparse
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn reconcile_topology(&mut self, graph: &Graph<N, E, ()>) {
    let graph_node_keys: std::collections::BTreeSet<GraphNodeKey> =
      graph.get_nodes().into_iter().map(|n| n.read_arc().key()).collect();
    let graph_edge_keys: std::collections::BTreeSet<GraphEdgeKey> =
      graph.get_edges().into_iter().map(|e| e.read_arc().key()).collect();

    // Add missing nodes with empty partition data
    for &key in &graph_node_keys {
      self
        .nodes
        .entry(key)
        .or_insert_with(|| SparseNodePartition::empty(&self.alphabet));
    }

    // Add missing edges with default partition data
    for &key in &graph_edge_keys {
      self.edges.entry(key).or_default();
    }

    // Remove stale entries for nodes/edges no longer in graph
    self.nodes.retain(|k, _| graph_node_keys.contains(k));
    self.edges.retain(|k, _| graph_edge_keys.contains(k));
  }
}

impl<N, E> PartitionMarginalOps<N, E> for PartitionMarginalSparse
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn attach_sequences(&mut self, _graph: &Graph<N, E, ()>, _aln: &[FastaRecord]) -> Result<(), Report> {
    // Sparse partitions get sequences attached during compression phase
    Ok(())
  }

  fn process_node_backward(
    &mut self,
    graph: &Graph<N, E, ()>,
    node: &GraphNodeBackward<N, E, ()>,
    cache: &mut ExactStateCache,
  ) -> Result<(), Report> {
    marginal_passes::process_node_backward(self, graph, node, cache)
  }

  fn process_node_forward(
    &mut self,
    graph: &Graph<N, E, ()>,
    node: &GraphNodeForward<N, E, ()>,
    cache: &mut ExactStateCache,
  ) -> Result<(), Report> {
    marginal_passes::process_node_forward(self, graph, node, cache)
  }

  fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    if let Some(node_data) = self.nodes.get(&node_key) {
      if !node_data.seq.sequence.is_empty() {
        node_data.seq.sequence.clone()
      } else {
        seq![]
      }
    } else {
      seq![]
    }
  }

  fn reconstruct_node_sequence(
    &self,
    graph: &Graph<N, E, ()>,
    node: &GraphNodeForward<N, E, ()>,
    include_leaves: bool,
    cache: &mut SequenceReconstructionCache,
  ) -> Option<Seq> {
    if !include_leaves && node.is_leaf {
      return None;
    }

    self.reconstruct_sequence_cached(graph, node.key, cache).ok()
  }

  fn get_sequence_length(&self) -> usize {
    self.length
  }
}
