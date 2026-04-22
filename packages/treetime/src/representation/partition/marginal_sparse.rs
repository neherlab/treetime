use crate::alphabet::alphabet::Alphabet;
use crate::commands::clock::reroot::RerootChanges;
use crate::commands::optimize::partition_ops::PartitionOptimizeOps;
use crate::commands::timetree::partition_ops::PartitionRerootOps;
use crate::gtr::gtr::GTR;
use crate::make_internal_report;
use crate::representation::partition::marginal_passes;
use crate::representation::partition::traits::BranchTopology;
use crate::representation::partition::traits::HasLogLh;
use crate::representation::partition::traits::PartitionBranchOps;
use crate::representation::partition::traits::PartitionCompressed;
use crate::representation::partition::traits::{PartitionMarginal, PartitionMarginalOps};
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
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range::range_contains;
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug)]
pub struct PartitionMarginalSparse {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub root_sequence: Seq,
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

  /// Extract root sequence from the root node and clear internal node sequences.
  ///
  /// After Fitch compression, every node carries a full resolved sequence. Only
  /// the root sequence is needed for downstream reconstruction (node_state_at,
  /// reconstruct_node_sequence cascade from it via edge subs). Clearing internal
  /// sequences saves memory and removes stale data that would otherwise persist
  /// across reroots.
  pub fn extract_root_sequence<N, E>(&mut self, graph: &Graph<N, E, ()>) -> Result<(), Report>
  where
    N: GraphNode,
    E: GraphEdge,
  {
    let root_key = graph.get_exactly_one_root()?.read_arc().key();
    self.root_sequence = self.nodes[&root_key].seq.sequence.clone();
    for (key, node_data) in &mut self.nodes {
      if *key != root_key && !graph.is_leaf(*key) {
        node_data.seq.sequence = seq![];
      }
    }
    Ok(())
  }

  /// Return positions that can change the reconstructed mutation set for one edge.
  ///
  /// In sparse mode, only a small set of sites can change the branch result:
  /// sites changed on this edge, or sites where the parent or child has a
  /// non-default marginal state. Everything else stays the same on both ends.
  fn edge_candidate_positions(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<usize>, Report> {
    let Some(edge) = self.edges.get(&edge_key) else {
      return Ok(vec![]);
    };
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;

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

  /// Reconstruct one node state at one site from sparse data.
  ///
  /// This is the single-site version of `reconstruct_node_sequence()`: start
  /// from the parent state, apply edge changes, mask unknowns, then apply the
  /// node's current best state for variable sites.
  ///
  /// After marginal inference, `edge.subs` alone is not enough. The final state
  /// also depends on the node's current marginal result.
  fn node_state_at(
    &self,
    graph: &dyn BranchTopology,
    node_key: GraphNodeKey,
    pos: usize,
    cache: &mut BTreeMap<(GraphNodeKey, usize), AsciiChar>,
  ) -> Result<AsciiChar, Report> {
    if let Some(state) = cache.get(&(node_key, pos)) {
      return Ok(*state);
    }

    let node_data = self.nodes.get(&node_key);

    debug_assert!(
      !self.root_sequence.is_empty(),
      "root_sequence is empty: call extract_root_sequence() after compress_sequences()"
    );

    let base_state = match graph.node_parent(node_key)? {
      None => self
        .root_sequence
        .get(pos)
        .copied()
        .unwrap_or_else(|| self.alphabet.char(0)),
      Some((parent_key, parent_edge_key)) => {
        let parent_state = self.node_state_at(graph, parent_key, pos, cache)?;

        // Apply edge changes if edge partition data exists
        self.edges.get(&parent_edge_key).map_or(parent_state, |edge_data| {
          // Indels have highest precedence, then substitutions, then parent state
          if let Some(indel) = edge_data
            .indels
            .iter()
            .find(|indel| indel.range.0 <= pos && pos < indel.range.1)
          {
            if indel.deletion {
              self.alphabet.gap()
            } else {
              indel.seq[pos - indel.range.0]
            }
          } else if let Some(sub) = edge_data.subs.iter().find(|sub| sub.pos() == pos) {
            sub.qry()
          } else {
            parent_state
          }
        })
      },
    };

    // Variable sites have highest precedence, then unknown, then base state
    let state = node_data.map_or(base_state, |d| {
      if let Some(states) = d.profile.variable.get(&pos) {
        self.alphabet.char(argmax_first(&states.dis.view()).unwrap_or(0))
      } else if range_contains(&d.seq.unknown, pos) {
        self.alphabet.unknown()
      } else {
        base_state
      }
    });

    cache.insert((node_key, pos), state);
    Ok(state)
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

    // Handle edge merge: compose mutations from two edges.
    // When the removed node is the old root (trivial root removal), update
    // root_sequence using the parent_edge subs before discarding them.
    // The parent_edge goes from the new-root side toward the removed node:
    // sub.reff() = new-root-side state, sub.qry() = removed-node state.
    if let Some(info) = &changes.edge_merge {
      let parent_edge = self
        .edges
        .remove(&info.parent_edge_key)
        .ok_or_else(|| make_internal_report!("Parent edge {:?} must exist for merge", info.parent_edge_key))?;
      let child_edge = self
        .edges
        .remove(&info.child_edge_key)
        .ok_or_else(|| make_internal_report!("Child edge {:?} must exist for merge", info.child_edge_key))?;

      if changes.inverted_edge_keys.is_empty() {
        for sub in &parent_edge.subs {
          if sub.pos() < self.root_sequence.len() {
            self.root_sequence[sub.pos()] = sub.reff();
          }
        }
        for indel in &parent_edge.indels {
          if indel.range.0 < self.root_sequence.len() && indel.range.1 <= self.root_sequence.len() {
            if indel.deletion {
              self.root_sequence[indel.range.0..indel.range.1].fill(self.alphabet.gap());
            } else {
              self.root_sequence[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
            }
          }
        }
      }

      self.nodes.remove(&info.removed_node_key);

      let merged_subs = compose_substitutions(&parent_edge.subs, &child_edge.subs)?;

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

    // Derive root sequence for the new root from the old root sequence.
    // Walk inverted edges (listed in old_root->new_root order). After sub
    // inversion, each sub's ref() holds the new-root-ward state (was the
    // original child qry). Applying ref() at each position advances the
    // sequence one step toward the new root.
    if !changes.inverted_edge_keys.is_empty() {
      let mut new_root_seq = self.root_sequence.clone();
      for edge_key in &changes.inverted_edge_keys {
        if let Some(edge_data) = self.edges.get(edge_key) {
          for sub in &edge_data.subs {
            if sub.pos() < new_root_seq.len() {
              new_root_seq[sub.pos()] = sub.reff();
            }
          }
          for indel in &edge_data.indels {
            if indel.range.0 < new_root_seq.len() && indel.range.1 <= new_root_seq.len() {
              // After inversion, deletion flag is flipped. To go in the
              // original (old_root->new_root) direction:
              // - current deletion (was insertion): child had `seq`
              // - current insertion (was deletion): child had gaps
              if indel.deletion {
                new_root_seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
              } else {
                new_root_seq[indel.range.0..indel.range.1].fill(self.alphabet.gap());
              }
            }
          }
        }
      }
      self.root_sequence = new_root_seq;
    }

    Ok(())
  }
}

impl PartitionBranchOps for PartitionMarginalSparse {
  fn sequence_length(&self) -> usize {
    self.length
  }

  fn edge_subs(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    let mut cache = BTreeMap::new();
    let mut subs = Vec::new();

    for pos in self.edge_candidate_positions(graph, edge_key)? {
      let parent_state = self.node_state_at(graph, parent_key, pos, &mut cache)?;
      let child_state = self.node_state_at(graph, child_key, pos, &mut cache)?;

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

  fn edge_effective_length(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
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
    edge_key: GraphEdgeKey,
  ) -> Result<crate::commands::optimize::optimize_unified::OptimizationContribution, Report> {
    crate::commands::optimize::optimize_unified::OptimizationContribution::from_sparse(edge_key, self)
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

  fn process_node_backward(&mut self, node: &GraphNodeBackward<N, E, ()>) -> Result<(), Report> {
    marginal_passes::process_node_backward(self, node)
  }

  fn process_node_forward(&mut self, graph: &Graph<N, E, ()>, node: &GraphNodeForward<N, E, ()>) -> Result<(), Report> {
    marginal_passes::process_node_forward(self, graph, node)
  }

  fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    if let Some(node_data) = self.nodes.get(&node_key) {
      if !node_data.seq.sequence.is_empty() {
        node_data.seq.sequence.clone()
      } else {
        // Return empty seq to be filled later by the reconstruction algorithm
        seq![]
      }
    } else {
      seq![]
    }
  }

  fn reconstruct_node_sequence(&mut self, node: &GraphNodeForward<N, E, ()>, include_leaves: bool) -> Option<Seq> {
    if !include_leaves && node.is_leaf {
      return None;
    }

    let mut node_data = self.nodes.remove(&node.key)?;

    let mut seq = if node.is_root {
      self.root_sequence.clone()
    } else {
      let (parent_key, edge_key) = get_exactly_one(&node.parent_keys).ok()?;
      let parent_data = self.nodes.get(parent_key)?;
      let edge_data = self.edges.get(edge_key)?;

      let mut seq = parent_data.seq.sequence.clone();

      // Implant mutations
      for m in &edge_data.subs {
        seq[m.pos()] = m.qry();
      }

      // Implant indels
      for indel in &edge_data.indels {
        if indel.deletion {
          seq[indel.range.0..indel.range.1].fill(self.alphabet.gap());
        } else {
          seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
        }
      }

      seq
    };

    // At the node itself, mask whatever is unknown in the node.
    let alphabet = &self.alphabet;
    for r in &node_data.seq.unknown {
      let ambig_char = alphabet.unknown();
      seq[r.0..r.1].fill(ambig_char);
    }

    // Change variable sites to their most likely state.
    // Uses argmax_first for NumPy-compatible tie-breaking.
    for (pos, states) in &node_data.profile.variable {
      seq[*pos] = alphabet.char(argmax_first(&states.dis.view()).unwrap_or(0));
    }

    node_data.seq.sequence = seq.clone();
    self.nodes.insert(node.key, node_data);

    Some(seq)
  }

  fn get_sequence_length(&self) -> usize {
    self.length
  }
}
