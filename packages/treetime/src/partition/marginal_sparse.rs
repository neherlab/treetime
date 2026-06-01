use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::sample::{SampleMode, resolve_profile};
use crate::constants::SUPERTINY_NUMBER;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::{MutationCounts, is_profile_informative};
use crate::partition::marginal_passes;
use crate::partition::optimization_contribution::OptimizationContribution;
use crate::partition::sparse::{SparseEdgePartition, SparseNodePartition, SparseSeqDistribution, VarPos};
use crate::partition::traits::{
  BranchTopology, HasGtr, HasLogLh, PartitionBranchOps, PartitionMarginalOps, PartitionMarginalPasses,
  PartitionOptimizeOps, PartitionRerootOps, PartitionTimetreeOps, TransitionCounting,
};
use crate::seq::mutation::Sub;
use crate::{make_error, make_internal_report};
use eyre::Report;
use ndarray::{Array1, Array2};
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_graph::reroot::{EdgeMergeInfo, RerootChanges};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{Seq, seq};
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::collections::container::get_exactly_one;
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

impl HasGtr for PartitionMarginalSparse {
  fn gtr(&self) -> &GTR {
    &self.gtr
  }
  fn gtr_mut(&mut self) -> &mut GTR {
    &mut self.gtr
  }
  fn sequence_length(&self) -> usize {
    self.length
  }
}

pub(crate) fn reconstruct_map_seq(
  base_seq: &Seq,
  edge: Option<&SparseEdgePartition>,
  node: &SparseNodePartition,
  alphabet: &Alphabet,
) -> Seq {
  let mut rng = rand::thread_rng();
  reconstruct_map_seq_sampled(base_seq, edge, node, alphabet, false, &mut rng)
}

pub(crate) fn reconstruct_map_seq_sampled(
  base_seq: &Seq,
  edge: Option<&SparseEdgePartition>,
  node: &SparseNodePartition,
  alphabet: &Alphabet,
  sample: bool,
  rng: &mut dyn rand::RngCore,
) -> Seq {
  let mut seq = if let Some(edge) = edge {
    let mut seq = base_seq.clone();
    for m in edge.fitch_subs() {
      seq[m.pos()] = m.qry();
    }
    for indel in &edge.indels {
      if indel.deletion {
        seq[indel.range.0..indel.range.1].fill(alphabet.gap());
      } else {
        seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
      }
    }
    seq
  } else {
    base_seq.clone()
  };

  for r in &node.seq.unknown {
    seq[r.0..r.1].fill(alphabet.unknown());
  }

  for (pos, states) in &node.profile.variable {
    seq[*pos] = alphabet.char(resolve_profile(states.dis.view(), sample, rng));
  }

  if sample {
    let variable_positions = &node.profile.variable;
    let fixes: Vec<(usize, usize)> = seq
      .iter()
      .enumerate()
      .filter(|(pos, _)| !variable_positions.contains_key(pos))
      .filter_map(|(pos, &ch)| {
        node
          .profile
          .fixed
          .get(&ch)
          .map(|fp| (pos, resolve_profile(fp.view(), true, rng)))
      })
      .collect();
    for (pos, idx) in fixes {
      seq[pos] = alphabet.char(idx);
    }
  }

  seq
}

impl HasLogLh for PartitionMarginalSparse {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.nodes.get(&node_key).map_or(0.0, |node| node.profile.log_lh)
  }

  fn reset_node_log_likelihoods(&mut self) {
    for node_data in self.nodes.values_mut() {
      node_data.profile.log_lh = 0.0;
    }
  }
}

impl PartitionMarginalSparse {
  #[allow(clippy::same_name_method)]
  pub fn get_sequence_length(&self) -> usize {
    self.length
  }

  fn count_transitions_impl<N, E>(&self, graph: &Graph<N, E, ()>) -> Result<MutationCounts, Report>
  where
    N: GraphNode,
    E: EdgeOptimizeOps,
  {
    let n_states = self.gtr.pi.len();
    let min_bl = crate::constants::MIN_BRANCH_LENGTH_FRACTION / self.length as f64;
    let mut nij = Array2::zeros((n_states, n_states));
    let mut Ti = Array1::zeros(n_states);

    for edge in graph.get_edges() {
      let edge_arc = edge.read_arc();
      let branch_length = edge_arc.payload().read_arc().branch_length().unwrap_or(0.0).max(min_bl);
      let edge_key = edge_arc.key();
      let edge_data = &self.edges[&edge_key];

      let exp_qt = self.gtr.expQt(branch_length) + SUPERTINY_NUMBER;

      Self::accumulate_sparse_transitions(
        &edge_data.msg_to_child,
        &edge_data.msg_to_parent,
        &exp_qt,
        branch_length,
        n_states,
        &mut nij,
        &mut Ti,
      );
    }

    let root = graph.get_exactly_one_root()?;
    let root_key = root.read_arc().key();
    let root_profile = &self.nodes[&root_key].profile;
    let mut root_state = Array1::zeros(n_states);
    let root_dis = Self::aggregate_sparse_profile(root_profile, n_states);
    if is_profile_informative(root_dis.view(), n_states) {
      if let Some(root_idx) = argmax_first(&root_dis.view()) {
        root_state[root_idx] = 1.0;
      }
    }

    nij.diag_mut().fill(0.0);

    Ok(MutationCounts { nij, Ti, root_state })
  }

  fn accumulate_sparse_transitions(
    msg_to_child: &SparseSeqDistribution,
    msg_to_parent: &SparseSeqDistribution,
    exp_qt: &Array2<f64>,
    branch_length: f64,
    n_states: usize,
    nij: &mut Array2<f64>,
    Ti: &mut Array1<f64>,
  ) {
    for (pos, pp_var) in &msg_to_parent.variable {
      let pc = msg_to_child
        .variable
        .get(pos)
        .map_or_else(|| Self::fixed_profile_for_var(msg_to_child, pp_var), |v| &v.dis);
      Self::accumulate_site_transition_weighted(&pp_var.dis, pc, exp_qt, branch_length, n_states, nij, Ti, 1);
    }

    for (ch, parent_fixed) in &msg_to_parent.fixed {
      let count = msg_to_parent.fixed_counts.get(*ch).unwrap_or(0);
      if count == 0 {
        continue;
      }
      let child_fixed = msg_to_child.fixed.get(ch).unwrap_or(parent_fixed);
      Self::accumulate_site_transition_weighted(
        parent_fixed,
        child_fixed,
        exp_qt,
        branch_length,
        n_states,
        nij,
        Ti,
        count,
      );
    }
  }

  fn fixed_profile_for_var<'a>(dist: &'a SparseSeqDistribution, var: &'a VarPos) -> &'a Array1<f64> {
    dist.fixed.get(&var.state).unwrap_or(&var.dis)
  }

  fn accumulate_site_transition_weighted(
    pp: &Array1<f64>,
    pc: &Array1<f64>,
    exp_qt: &Array2<f64>,
    branch_length: f64,
    n_states: usize,
    nij: &mut Array2<f64>,
    Ti: &mut Array1<f64>,
    weight: usize,
  ) {
    let weight = weight as f64;
    let mut site_sum = 0.0;
    let mut joint = Array2::zeros((n_states, n_states));

    for i in 0..n_states {
      for j in 0..n_states {
        let val = pp[i] * exp_qt[[i, j]] * pc[j];
        joint[[i, j]] = val;
        site_sum += val;
      }
    }

    if site_sum > 0.0 {
      joint /= site_sum;
    }

    for i in 0..n_states {
      for j in 0..n_states {
        nij[[i, j]] += weight * joint[[i, j]];
      }
    }

    for k in 0..n_states {
      let mut parent_sum = 0.0;
      let mut child_sum = 0.0;
      for s in 0..n_states {
        parent_sum += joint[[s, k]];
        child_sum += joint[[k, s]];
      }
      Ti[k] += weight * 0.5 * branch_length * (parent_sum + child_sum);
    }
  }

  fn aggregate_sparse_profile(profile: &SparseSeqDistribution, n_states: usize) -> Array1<f64> {
    let mut result = Array1::zeros(n_states);
    for var in profile.variable.values() {
      if is_profile_informative(var.dis.view(), n_states) {
        if let Some(idx) = argmax_first(&var.dis.view()) {
          result[idx] += 1.0;
        }
      }
    }
    for (ch, fixed_profile) in &profile.fixed {
      let count = profile.fixed_counts.get(*ch).unwrap_or(0) as f64;
      if count > 0.0 && is_profile_informative(fixed_profile.view(), n_states) {
        if let Some(idx) = argmax_first(&fixed_profile.view()) {
          result[idx] += count;
        }
      }
    }
    result
  }

  // Phase 1: topology changes + root_sequence derivation + new node init.
  // Must complete before remove_trivial_root (phase 2) because
  // derive_root_sequence reads the parent_edge that phase 2 deletes.
  //
  // Note: root_composition + child-side fitch_subs + indels != child_composition.
  // Non-char (N, gap) differences between nodes are not encoded as Fitch subs -
  // they are tracked through non_char ranges on each node instead.
  fn apply_reroot_changes(&mut self, changes: &RerootChanges) -> Result<(), Report> {
    // Split edge: child-side gets all mutations, parent-side is empty
    if let Some(info) = &changes.edge_split {
      let mut old_edge_data = self
        .edges
        .remove(&info.old_edge_key)
        .ok_or_else(|| make_internal_report!("Old edge {:?} must exist for split", info.old_edge_key))?;
      old_edge_data.clear_ml_subs();
      self.edges.insert(info.child_side_edge_key, old_edge_data);
      self
        .edges
        .insert(info.parent_side_edge_key, SparseEdgePartition::default());
    }

    // Invert edges on the reroot path
    for edge_key in &changes.inverted_edge_keys {
      let edge_data = self
        .edges
        .get_mut(edge_key)
        .ok_or_else(|| make_internal_report!("Edge {edge_key:?} must exist on reroot path"))?;
      edge_data.invert_for_reroot();
    }

    // Derive root_sequence for the new root
    self.derive_root_sequence(changes);

    // Initialize the new split node from the finalized root_sequence
    if let Some(info) = &changes.edge_split {
      self.nodes.insert(
        info.new_node_key,
        SparseNodePartition::new(&self.root_sequence, &self.alphabet)?,
      );
    }

    Ok(())
  }

  fn remove_trivial_root(&mut self, info: &EdgeMergeInfo) -> Result<(), Report> {
    let parent_edge = self
      .edges
      .remove(&info.parent_edge_key)
      .ok_or_else(|| make_internal_report!("Parent edge {:?} must exist for merge", info.parent_edge_key))?;
    let child_edge = self
      .edges
      .remove(&info.child_edge_key)
      .ok_or_else(|| make_internal_report!("Child edge {:?} must exist for merge", info.child_edge_key))?;

    self.nodes.remove(&info.removed_node_key);

    let merged_subs = parent_edge.chain_fitch_subs(child_edge.fitch_subs())?;
    let merged_indels = parent_edge.chain_fitch_indels(&child_edge.indels);

    let mut merged_edge = SparseEdgePartition::default();
    merged_edge.set_fitch_subs(merged_subs);
    merged_edge.indels = merged_indels;

    self.edges.insert(info.merged_edge_key, merged_edge);
    Ok(())
  }

  fn derive_root_sequence(&mut self, changes: &RerootChanges) {
    if !changes.inverted_edge_keys.is_empty() {
      let mut new_root_seq = self.root_sequence.clone();
      for edge_key in &changes.inverted_edge_keys {
        if let Some(edge_data) = self.edges.get(edge_key) {
          Self::apply_edge_to_sequence(&mut new_root_seq, edge_data, &self.alphabet);
        }
      }
      self.root_sequence = new_root_seq;
    } else if let Some(info) = &changes.edge_merge {
      if let Some(parent_edge) = self.edges.get(&info.parent_edge_key) {
        let parent_edge = parent_edge.clone();
        Self::apply_edge_to_sequence(&mut self.root_sequence, &parent_edge, &self.alphabet);
      }
    }
  }

  fn apply_edge_to_sequence(seq: &mut Seq, edge: &SparseEdgePartition, alphabet: &Alphabet) {
    for sub in edge.fitch_subs() {
      if sub.pos() < seq.len() {
        seq[sub.pos()] = sub.reff();
      }
    }
    for indel in &edge.indels {
      if indel.range.0 < seq.len() && indel.range.1 <= seq.len() {
        if indel.deletion {
          seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
        } else {
          seq[indel.range.0..indel.range.1].fill(alphabet.gap());
        }
      }
    }
  }
}

impl<N, E> TransitionCounting<N, E> for PartitionMarginalSparse
where
  N: GraphNode,
  E: EdgeOptimizeOps,
{
  fn count_transitions(&self, graph: &Graph<N, E, ()>) -> Result<MutationCounts, Report> {
    self.count_transitions_impl(graph)
  }
}

impl PartitionRerootOps for PartitionMarginalSparse {
  fn apply_reroot(&mut self, changes: &RerootChanges) -> Result<(), Report> {
    self.apply_reroot_changes(changes)?;

    if let Some(info) = &changes.edge_merge {
      self.remove_trivial_root(info)?;
    }

    Ok(())
  }
}

impl PartitionBranchOps for PartitionMarginalSparse {
  fn sequence_length(&self) -> usize {
    self.length
  }

  fn edge_subs(&self, _graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    let edge = &self.edges[&edge_key];
    match edge.ml_subs() {
      Some(subs) => Ok(subs.to_vec()),
      None => make_error!("edge_subs() called before marginal inference populated subs_ml for edge {edge_key:?}"),
    }
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
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    OptimizationContribution::from_sparse(edge_key, self)
  }

  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.edges[&edge_key].indels.len()
  }
}

impl<N, E> PartitionTimetreeOps<N, E> for PartitionMarginalSparse
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn reconcile_topology(&mut self, graph: &Graph<N, E, ()>) {
    let graph_node_keys: BTreeSet<GraphNodeKey> = graph.get_nodes().into_iter().map(|n| n.read_arc().key()).collect();
    let graph_edge_keys: BTreeSet<GraphEdgeKey> = graph.get_edges().into_iter().map(|e| e.read_arc().key()).collect();

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

impl<N, E> PartitionMarginalPasses<N, E> for PartitionMarginalSparse
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn process_node_backward(&mut self, node: &GraphNodeBackward<N, E, ()>) -> Result<(), Report> {
    marginal_passes::process_node_backward(self, node)
  }

  fn process_node_forward(&mut self, graph: &Graph<N, E, ()>, node: &GraphNodeForward<N, E, ()>) -> Result<(), Report> {
    marginal_passes::process_node_forward(self, graph, node)
  }

  fn get_sequence_length(&self) -> usize {
    self.length
  }
}

impl<N, E> PartitionMarginalOps<N, E> for PartitionMarginalSparse
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn attach_sequences(&mut self, _graph: &Graph<N, E, ()>, _aln: &[FastaRecord]) -> Result<(), Report> {
    Ok(())
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
    &mut self,
    node: &GraphNodeForward<N, E, ()>,
    include_leaves: bool,
    sample_mode: SampleMode,
    rng: &mut dyn rand::RngCore,
  ) -> Option<Seq> {
    if !include_leaves && node.is_leaf {
      return None;
    }

    let mut node_data = self.nodes.remove(&node.key)?;

    let (base_seq, edge) = if node.is_root {
      (&self.root_sequence, None)
    } else {
      let (parent_key, edge_key) = get_exactly_one(&node.parent_keys).ok()?;
      let parent_data = self.nodes.get(parent_key)?;
      let edge_data = self.edges.get(edge_key)?;
      (&parent_data.seq.sequence, Some(edge_data))
    };

    let sample = sample_mode.samples_node(node.is_root);
    let seq = reconstruct_map_seq_sampled(base_seq, edge, &node_data, &self.alphabet, sample, rng);
    node_data.seq.sequence = seq.clone();
    self.nodes.insert(node.key, node_data);

    Some(seq)
  }
}
