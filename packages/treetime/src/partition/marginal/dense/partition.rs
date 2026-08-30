use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::sample::{SampleMode, resolve_profile};
use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::make_report;
use crate::partition::marginal::shared::MarginalData;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDistribution, DenseSeqInfo};
use crate::partition::traits::{
  BranchTopology, HasGtr, HasLogLh, PartitionBranchOps, PartitionMarginalOps, PartitionOptimizeOps, PartitionRerootOps,
  PartitionTimetreeOps, TransitionCounting,
};
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::izip;
use maplit::btreemap;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named, NodeAncestralOps};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{LogLh, Seq, seq};
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::interval::range::range_contains;
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug, Serialize)]
pub struct PartitionMarginalDense {
  pub data: MarginalData,
  pub index: usize,
  pub alphabet: Alphabet,
  pub length: usize,
}

impl PartitionMarginalDense {
  pub fn new(index: usize, gtr: GTR, alphabet: Alphabet, length: usize) -> Self {
    let min_branch_length = MIN_BRANCH_LENGTH_FRACTION / length as f64;
    Self {
      data: MarginalData {
        gtr,
        nodes: btreemap! {},
        edges: btreemap! {},
        min_branch_length,
        // Nucleotide ancestral inference filters signal-free (gap-only) root
        // columns out of the equilibrium-frequency prior.
        filter_uninformative_root: true,
      },
      index,
      alphabet,
      length,
    }
  }

  #[allow(clippy::same_name_method)]
  pub fn get_sequence_length(&self) -> usize {
    self.length
  }
}

impl HasGtr for PartitionMarginalDense {
  fn gtr(&self) -> &GTR {
    &self.data.gtr
  }
  fn gtr_mut(&mut self) -> &mut GTR {
    &mut self.data.gtr
  }
  fn sequence_length(&self) -> usize {
    self.length
  }
}

impl HasLogLh for PartitionMarginalDense {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> LogLh {
    self
      .data
      .nodes
      .get(&node_key)
      .map_or(LogLh::ZERO, |node| node.profile.log_lh)
  }

  fn reset_node_log_likelihoods(&mut self) {
    for node_data in self.data.nodes.values_mut() {
      node_data.profile.log_lh = LogLh::ZERO;
    }
  }
}

impl PartitionRerootOps for PartitionMarginalDense {
  fn apply_reroot(&mut self, changes: &treetime_graph::reroot::RerootChanges) -> Result<(), Report> {
    if let Some(info) = &changes.edge_split {
      self.data.edges.remove(&info.old_edge_key);
      self.data.edges.entry(info.parent_side_edge_key).or_default();
      self.data.edges.entry(info.child_side_edge_key).or_default();
      self
        .data
        .nodes
        .entry(info.new_node_key)
        .or_insert_with(|| DenseNodePartition {
          seq: DenseSeqInfo::default(),
          profile: DenseSeqDistribution::default(),
        });
    }

    if let Some(info) = &changes.edge_merge {
      self.data.edges.remove(&info.parent_edge_key);
      self.data.edges.remove(&info.child_edge_key);
      self.data.nodes.remove(&info.removed_node_key);
      self.data.edges.entry(info.merged_edge_key).or_default();
    }

    Ok(())
  }
}

impl<N, E> TransitionCounting<N, E> for PartitionMarginalDense
where
  N: GraphNode,
  E: EdgeOptimizeOps,
{
  fn count_transitions(&self, graph: &Graph<N, E, ()>) -> Result<MutationCounts, Report> {
    self.data.count_transitions(graph)
  }
}

impl PartitionBranchOps for PartitionMarginalDense {
  fn sequence_length(&self) -> usize {
    self.length
  }

  fn edge_subs(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    let parent_non_char = &self.data.nodes[&parent_key].seq.non_char;
    let child_non_char = &self.data.nodes[&child_key].seq.non_char;

    let parent_profile = &self.data.nodes[&parent_key].profile.dis;
    let child_profile = &self.data.nodes[&child_key].profile.dis;
    let mut subs = Vec::new();

    for (pos, parent, child) in izip!(0..parent_profile.nrows(), parent_profile.rows(), child_profile.rows()) {
      let parent_state = self.alphabet.char(argmax_first(&parent).unwrap_or(0));
      let child_state = self.alphabet.char(argmax_first(&child).unwrap_or(0));
      if parent_state == child_state {
        continue;
      }
      if !self.alphabet.is_canonical(parent_state) || !self.alphabet.is_canonical(child_state) {
        continue;
      }
      if range_contains(parent_non_char, pos) || range_contains(child_non_char, pos) {
        continue;
      }

      subs.push(Sub::new(parent_state, pos, child_state)?);
    }

    Ok(subs)
  }

  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<crate::seq::indel::InDel> {
    self.data.edges[&edge_key].indels.clone()
  }

  fn root_sequence(&self, graph: &dyn BranchTopology) -> Result<Seq, Report> {
    Ok(assign_sequence(&self.data.nodes[&graph.root_key()?], &self.alphabet))
  }

  fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    assign_sequence(&self.data.nodes[&node_key], &self.alphabet)
  }

  fn edge_effective_length(&self, graph: &dyn BranchTopology, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    let (parent_key, child_key) = graph.edge_endpoints(edge_key)?;
    let parent_non_char = &self.data.nodes[&parent_key].seq.non_char;
    let child_non_char = &self.data.nodes[&child_key].seq.non_char;

    let non_char_positions: usize = range_union(&[parent_non_char.clone(), child_non_char.clone()])
      .iter()
      .map(|(start, end)| end - start)
      .sum();

    Ok(self.length.saturating_sub(non_char_positions))
  }
}

impl PartitionOptimizeOps for PartitionMarginalDense {
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    Ok(OptimizationContribution::from_dense(edge_key, self))
  }

  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.data.edges[&edge_key].indels.len()
  }
}

impl<N, E> PartitionTimetreeOps<N, E> for PartitionMarginalDense
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn reconcile_topology(&mut self, graph: &Graph<N, E, ()>) {
    let graph_node_keys: BTreeSet<GraphNodeKey> = graph.get_nodes().into_iter().map(|n| n.read_arc().key()).collect();
    let graph_edge_keys: BTreeSet<GraphEdgeKey> = graph.get_edges().into_iter().map(|e| e.read_arc().key()).collect();

    for &key in &graph_node_keys {
      self.data.nodes.entry(key).or_insert_with(|| DenseNodePartition {
        seq: DenseSeqInfo::default(),
        profile: DenseSeqDistribution::default(),
      });
    }

    for &key in &graph_edge_keys {
      self.data.edges.entry(key).or_default();
    }

    self.data.nodes.retain(|k, _| graph_node_keys.contains(k));
    self.data.edges.retain(|k, _| graph_edge_keys.contains(k));
  }
}

impl<N, E> PartitionMarginalOps<N, E> for PartitionMarginalDense
where
  N: NodeAncestralOps,
  E: EdgeOptimizeOps,
{
  fn attach_sequences(&mut self, graph: &Graph<N, E, ()>, aln: &[FastaRecord]) -> Result<(), Report> {
    let aln_by_name = aln.iter().fold(BTreeMap::new(), |mut records, record| {
      records.entry(record.seq_name.as_str()).or_insert(record);
      records
    });
    for leaf in graph.get_leaves() {
      let leaf_key = leaf.read_arc().key();
      let mut leaf = leaf.read_arc().payload().write_arc();

      let leaf_name = {
        let name = leaf.name().ok_or_else(|| {
          make_report!("Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name.")
        })?;
        name.as_ref().to_owned()
      };

      let leaf_fasta = aln_by_name
        .get(leaf_name.as_str())
        .copied()
        .ok_or_else(|| make_report!("Leaf sequence not found: '{leaf_name}'"))?;

      leaf.set_desc(leaf_fasta.desc.clone());

      let alphabet = &self.alphabet;
      self
        .data
        .nodes
        .insert(leaf_key, DenseNodePartition::new(&leaf_fasta.seq, alphabet)?);
    }

    for edge in graph.get_edges() {
      let edge_key = edge.read_arc().key();
      self.data.edges.insert(edge_key, DenseEdgePartition::default());
    }

    Ok(())
  }

  fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    if let Some(seq_info) = self.data.nodes.get(&node_key) {
      // Convergence checks and other intermediate reads always use the deterministic most-likely
      // state, independent of the user's output sampling mode.
      assign_sequence(seq_info, &self.alphabet)
    } else {
      seq! {}
    }
  }

  fn reconstruct_node_sequence(
    &mut self,
    node: &GraphNodeForward<N, E>,
    include_leaves: bool,
    impute: bool,
    sample_mode: SampleMode,
    rng: &mut dyn rand::RngCore,
  ) -> Option<Seq> {
    let seq = {
      let seq_info = self.data.nodes.get(&node.key)?;
      if node.is_leaf {
        // The dense tip keeps its observed input (gaps and unknowns already stamped), so it never
        // suffers the sparse tip corruption. Imputation resolves ambiguous/unknown positions (N and
        // IUPAC codes, not gaps) to the argmax of the leaf marginal posterior. The forward pass folds
        // the observed ambiguity mask into the leaf profile, so its argmax is the parent-informed
        // most likely state, matching v0.
        let mut seq = seq_info.seq.sequence.clone();
        if impute && seq_info.profile.dis.nrows() == seq.len() {
          for pos in 0..seq.len() {
            let ch = seq[pos];
            if !self.alphabet.is_canonical(ch) && !self.alphabet.is_gap(ch) {
              if let Some(idx) = argmax_first(&seq_info.profile.dis.row(pos)) {
                seq[pos] = self.alphabet.char(idx);
              }
            }
          }
        }
        seq
      } else {
        let sample = sample_mode.samples_node(node.is_root);
        assign_sequence_sampled(seq_info, &self.alphabet, sample, rng)
      }
    };

    // Persist the reconstruction so the node-data serializer reads back this flag-aware sequence
    // (the dense augur path returns `seq.sequence`), keeping the JSON and the reconstructed FASTA
    // consistent and matching the sparse backend.
    if let Some(node_data) = self.data.nodes.get_mut(&node.key) {
      node_data.seq.sequence = seq.clone();
    }

    // A suppressed tip is still reconstructed above (so the node-data serializer reads the corrected
    // sequence), but is not emitted to the reconstructed-FASTA visitor.
    if !include_leaves && node.is_leaf {
      return None;
    }

    Some(seq)
  }
}

/// Deterministic most-likely-state sequence assignment. Used by the forward pass and convergence
/// reads, which must stay reproducible regardless of the user's output sampling mode.
pub(crate) fn assign_sequence(seq_info: &DenseNodePartition, alphabet: &Alphabet) -> Seq {
  assign_sequence_sampled(seq_info, alphabet, false, &mut rand::thread_rng())
}

fn assign_sequence_sampled(
  seq_info: &DenseNodePartition,
  alphabet: &Alphabet,
  sample: bool,
  rng: &mut dyn rand::RngCore,
) -> Seq {
  let mut seq = prof2seq_sampled(&seq_info.profile, alphabet, sample, rng);
  for gap in &seq_info.seq.gaps {
    seq[gap.0..gap.1].fill(alphabet.gap());
  }
  for unk in &seq_info.seq.unknown {
    seq[unk.0..unk.1].fill(alphabet.unknown());
  }
  seq
}

fn prof2seq_sampled(
  profile: &DenseSeqDistribution,
  alphabet: &Alphabet,
  sample: bool,
  rng: &mut dyn rand::RngCore,
) -> Seq {
  let mut seq = seq! {};
  for row in profile.dis.rows() {
    seq.push(alphabet.char(resolve_profile(row, sample, rng)));
  }
  seq
}
