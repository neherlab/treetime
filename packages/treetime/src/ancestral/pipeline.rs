use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::{ancestral_reconstruction_fitch, create_fitch_partition};
use crate::ancestral::marginal::{ancestral_reconstruction, branch_lengths_or_zero};
use crate::ancestral::params::MethodAncestral;
use crate::ancestral::sample::SampleMode;
use crate::cancel::Cancel;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::gtr::refinement::refine_gtr_model;
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::fitch::partition::PartitionFitch;
use crate::partition::marginal::dense::partition::{DenseMarginalEdges, PartitionMarginalDense};
use crate::partition::marginal::shared::update::{MarginalPasses, MarginalStates, MarginalUpdate};
use crate::partition::marginal::sparse::partition::{PartitionMarginalSparse, SparseMarginalEdges};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::dense::DenseNodeState;
use crate::partition::storage::sparse::SparseNodeState;
use crate::progress::ProgressSink;
use crate::seq::alignment::ReconstructionInput;
use crate::seq::indel::InDel;
use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use strum::VariantNames;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::AsciiChar;
use treetime_primitives::LogLh;
use treetime_primitives::Seq;
use treetime_utils::make_error;
use treetime_utils::sync::random::get_random_number_generator;

pub struct AncestralParams {
  pub method: MethodAncestral,
  pub model: GtrModelName,
  pub dense: Option<bool>,
  /// Emit reconstructed leaf sequences in addition to internal nodes.
  pub include_leaves: bool,
  /// Resolve ambiguous and unknown tip states (`N` and IUPAC codes) to the most likely inferred
  /// state. Only defined for marginal reconstruction; a no-op for Fitch parsimony.
  pub impute_missing_data: bool,
  pub gtr_iterations: usize,
  pub site_specific_gtr: bool,
  pub seed: Option<u64>,
  pub sample_from_profile: SampleMode,
  pub ignore_missing_alns: bool,
}

/// A sparse reconstruction: the durable partition inputs, the node states carried between passes, and
/// the per-edge results of the last pass, as distinct owned values. The output writers build a
/// short-lived read view over these.
#[derive(Clone, Debug, Serialize)]
pub struct SparseReconstruction {
  pub partition: PartitionMarginalSparse,
  /// The substitution model this reconstruction was produced under, carried as a value alongside the
  /// immutable partition.
  pub gtr: GTR,
  pub node_states: BTreeMap<GraphNodeKey, SparseNodeState>,
  pub edges: SparseMarginalEdges,
}

impl SparseReconstruction {
  /// A reconstruction seeded from the Fitch handoff, before any marginal pass has run: durable
  /// observations and leaf node states, with no per-edge results yet.
  pub fn seeded(
    partition: PartitionMarginalSparse,
    gtr: GTR,
    node_states: BTreeMap<GraphNodeKey, SparseNodeState>,
  ) -> Self {
    Self {
      partition,
      gtr,
      node_states,
      edges: SparseMarginalEdges::default(),
    }
  }

  /// The sequence length this reconstruction represents.
  pub fn sequence_length(&self) -> usize {
    self.partition.length
  }

  /// MAP-derived nucleotide substitutions for one edge, read from the forward-pass estimates.
  pub fn edge_subs(&self, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    self.partition.edge_subs(&self.edges.estimates, edge_key)
  }

  /// The number of alignment positions where both endpoints carry canonical states for one edge.
  pub fn edge_effective_length(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    self.partition.edge_effective_length(graph, edge_key)
  }

  /// The per-edge branch-length optimization contribution, built from the last update's messages.
  pub fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    self
      .partition
      .create_edge_contribution(&self.gtr, &self.edges.backward, &self.edges.forward, edge_key)
  }

  /// The number of indel events on one edge.
  pub fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.partition.edge_indel_count(edge_key)
  }

  /// The reconstructed sequence for one node, resolved against the posterior (MAP or sampled draw).
  pub fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    self.partition.node_sequence(&self.node_states, node_key)
  }

  /// The reconstructed root sequence.
  pub fn root_sequence(&self, _graph: &Graph) -> Result<Seq, Report> {
    Ok(self.partition.root_sequence())
  }

  /// Grouped aligned insertions and deletions for one edge.
  pub fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    self.partition.edge_indels(edge_key)
  }

  /// The substitutions and indels on one edge as one mutation list on the given track.
  pub fn edge_mutations(&self, edge_key: GraphEdgeKey, track: &MutationTrack) -> Result<Vec<Mutation>, Report> {
    combine_edge_mutations(self.edge_subs(edge_key)?, &self.edge_indels(edge_key), track)
  }

  /// The node sequence written into the augur node-data JSON. For the sparse representation this equals
  /// [`Self::node_sequence`]: both resolve the parsimony chain against the posterior, so the JSON and
  /// the reconstructed FASTA carry the same MAP states.
  pub fn augur_node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    self.node_sequence(node_key)
  }

  /// The alphabet's ambiguous (unknown) character, used to fill masked positions in output sequences.
  pub fn ambiguous_char(&self) -> AsciiChar {
    self.partition.alphabet.unknown()
  }

  /// Run a full marginal update, returning the reconstruction at the refreshed node states and per-edge
  /// results together with the substitution log likelihood.
  ///
  /// The reconstruction is consumed and a new one returned, so a failed pass produces no reconstruction
  /// at all rather than one whose maps come from different passes.
  pub fn marginal_update(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<(Self, LogLh), Report> {
    let Self {
      partition,
      gtr,
      node_states,
      ..
    } = self;
    let MarginalUpdate {
      node_states,
      edges,
      log_lh,
    } = partition.marginal_update(&gtr, graph, branch_lengths, node_states)?;
    Ok((
      Self {
        partition,
        gtr,
        node_states,
        edges,
      },
      log_lh,
    ))
  }
}

/// A dense reconstruction: the durable partition inputs, the node states carried between passes, and
/// the per-edge results of the last pass, as distinct owned values.
#[derive(Clone, Debug, Serialize)]
pub struct DenseReconstruction {
  pub partition: PartitionMarginalDense,
  /// The substitution model this reconstruction was produced under, carried as a value alongside the
  /// immutable partition.
  pub gtr: GTR,
  pub node_states: BTreeMap<GraphNodeKey, DenseNodeState>,
  pub edges: DenseMarginalEdges,
}

impl DenseReconstruction {
  /// A reconstruction seeded from the alignment, before any marginal pass has run: durable inputs and
  /// leaf node states, with no per-edge results yet.
  pub fn seeded(
    partition: PartitionMarginalDense,
    gtr: GTR,
    node_states: BTreeMap<GraphNodeKey, DenseNodeState>,
  ) -> Self {
    Self {
      partition,
      gtr,
      node_states,
      edges: DenseMarginalEdges::default(),
    }
  }

  /// The sequence length this reconstruction represents.
  pub fn sequence_length(&self) -> usize {
    self.partition.length
  }

  /// MAP-derived nucleotide substitutions for one edge, read from the node-state profiles.
  pub fn edge_subs(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    self.partition.edge_subs(&self.node_states, graph, edge_key)
  }

  /// The number of alignment positions where both endpoints carry canonical states for one edge.
  pub fn edge_effective_length(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    self.partition.edge_effective_length(&self.node_states, graph, edge_key)
  }

  /// The per-edge branch-length optimization contribution, built from the last update's messages.
  pub fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> OptimizationContribution {
    self
      .partition
      .create_edge_contribution(&self.gtr, &self.edges.backward, &self.edges.forward, edge_key)
  }

  /// The number of indel events on one edge.
  pub fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.partition.edge_indel_count(&self.edges.estimates, edge_key)
  }

  /// The reconstructed most-likely-state sequence for one node (deterministic MAP).
  pub fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    self.partition.node_sequence(&self.node_states, node_key)
  }

  /// The reconstructed root sequence (deterministic MAP).
  pub fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    self.partition.root_sequence(&self.node_states, graph)
  }

  /// Grouped aligned insertions and deletions for one edge.
  pub fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    self.partition.edge_indels(&self.edges.estimates, edge_key)
  }

  /// The substitutions and indels on one edge as one mutation list on the given track.
  pub fn edge_mutations(
    &self,
    graph: &Graph,
    edge_key: GraphEdgeKey,
    track: &MutationTrack,
  ) -> Result<Vec<Mutation>, Report> {
    combine_edge_mutations(self.edge_subs(graph, edge_key)?, &self.edge_indels(edge_key), track)
  }

  /// The node sequence written into the augur node-data JSON. Unlike [`Self::node_sequence`] (which
  /// re-derives the MAP state from the profile), this reads back the flag-aware sequence the marginal
  /// reconstruction pass stored in `seq.sequence` (observed echo or imputation), keeping the JSON
  /// consistent with the reconstructed FASTA and with the sparse backend.
  pub fn augur_node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    self.node_states[&node_key].seq.sequence.clone()
  }

  /// The alphabet's ambiguous (unknown) character, used to fill masked positions in output sequences.
  pub fn ambiguous_char(&self) -> AsciiChar {
    self.partition.alphabet.unknown()
  }

  /// Run a full marginal update, returning the reconstruction at the refreshed node states and per-edge
  /// results together with the substitution log likelihood.
  ///
  /// The reconstruction is consumed and a new one returned, so a failed pass produces no reconstruction
  /// at all rather than one whose maps come from different passes.
  pub fn marginal_update(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<(Self, LogLh), Report> {
    let Self {
      partition,
      gtr,
      node_states,
      ..
    } = self;
    let MarginalUpdate {
      node_states,
      edges,
      log_lh,
    } = partition.marginal_update(&gtr, graph, branch_lengths, node_states)?;
    Ok((
      Self {
        partition,
        gtr,
        node_states,
        edges,
      },
      log_lh,
    ))
  }
}

#[derive(Clone, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum AncestralPartition {
  Fitch(PartitionFitch),
  Sparse(SparseReconstruction),
  Dense(DenseReconstruction),
}

/// Output-side read access over a completed reconstruction, dispatching each operation to the concrete
/// representation. The tree writers read `node_sequence`/`root_sequence`/`edge_mutations`; the augur
/// node-data writer reads `augur_node_sequence`/`augur_root_sequence`/`edge_subs` and the alphabet's
/// ambiguous character. The tree and augur node/root sequences differ for the dense representation (MAP
/// states versus the stored flag-aware reconstruction) and are kept as distinct accessors.
impl AncestralPartition {
  /// The alignment length (number of sites).
  pub fn sequence_length(&self) -> usize {
    match self {
      Self::Fitch(partition) => partition.sequence_length(),
      Self::Sparse(partition) => partition.sequence_length(),
      Self::Dense(partition) => partition.sequence_length(),
    }
  }

  /// The alphabet's ambiguous (unknown) character.
  pub fn ambiguous_char(&self) -> AsciiChar {
    match self {
      Self::Fitch(partition) => partition.ambiguous_char(),
      Self::Sparse(partition) => partition.ambiguous_char(),
      Self::Dense(partition) => partition.ambiguous_char(),
    }
  }

  /// The reconstructed sequence for one node, as read by the tree writers.
  pub fn node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    match self {
      Self::Fitch(partition) => partition.node_sequence(node_key),
      Self::Sparse(partition) => partition.node_sequence(node_key),
      Self::Dense(partition) => partition.node_sequence(node_key),
    }
  }

  /// The reconstructed sequence for one node, as written into the augur node-data JSON.
  pub fn augur_node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    match self {
      Self::Fitch(partition) => partition.node_sequence(node_key),
      Self::Sparse(partition) => partition.augur_node_sequence(node_key),
      Self::Dense(partition) => partition.augur_node_sequence(node_key),
    }
  }

  /// The reconstructed root sequence, as read by the tree writers.
  pub fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    match self {
      Self::Fitch(partition) => partition.root_sequence(graph),
      Self::Sparse(partition) => partition.root_sequence(graph),
      Self::Dense(partition) => partition.root_sequence(graph),
    }
  }

  /// The reconstructed root sequence used as the augur JSON reference, matching the root node's augur
  /// sequence.
  pub fn augur_root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    Ok(self.augur_node_sequence(graph.root_key()?))
  }

  /// MAP-derived nucleotide substitutions on one edge (parent -> child).
  pub fn edge_subs(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    match self {
      Self::Fitch(partition) => partition.edge_subs(graph, edge_key),
      Self::Sparse(partition) => partition.edge_subs(edge_key),
      Self::Dense(partition) => partition.edge_subs(graph, edge_key),
    }
  }

  /// Grouped aligned insertions and deletions for one edge.
  pub fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    match self {
      Self::Fitch(partition) => partition.edge_indels(edge_key),
      Self::Sparse(partition) => partition.edge_indels(edge_key),
      Self::Dense(partition) => partition.edge_indels(edge_key),
    }
  }

  /// The substitutions and indels on one edge as one mutation list on the given track.
  pub fn edge_mutations(
    &self,
    graph: &Graph,
    edge_key: GraphEdgeKey,
    track: &MutationTrack,
  ) -> Result<Vec<Mutation>, Report> {
    combine_edge_mutations(self.edge_subs(graph, edge_key)?, &self.edge_indels(edge_key), track)
  }
}

#[derive(Debug, Serialize)]
pub struct AncestralOutput {
  #[serde(skip)]
  pub gtr: Option<GTR>,
  pub model_name: GtrModelName,
  #[serde(skip)]
  pub mask: Vec<bool>,
  /// Node ids in the order the reconstruction walk emits them (depth-first preorder, suppressed tips
  /// excluded). The reconstructed-FASTA writer replays this order and reads each sequence back off the
  /// partition, so the streamed records match the walk without holding every sequence in memory.
  #[serde(skip)]
  pub emitted_nodes: Vec<GraphNodeKey>,
}

pub struct AncestralOutputFull {
  pub output: AncestralOutput,
  pub partition: Option<AncestralPartition>,
}

pub fn run(
  params: &AncestralParams,
  input: &ReconstructionInput,
  alphabet: Alphabet,
  mask: Vec<bool>,
  cancel: &dyn Cancel,
  progress: &dyn ProgressSink,
) -> Result<AncestralOutputFull, Report> {
  let branch_lengths = input.branch_lengths();
  let profile_lengths = branch_lengths_or_zero(&branch_lengths);
  if params.site_specific_gtr {
    return make_error!(
      "--site-specific-gtr is not yet integrated into the ancestral reconstruction pipeline. \
       The mathematical core (GTRSiteSpecific) is implemented but partition system wiring is pending."
    );
  }

  if params.sample_from_profile != SampleMode::Argmax && params.method != MethodAncestral::Marginal {
    return make_error!(
      "--sample-from-profile={:?} requires --method-anc=marginal. Posterior sampling is only defined \
       for marginal reconstruction; {:?} has no posterior profile to sample.",
      params.sample_from_profile,
      params.method
    );
  }

  // The caller completes the alignment (fills fully-ambiguous sequences for tips absent from it) and
  // computes the mask before building the merged input, so every leaf's node input carries a
  // sequence and attachment always finds one by node key.
  let graph = &input.graph;
  let node_inputs = &input.nodes;
  let mut rng = get_random_number_generator(params.seed);

  match params.method {
    MethodAncestral::Parsimony => {
      cancel.check()?;
      progress.report("Fitch parsimony", 0.3, "");
      let partition = create_fitch_partition(graph, 0, alphabet, node_inputs)?;
      let mut partitions_parsimony = vec![partition];

      if params.impute_missing_data {
        log::warn!(
          "--impute-missing-data has no effect with --method-anc=parsimony: Fitch parsimony produces no \
           posterior profile to impute missing tip states from. Leaf states are emitted as observed."
        );
      }

      let emitted_nodes = ancestral_reconstruction_fitch(graph, params.include_leaves, &mut partitions_parsimony)?;

      let partition = partitions_parsimony
        .into_iter()
        .next()
        .expect("partition vec not empty");
      progress.report("Done", 1.0, "");
      Ok(AncestralOutputFull {
        output: AncestralOutput {
          gtr: None,
          model_name: params.model,
          mask,
          emitted_nodes,
        },
        partition: Some(AncestralPartition::Fitch(partition)),
      })
    },
    MethodAncestral::Marginal => {
      cancel.check()?;
      progress.report("Inferring GTR model", 0.2, "");

      let created = create_marginal_partition(
        graph,
        0,
        alphabet,
        node_inputs,
        params.model,
        params.dense,
        &profile_lengths,
      )?;
      let model_name = created.model_name;
      let gtr = created.gtr;
      let refine = params.gtr_iterations > 0 && params.model == GtrModelName::Infer;

      match created.partition {
        MarginalPartition::Sparse(partition, node_states) => {
          cancel.check()?;
          progress.report("Marginal reconstruction", 0.4, "");
          let update = partition.marginal_update(&gtr, graph, &profile_lengths, node_states)?;

          let (gtr, update) = if refine {
            refine_gtr_model(
              &partition,
              gtr,
              update,
              params.gtr_iterations,
              1.0,
              graph,
              &profile_lengths,
            )?
          } else {
            (gtr, update)
          };
          let MarginalUpdate {
            mut node_states, edges, ..
          } = update;

          cancel.check()?;
          progress.report("Reconstructing sequences", 0.6, "");
          let emitted_nodes = ancestral_reconstruction(graph, |node| {
            partition.advance_node_state(
              &mut node_states,
              &edges.forward,
              node,
              params.include_leaves,
              params.impute_missing_data,
              params.sample_from_profile,
              &mut rng,
            )
          })?;

          progress.report("Done", 1.0, "");
          Ok(AncestralOutputFull {
            output: AncestralOutput {
              gtr: Some(gtr.clone()),
              model_name,
              mask,
              emitted_nodes,
            },
            partition: Some(AncestralPartition::Sparse(SparseReconstruction {
              partition,
              gtr,
              node_states,
              edges,
            })),
          })
        },
        MarginalPartition::Dense(partition) => {
          cancel.check()?;
          progress.report("Marginal reconstruction", 0.4, "");
          let node_states = partition.attach_sequences(graph, node_inputs)?;
          // Dense gap classification is non-idempotent, so the baseline ran two marginal passes
          // after attachment (its `initialize_marginal` attached and updated once, then a separate
          // `marginal_update` ran again) before GTR refinement. Pass the node states through both
          // passes so internal-node gap states settle exactly as they did before.
          let MarginalStates { node_states, .. } =
            partition.marginal_states(&gtr, graph, &profile_lengths, node_states)?;
          let update = partition.marginal_update(&gtr, graph, &profile_lengths, node_states)?;

          let (gtr, update) = if refine {
            refine_gtr_model(
              &partition,
              gtr,
              update,
              params.gtr_iterations,
              1.0,
              graph,
              &profile_lengths,
            )?
          } else {
            (gtr, update)
          };
          let MarginalUpdate {
            mut node_states, edges, ..
          } = update;

          cancel.check()?;
          progress.report("Reconstructing sequences", 0.6, "");
          // Dense reconstruction persists each node's flag-aware sequence into `seq.sequence` (read
          // back by the augur node-data path), so the walk still materializes it and only discards the
          // returned value; the emitted-node order is all the FASTA writer needs.
          let emitted_nodes = ancestral_reconstruction(graph, |node| {
            partition
              .reconstruct_node_sequence(
                &mut node_states,
                node,
                params.include_leaves,
                params.impute_missing_data,
                params.sample_from_profile,
                &mut rng,
              )
              .map(|_| ())
          })?;

          progress.report("Done", 1.0, "");
          Ok(AncestralOutputFull {
            output: AncestralOutput {
              gtr: Some(gtr.clone()),
              model_name,
              mask,
              emitted_nodes,
            },
            partition: Some(AncestralPartition::Dense(DenseReconstruction {
              partition,
              gtr,
              node_states,
              edges,
            })),
          })
        },
      }
    },
    MethodAncestral::Joint => {
      let available = MethodAncestral::VARIANTS
        .iter()
        .filter(|v| **v != "joint")
        .copied()
        .collect::<Vec<_>>()
        .join(", ");
      make_error!("Joint ancestral reconstruction has been removed. Available methods: {available}")
    },
  }
}
