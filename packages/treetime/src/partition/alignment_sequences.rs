use crate::alphabet::alphabet::Alphabet;
use crate::error::OperationError;
use crate::partition::assembly::{AlignmentGroup, AlignmentLayout};
use crate::partition::marginal::posterior_states::AmbiguityThreshold;
use crate::partition::marginal::reconstruction::ReconstructionRef;
use crate::partition::result::AlignmentSequences;
use crate::seq::indel::InDel;
use crate::seq::mutation::{MutationEvent, MutationTrack, SequenceMutations, Sub, emit_sequences};
use crate::seq::sink::{AmbiguousSites, SeqStream, SeqTrack};
use eyre::Report;
use itertools::Itertools;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;
use treetime_utils::{make_internal_error, make_internal_report};

pub(crate) fn stream_alignment_sequences<P: PartitionSequences>(
  graph: &Graph,
  group: &AlignmentGroup,
  partitions: &[P],
  options: StreamOptions,
  sink: Option<SeqStream<'_>>,
) -> Result<AlignmentSequences, OperationError> {
  let track = SeqTrack::Nuc(group.alignment());
  let StreamOptions {
    include_leaves,
    report_unknown,
  } = options;
  let alphabet = partitions
    .first()
    .ok_or_else(no_partitions)
    .map_err(OperationError::InferenceFailed)?
    .alphabet();
  if group.layout.is_whole() {
    let mutations = partitions[0].stream_sequences(graph, track, include_leaves, report_unknown, sink)?;
    return AlignmentSequences::new(mutations, alphabet).map_err(OperationError::InferenceFailed);
  }
  let mutations = partitions
    .iter()
    .map(|partition| partition.stream_sequences(graph, track, include_leaves, report_unknown, None))
    .try_collect()?;
  let mutations = assemble_sequence_mutations(
    &group.layout,
    mutations,
    |member, edge_key| partitions[member].edge_indels(edge_key),
    &track.mutation_track(),
  )
  .map_err(OperationError::classify)?;
  if let Some(sink) = sink {
    emit_sequences(
      graph,
      track,
      include_leaves,
      |node_key| {
        let sequences: Vec<Seq> = partitions
          .iter()
          .map(|partition| partition.node_sequence(graph, node_key))
          .try_collect()?;
        group.layout.assemble_sequence(&sequences)
      },
      |node_key, threshold| {
        let sites: Vec<AmbiguousSites> = partitions
          .iter()
          .map(|partition| partition.node_ambiguity(graph, node_key, threshold))
          .try_collect()?;
        group.layout.assemble_ambiguity(&sites)
      },
      sink,
    )?;
  }
  AlignmentSequences::new(mutations, alphabet).map_err(OperationError::InferenceFailed)
}

pub(crate) fn assemble_sequence_mutations(
  layout: &AlignmentLayout,
  mutations: Vec<SequenceMutations>,
  edge_indels: impl Fn(usize, GraphEdgeKey) -> Vec<InDel>,
  track: &MutationTrack,
) -> Result<SequenceMutations, Report> {
  if layout.is_whole() {
    let [mutations] = <[SequenceMutations; 1]>::try_from(mutations)
      .map_err(|mutations| make_internal_report!("A whole alignment has {} partition results", mutations.len()))?;
    return Ok(mutations);
  }
  let root_sequences = mutations
    .iter()
    .map(|mutations| mutations.root_sequence.clone())
    .collect_vec();
  let root_sequence = layout.assemble_sequence(&root_sequences)?;
  let edge_mutations = mutations
    .first()
    .ok_or_else(no_partitions)?
    .edge_mutations
    .keys()
    .map(|&edge_key| -> Result<_, Report> {
      let substitutions = mutations
        .iter()
        .map(|mutations| edge_substitutions(mutations, edge_key))
        .try_collect::<_, Vec<_>, Report>()?;
      let indels = (0..mutations.len())
        .map(|member| edge_indels(member, edge_key))
        .collect_vec();
      Ok((
        edge_key,
        layout.assemble_edge_mutations(&substitutions, &indels, track)?,
      ))
    })
    .try_collect::<_, _, Report>()?;
  Ok(SequenceMutations {
    root_sequence,
    edge_mutations,
  })
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct StreamOptions {
  pub(crate) include_leaves: bool,
  pub(crate) report_unknown: bool,
}

impl<P: PartitionSequences + ?Sized> PartitionSequences for &P {
  fn alphabet(&self) -> &Alphabet {
    (**self).alphabet()
  }

  fn node_sequence(&self, graph: &Graph, node_key: GraphNodeKey) -> Result<Seq, Report> {
    (**self).node_sequence(graph, node_key)
  }

  fn node_ambiguity(
    &self,
    graph: &Graph,
    node_key: GraphNodeKey,
    threshold: AmbiguityThreshold,
  ) -> Result<AmbiguousSites, Report> {
    (**self).node_ambiguity(graph, node_key, threshold)
  }

  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    (**self).edge_indels(edge_key)
  }

  fn stream_sequences(
    &self,
    graph: &Graph,
    track: SeqTrack<'_>,
    include_leaves: bool,
    report_unknown: bool,
    sink: Option<SeqStream<'_>>,
  ) -> Result<SequenceMutations, OperationError> {
    (**self).stream_sequences(graph, track, include_leaves, report_unknown, sink)
  }
}

pub(crate) trait PartitionSequences {
  fn alphabet(&self) -> &Alphabet;

  fn node_sequence(&self, graph: &Graph, node_key: GraphNodeKey) -> Result<Seq, Report>;

  fn node_ambiguity(
    &self,
    graph: &Graph,
    node_key: GraphNodeKey,
    threshold: AmbiguityThreshold,
  ) -> Result<AmbiguousSites, Report>;

  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel>;

  fn stream_sequences(
    &self,
    graph: &Graph,
    track: SeqTrack<'_>,
    include_leaves: bool,
    report_unknown: bool,
    sink: Option<SeqStream<'_>>,
  ) -> Result<SequenceMutations, OperationError>;
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct ImputedReconstruction<'a> {
  pub(crate) reconstruction: ReconstructionRef<'a>,
  pub(crate) impute: bool,
}

impl PartitionSequences for ImputedReconstruction<'_> {
  fn alphabet(&self) -> &Alphabet {
    self.reconstruction.alphabet()
  }

  fn node_sequence(&self, graph: &Graph, node_key: GraphNodeKey) -> Result<Seq, Report> {
    self.reconstruction.node_sequence(graph, self.impute, node_key)
  }

  fn node_ambiguity(
    &self,
    graph: &Graph,
    node_key: GraphNodeKey,
    threshold: AmbiguityThreshold,
  ) -> Result<AmbiguousSites, Report> {
    self
      .reconstruction
      .node_ambiguity(graph, self.impute, node_key, threshold)
  }

  fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    self.reconstruction.edge_indels(edge_key)
  }

  fn stream_sequences(
    &self,
    graph: &Graph,
    track: SeqTrack<'_>,
    include_leaves: bool,
    report_unknown: bool,
    sink: Option<SeqStream<'_>>,
  ) -> Result<SequenceMutations, OperationError> {
    self
      .reconstruction
      .stream_sequences(graph, self.impute, track, include_leaves, report_unknown, sink)
  }
}

fn no_partitions() -> Report {
  make_internal_report!("An alignment has no partition results to assemble")
}

fn edge_substitutions(mutations: &SequenceMutations, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
  let Some(edge_mutations) = mutations.edge_mutations.get(&edge_key) else {
    return make_internal_error!("A partition of the alignment has no mutations on edge {edge_key}");
  };
  Ok(
    edge_mutations
      .iter()
      .filter_map(|mutation| match &mutation.event {
        MutationEvent::Substitution(substitution) => Some(substitution.clone()),
        MutationEvent::Insertion(_) | MutationEvent::Deletion(_) => None,
      })
      .collect(),
  )
}
