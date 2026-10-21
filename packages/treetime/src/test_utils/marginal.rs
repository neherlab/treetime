use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::attach::{MissingLeafLimit, complete_alignment_for_leaves};
use crate::ancestral::mask::create_mask;
use crate::branch_lengths::branch_lengths_or_zero;
use crate::gtr::get_gtr::{JC69Params, jc69};
use crate::gtr::gtr::{GTR, Mu};
use crate::partition::fitch::passes::create_fitch_partition;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::reconstruction::{
  DenseReconstruction, MarginalReconstruction, ReconstructionRef, SparseReconstruction,
};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::set::PartitionSet;
use crate::partition::storage::sparse::SparseNodeState;
use crate::progress::NoopProgress;
use crate::seq::alignment::{
  AncestralInput, EdgeSeqInput, NodeSeqInput, SequenceSource, get_common_length_of_node_inputs, pair_leaf_sequences,
};
use crate::seq::sink::{SeqItem, SeqSink};
use eyre::Report;
use itertools::Itertools;
use ndarray::array;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::LazyLock;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_io::fasta::fasta_read;
use treetime_io::nwk::{NwkParse, nwk_read};
use treetime_primitives::{AlignmentRecord, LogLh, Seq, seq};

pub(crate) static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

pub(crate) fn run_dense_marginal(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &NodeNames,
  aln: &[AlignmentRecord],
  gtr: GTR,
) -> Result<(f64, DenseReconstruction), Report> {
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;
  let partition = PartitionMarginalDense::new(alphabet, graph, &leaf_seq_inputs(graph, names, aln.to_vec()))?;
  let (reconstruction, log_lh) =
    DenseReconstruction::seeded(partition, gtr).marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
  Ok((log_lh.value(), reconstruction))
}

pub(crate) fn run_sparse_marginal(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &NodeNames,
  aln: &[AlignmentRecord],
  gtr: GTR,
) -> Result<(f64, SparseReconstruction), Report> {
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;
  let fitch = create_fitch_partition(graph, alphabet, leaf_seq_inputs(graph, names, aln.to_vec()))?;
  let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
  let (reconstruction, log_lh) = SparseReconstruction::seeded(partition, gtr, node_states)
    .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
  Ok((log_lh.value(), reconstruction))
}

pub(crate) fn run_dense_marginal_with_newick(
  newick: &str,
  aln_str: &str,
  gtr: &GTR,
) -> Result<(f64, DenseReconstruction), Report> {
  let NwkParse {
    graph,
    branch_lengths,
    names,
    ..
  } = nwk_read(newick.as_bytes())?;
  run_dense_marginal(&graph, &branch_lengths, &names, &nuc_records(aln_str)?, gtr.clone())
}

pub(crate) fn run_sparse_marginal_with_newick(
  newick: &str,
  aln_str: &str,
  gtr: &GTR,
) -> Result<(f64, SparseReconstruction), Report> {
  let NwkParse {
    graph,
    branch_lengths,
    names,
    ..
  } = nwk_read(newick.as_bytes())?;
  run_sparse_marginal(&graph, &branch_lengths, &names, &nuc_records(aln_str)?, gtr.clone())
}

pub(crate) fn make_nonuniform_gtr() -> Result<GTR, Report> {
  GTR::builder()
    .n_states(NUC_ALPHABET.n_canonical())
    .mu(Mu::FromMatrix(1.0))
    .pi(array![0.2, 0.3, 0.15, 0.35])
    .build()
}

fn nuc_records(aln_str: &str) -> Result<Vec<AlignmentRecord>, Report> {
  Ok(
    fasta_read(aln_str.as_bytes(), &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect(),
  )
}

pub(crate) fn marginal_update_reconstruction(
  reconstruction: MarginalReconstruction,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<(MarginalReconstruction, LogLh), Report> {
  Ok(match reconstruction {
    MarginalReconstruction::Dense(reconstruction) => {
      let (reconstruction, log_lh) = reconstruction.marginal_update(graph, branch_lengths)?;
      (MarginalReconstruction::Dense(reconstruction), log_lh)
    },
    MarginalReconstruction::Sparse(reconstruction) => {
      let (reconstruction, log_lh) = reconstruction.marginal_update(graph, branch_lengths)?;
      (MarginalReconstruction::Sparse(reconstruction), log_lh)
    },
  })
}

#[derive(Default)]
pub(crate) struct RecordingSeqSink {
  pub(crate) items: Vec<(GraphNodeKey, bool, Seq)>,
}

impl SeqSink for RecordingSeqSink {
  fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
    self.items.push((item.key, item.emitted, item.seq.clone()));
    Ok(())
  }
}

pub(crate) fn internal_node_keys(graph: &Graph) -> Vec<GraphNodeKey> {
  graph.get_internal_nodes().map(|node| node.key()).collect()
}

pub(crate) fn node_keys(graph: &Graph) -> Vec<GraphNodeKey> {
  graph.get_nodes().map(|node| node.key()).collect()
}

pub(crate) fn emitted_sequences_by_name(
  names: &NodeNames,
  emitted_nodes: &[GraphNodeKey],
  sampled: &BTreeMap<GraphNodeKey, Seq>,
  node_sequence: impl Fn(GraphNodeKey) -> Result<Seq, Report>,
) -> Result<BTreeMap<String, Seq>, Report> {
  emitted_nodes
    .iter()
    .map(|key| {
      let seq = match sampled.get(key) {
        Some(seq) => seq.clone(),
        None => node_sequence(*key)?,
      };
      Ok((names[*key].to_owned(), seq))
    })
    .collect()
}

pub(crate) fn sparse_reconstruction(reconstruction: &impl SingleReconstruction) -> &SparseReconstruction {
  match reconstruction.single() {
    ReconstructionRef::Sparse(sparse) => sparse,
    ReconstructionRef::Dense(_) => panic!("expected a sparse reconstruction, got a dense one"),
  }
}

pub(crate) fn sparse_reconstruction_mut(reconstruction: &mut impl SingleReconstruction) -> &mut SparseReconstruction {
  match reconstruction.single_mut() {
    ReconstructionMut::Sparse(sparse) => sparse,
    ReconstructionMut::Dense(_) => panic!("expected a sparse reconstruction, got a dense one"),
  }
}

pub(crate) fn dense_reconstruction_mut(reconstruction: &mut impl SingleReconstruction) -> &mut DenseReconstruction {
  match reconstruction.single_mut() {
    ReconstructionMut::Dense(dense) => dense,
    ReconstructionMut::Sparse(_) => panic!("expected a dense reconstruction, got a sparse one"),
  }
}

pub(crate) fn dense_reconstruction(reconstruction: &impl SingleReconstruction) -> &DenseReconstruction {
  match reconstruction.single() {
    ReconstructionRef::Dense(dense) => dense,
    ReconstructionRef::Sparse(_) => panic!("expected a dense reconstruction, got a sparse one"),
  }
}

pub(crate) trait SingleReconstruction {
  fn single(&self) -> ReconstructionRef<'_>;

  fn single_mut(&mut self) -> ReconstructionMut<'_>;

  fn gtr_mut(&mut self) -> &mut GTR {
    match self.single_mut() {
      ReconstructionMut::Dense(dense) => &mut dense.gtr,
      ReconstructionMut::Sparse(sparse) => &mut sparse.gtr,
    }
  }
}

pub(crate) enum ReconstructionMut<'a> {
  Dense(&'a mut DenseReconstruction),
  Sparse(&'a mut SparseReconstruction),
}

impl SingleReconstruction for MarginalReconstruction {
  fn single(&self) -> ReconstructionRef<'_> {
    self.view()
  }

  fn single_mut(&mut self) -> ReconstructionMut<'_> {
    match self {
      Self::Dense(dense) => ReconstructionMut::Dense(dense),
      Self::Sparse(sparse) => ReconstructionMut::Sparse(sparse),
    }
  }
}

impl SingleReconstruction for PartitionSet {
  fn single(&self) -> ReconstructionRef<'_> {
    self
      .iter()
      .exactly_one()
      .ok()
      .expect("expected a partition set with exactly one partition")
  }

  fn single_mut(&mut self) -> ReconstructionMut<'_> {
    let single = match self {
      Self::Dense(dense) => dense.iter_mut().exactly_one().ok().map(ReconstructionMut::Dense),
      Self::Sparse(sparse) => sparse.iter_mut().exactly_one().ok().map(ReconstructionMut::Sparse),
    };
    single.expect("expected a partition set with exactly one partition")
  }
}

pub(crate) fn dense_partition_with_constant_leaves(
  graph: &Graph,
  alphabet: Alphabet,
  length: usize,
) -> Result<PartitionMarginalDense, Report> {
  let fill = alphabet.char(0);
  let node_inputs = graph
    .get_leaves()
    .map(|leaf| {
      let input = NodeSeqInput {
        name: leaf.key().to_string(),
        seq: Some(seq![fill; length]),
      };
      (leaf.key(), input)
    })
    .collect();
  PartitionMarginalDense::new(alphabet, graph, &node_inputs)
}

pub(crate) fn leaf_seq_inputs(
  graph: &Graph,
  names: &NodeNames,
  aln: Vec<AlignmentRecord>,
) -> BTreeMap<GraphNodeKey, NodeSeqInput> {
  pair_leaf_sequences(graph, names, aln).sequences.nodes
}

pub(crate) struct CompletedSequences {
  pub(crate) nodes: BTreeMap<GraphNodeKey, NodeSeqInput>,
  pub(crate) mask: Vec<bool>,
}

pub(crate) fn complete_sequence_source(
  graph: &Graph,
  names: &NodeNames,
  aln: Vec<AlignmentRecord>,
  alphabet: &Alphabet,
  path: &str,
) -> Result<SequenceSource, Report> {
  let CompletedSequences { nodes, mask } = complete_leaf_sequences(graph, names, aln, alphabet)?;
  Ok(SequenceSource {
    path: PathBuf::from(path),
    alphabet: alphabet.clone(),
    length: mask.len(),
    nodes,
    mask,
  })
}

pub(crate) fn sequence_source(
  alphabet: Alphabet,
  nodes: BTreeMap<GraphNodeKey, NodeSeqInput>,
) -> Result<SequenceSource, Report> {
  let length = get_common_length_of_node_inputs(&nodes)?;
  let mask = create_mask(nodes.values().filter_map(|node| node.seq.as_ref()), length, &alphabet);
  Ok(SequenceSource {
    path: PathBuf::from("aln.fasta"),
    alphabet,
    length,
    nodes,
    mask,
  })
}

pub(crate) fn ancestral_input_for_tree(parse: NwkParse, sources: Vec<SequenceSource>) -> AncestralInput {
  AncestralInput {
    edges: parse
      .branch_lengths
      .into_iter()
      .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
      .collect(),
    graph: parse.graph,
    sources,
  }
}

pub(crate) fn complete_leaf_sequences(
  graph: &Graph,
  names: &NodeNames,
  aln: Vec<AlignmentRecord>,
  alphabet: &Alphabet,
) -> Result<CompletedSequences, Report> {
  let mut sequences = pair_leaf_sequences(graph, names, aln).sequences;
  let length = sequences.common_length()?;
  complete_alignment_for_leaves(
    graph,
    &mut sequences.nodes,
    length,
    alphabet,
    MissingLeafLimit::Enforced,
    "the alignment",
    &NoopProgress,
  )?;
  let mask = sequences.mask(length, alphabet);
  Ok(CompletedSequences {
    nodes: sequences.nodes,
    mask,
  })
}

impl From<MarginalReconstruction> for PartitionSet {
  fn from(reconstruction: MarginalReconstruction) -> Self {
    match reconstruction {
      MarginalReconstruction::Dense(dense) => Self::Dense(vec![dense]),
      MarginalReconstruction::Sparse(sparse) => Self::Sparse(vec![sparse]),
    }
  }
}

pub(crate) fn sparse_topology_partitions(partitions: Vec<PartitionMarginalSparse>) -> Vec<SparseReconstruction> {
  partitions
    .into_iter()
    .map(|partition| sparse_topology_partition(partition, BTreeMap::new()))
    .collect()
}

pub(crate) fn sparse_topology_partition(
  partition: PartitionMarginalSparse,
  node_states: BTreeMap<GraphNodeKey, SparseNodeState>,
) -> SparseReconstruction {
  SparseReconstruction::seeded(
    partition,
    jc69(JC69Params::default()).expect("JC69 builds"),
    node_states,
  )
}
