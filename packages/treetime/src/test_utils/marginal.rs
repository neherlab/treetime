use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::branch_lengths::branch_lengths_or_zero;
use crate::gtr::gtr::GTR;
use crate::partition::fitch::passes::create_fitch_partition;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction, SparseReconstruction};
use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate};
use crate::seq::alignment::{NodeSeqInput, node_seq_inputs};
use crate::seq::sink::{SeqItem, SeqSink};
use eyre::Report;
use std::collections::BTreeMap;
use std::sync::LazyLock;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::fasta_read;
use treetime_io::nwk::nwk_read;
use treetime_primitives::{AlignmentRecord, Seq, seq};

pub(crate) static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

pub(crate) fn run_dense_marginal_with_newick(newick: &str, aln_str: &str, gtr: &GTR) -> Result<f64, Report> {
  let nwk_parsed = nwk_read(newick.as_bytes())?;
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let aln: Vec<AlignmentRecord> = fasta_read(aln_str.as_bytes(), &*NUC_ALPHABET)?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;
  let partition = PartitionMarginalDense::new(0, alphabet, &graph, &node_seq_inputs(&graph, &names, aln))?;
  let MarginalUpdate { log_lh, .. } =
    partition.marginal_update(gtr, &graph, &branch_lengths_or_zero(&branch_lengths), &())?;
  Ok(log_lh.value())
}

pub(crate) fn run_sparse_marginal_with_newick(newick: &str, aln_str: &str, gtr: &GTR) -> Result<f64, Report> {
  let nwk_parsed = nwk_read(newick.as_bytes())?;
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let aln: Vec<AlignmentRecord> = fasta_read(aln_str.as_bytes(), &*NUC_ALPHABET)?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();
  let alphabet = Alphabet::new(AlphabetName::Nuc)?;

  let fitch = create_fitch_partition(&graph, 0, alphabet, node_seq_inputs(&graph, &names, aln))?;
  let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;

  let MarginalUpdate { log_lh, .. } =
    partition.marginal_update(gtr, &graph, &branch_lengths_or_zero(&branch_lengths), &node_states)?;
  Ok(log_lh.value())
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
  names: &BTreeMap<GraphNodeKey, Option<String>>,
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
      Ok((names[key].clone().expect("all test nodes are named"), seq))
    })
    .collect()
}

pub(crate) fn sparse_reconstruction(reconstruction: &MarginalReconstruction) -> &SparseReconstruction {
  match reconstruction {
    MarginalReconstruction::Sparse(sparse) => sparse,
    MarginalReconstruction::Dense(_) => panic!("expected a sparse reconstruction, got a dense one"),
  }
}

pub(crate) fn sparse_reconstruction_mut(reconstruction: &mut MarginalReconstruction) -> &mut SparseReconstruction {
  match reconstruction {
    MarginalReconstruction::Sparse(sparse) => sparse,
    MarginalReconstruction::Dense(_) => panic!("expected a sparse reconstruction, got a dense one"),
  }
}

pub(crate) fn dense_reconstruction_mut(reconstruction: &mut MarginalReconstruction) -> &mut DenseReconstruction {
  match reconstruction {
    MarginalReconstruction::Dense(dense) => dense,
    MarginalReconstruction::Sparse(_) => panic!("expected a dense reconstruction, got a sparse one"),
  }
}

pub(crate) fn dense_reconstruction(reconstruction: &MarginalReconstruction) -> &DenseReconstruction {
  match reconstruction {
    MarginalReconstruction::Dense(dense) => dense,
    MarginalReconstruction::Sparse(_) => panic!("expected a dense reconstruction, got a sparse one"),
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
        name: None,
        seq: Some(seq![fill; length]),
      };
      (leaf.key(), input)
    })
    .collect();
  PartitionMarginalDense::new(0, alphabet, graph, &node_inputs)
}
