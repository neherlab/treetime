use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::gtr::get_gtr::{GtrModel, log_gtr};
use crate::gtr::gtr::GTR;
use crate::partition::algo::infer_dense::infer_dense;
use crate::partition::fitch::gtr_inference::infer_gtr_fitch;
use crate::partition::fitch::partition::PartitionFitch;
use crate::partition::fitch::passes::{create_fitch_partition, create_inherited_fitch_partition};
use crate::partition::inferred_rate::InferredRate;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction, SparseReconstruction};
use crate::progress::LogSink;
use crate::seq::alignment::NodeSeqInput;
use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::scope::Inheritance;

#[expect(
  clippy::too_many_arguments,
  reason = "the partition takes its representation, model and rate rule, and its data (graph, alphabet, leaf inputs, branch lengths) and the log sink as separate values; no type holds them together"
)]
pub(crate) fn build_marginal_partition(
  representation: Representation,
  model: &GtrModel,
  rate: InferredRate,
  graph: &Graph,
  alphabet: Alphabet,
  node_inputs: NodeMap<NodeSeqInput>,
  branch_lengths: &EdgeMap<f64>,
  log: &dyn LogSink,
) -> Result<MarginalReconstruction, Report> {
  let data = PartitionSource {
    graph,
    inheritance: Inheritance::ALL,
    alphabet,
    node_inputs,
    branch_lengths,
  };
  build_inherited_marginal_partition(representation, model, rate, data, log)
}

pub(crate) struct PartitionSource<'a> {
  pub(crate) graph: &'a Graph,
  pub(crate) inheritance: Inheritance,
  pub(crate) alphabet: Alphabet,
  pub(crate) node_inputs: NodeMap<NodeSeqInput>,
  pub(crate) branch_lengths: &'a EdgeMap<f64>,
}

pub(crate) fn build_inherited_marginal_partition(
  representation: Representation,
  model: &GtrModel,
  rate: InferredRate,
  data: PartitionSource<'_>,
  log: &dyn LogSink,
) -> Result<MarginalReconstruction, Report> {
  let PartitionSource {
    graph,
    inheritance,
    alphabet,
    node_inputs,
    branch_lengths,
  } = data;
  Ok(match representation {
    Representation::Sparse => {
      let fitch = create_inherited_fitch_partition(graph, inheritance, alphabet, node_inputs)?;
      MarginalReconstruction::Sparse(sparse_from_fitch(model, rate, fitch, graph, branch_lengths, log)?)
    },
    Representation::Dense if model.is_infer() => {
      let fitch = create_inherited_fitch_partition(graph, inheritance, alphabet, node_inputs.clone())?;
      let gtr = fitch_gtr(model, rate, &fitch, graph, branch_lengths, log)?;
      let partition = fitch.into_marginal_dense(graph, &node_inputs)?;
      MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, gtr))
    },
    Representation::Dense => {
      let alphabet_name = alphabet.name();
      let partition = PartitionMarginalDense::new(alphabet, graph, &node_inputs)?.with_inheritance(inheritance);
      let gtr = named_gtr(model, alphabet_name, log)?;
      MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, gtr))
    },
  })
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum Representation {
  Dense,
  Sparse,
}

impl Representation {
  pub(crate) fn resolve(dense: Option<bool>) -> Self {
    if dense.unwrap_or_else(infer_dense) {
      Self::Dense
    } else {
      Self::Sparse
    }
  }
}

pub(crate) fn build_sparse_reconstruction(
  model: &GtrModel,
  rate: InferredRate,
  graph: &Graph,
  alphabet: Alphabet,
  node_inputs: NodeMap<NodeSeqInput>,
  branch_lengths: &EdgeMap<f64>,
  log: &dyn LogSink,
) -> Result<SparseReconstruction, Report> {
  let fitch = create_fitch_partition(graph, alphabet, node_inputs)?;
  sparse_from_fitch(model, rate, fitch, graph, branch_lengths, log)
}

fn sparse_from_fitch(
  model: &GtrModel,
  rate: InferredRate,
  fitch: PartitionFitch,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  log: &dyn LogSink,
) -> Result<SparseReconstruction, Report> {
  let gtr = fitch_gtr(model, rate, &fitch, graph, branch_lengths, log)?;
  let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
  Ok(SparseReconstruction::seeded(partition, gtr, node_states))
}

fn fitch_gtr(
  model: &GtrModel,
  rate: InferredRate,
  fitch: &PartitionFitch,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  log: &dyn LogSink,
) -> Result<GTR, Report> {
  if !model.is_infer() {
    return named_gtr(model, fitch.alphabet.name(), log);
  }
  let mut gtr = infer_gtr_fitch(fitch, graph, branch_lengths, log)?;
  rate.apply(&mut gtr)?;
  log_gtr(&gtr, model.name(), log)?;
  Ok(gtr)
}

fn named_gtr(model: &GtrModel, alphabet: AlphabetName, log: &dyn LogSink) -> Result<GTR, Report> {
  let gtr = model.build(alphabet)?;
  log_gtr(&gtr, model.name(), log)?;
  Ok(gtr)
}
