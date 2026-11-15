use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::gtr::gamma::{GammaRates, INITIAL_GAMMA_SHAPE};
use crate::gtr::get_gtr::{GtrModel, log_gtr};
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::{GtrSharing, InferGtrOptions, infer_gtr_group};
use crate::gtr::partition_model::PartitionModel;
use crate::gtr::shared_params::SharedParams;
use crate::partition::algo::infer_dense::infer_dense;
use crate::partition::fitch::gtr_inference::{get_mutation_counts_fitch, gtr_from_inference, infer_gtr_fitch};
use crate::partition::fitch::partition::PartitionFitch;
use crate::partition::fitch::passes::{create_fitch_partition, create_inherited_fitch_partition};
use crate::partition::inferred_rate::InferredRate;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction, SparseReconstruction};
use crate::partition::resolve::PartitionLayout;
use crate::partition::split::ModelSettings;
use crate::progress::LogSink;
use crate::seq::alignment::NodeSeqInput;
use eyre::Report;
use itertools::{Itertools, izip};
use std::sync::Arc;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::scope::{EdgeScopes, Inheritance};
use treetime_utils::{make_error, make_internal_report};

#[expect(
  clippy::too_many_arguments,
  reason = "the partition takes its representation, model, rate rule and gamma category count, and its data (graph, alphabet, leaf inputs, branch lengths) and the log sink as separate values; no type holds them together"
)]
pub(crate) fn build_marginal_partition(
  representation: Representation,
  model: &GtrModel,
  rate: InferredRate,
  gamma_categories: usize,
  graph: &Graph,
  alphabet: Alphabet,
  node_inputs: NodeMap<NodeSeqInput>,
  branch_lengths: &EdgeMap<f64>,
  log: &dyn LogSink,
) -> Result<MarginalReconstruction, Report> {
  let input = PartitionInput {
    model,
    rate,
    inheritance: Inheritance::ALL,
    alphabet,
    nodes: node_inputs,
  };
  let models = ModelSettings {
    gamma_categories,
    shared: SharedParams::default(),
  };
  build_marginal_partitions(representation, vec![input], models, graph, branch_lengths, log)?
    .pop()
    .ok_or_else(|| make_internal_report!("Building one partition returned no reconstruction"))
}

pub(crate) fn build_marginal_partitions(
  representation: Representation,
  inputs: Vec<PartitionInput<'_>>,
  models: ModelSettings,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  log: &dyn LogSink,
) -> Result<Vec<MarginalReconstruction>, Report> {
  let gamma = GammaRates::new(models.gamma_categories, INITIAL_GAMMA_SHAPE)?;
  let prepared = inputs
    .into_iter()
    .map(|input| PreparedPartition::new(representation, input, graph))
    .collect::<Result<Vec<_>, Report>>()?;
  let gtrs = fit_models(&prepared, models.shared, graph, branch_lengths, log)?;
  izip!(prepared, gtrs)
    .map(|(prepared, gtr)| {
      let shared_shape = models.shared.gamma_shape && prepared.model.is_infer();
      let model = PartitionModel::new(gtr, gamma.clone()).with_shared_shape(shared_shape);
      prepared.seed(model, graph)
    })
    .collect()
}

pub(crate) struct PartitionInput<'a> {
  pub(crate) model: &'a GtrModel,
  pub(crate) rate: InferredRate,
  pub(crate) inheritance: Inheritance,
  pub(crate) alphabet: Alphabet,
  pub(crate) nodes: NodeMap<NodeSeqInput>,
}

pub(crate) fn segment_inheritances<'l>(
  scopes: Option<&Arc<EdgeScopes>>,
  alignment_count: usize,
  layouts: impl ExactSizeIterator<Item = &'l PartitionLayout>,
) -> Result<Vec<Inheritance>, Report> {
  let Some(scopes) = scopes else {
    return Ok(vec![Inheritance::ALL; layouts.len()]);
  };
  if scopes.segment_count() != alignment_count {
    return make_error!(
      "The graph has {} segment trees, but the run has {} alignments: give one alignment per segment tree, in the same order",
      scopes.segment_count(),
      alignment_count
    );
  }
  layouts
    .map(|layout| Inheritance::segment(Arc::clone(scopes), layout.alignment - 1))
    .collect()
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

fn fit_models(
  prepared: &[PreparedPartition<'_>],
  shared: SharedParams,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  log: &dyn LogSink,
) -> Result<Vec<GTR>, Report> {
  let inferred = prepared
    .iter()
    .enumerate()
    .filter_map(|(index, prepared)| Some((index, prepared.inferred_fitch()?)))
    .collect_vec();
  let sharing = GtrSharing {
    exchange_rates: shared.exchange_rates,
    frequencies: shared.frequencies,
  };
  let mut fits: Vec<Option<GTR>> = vec![None; prepared.len()];
  if inferred.len() >= 2 && shared.shares_rate_matrix() {
    let counts = inferred
      .iter()
      .map(|(_, fitch)| get_mutation_counts_fitch(graph, fitch, branch_lengths))
      .collect::<Result<Vec<_>, Report>>()?;
    let counts = counts.iter().collect_vec();
    let results = infer_gtr_group(&counts, sharing, &InferGtrOptions::default(), log)?;
    for ((index, fitch), result) in izip!(&inferred, results) {
      fits[*index] = Some(gtr_from_inference(fitch, result)?);
    }
  } else {
    for (index, fitch) in &inferred {
      fits[*index] = Some(infer_gtr_fitch(fitch, graph, branch_lengths, log)?);
    }
  }
  izip!(prepared, fits)
    .map(|(prepared, fit)| match fit {
      Some(mut gtr) => {
        prepared.rate.apply(&mut gtr)?;
        log_gtr(&gtr, prepared.model.name(), log)?;
        Ok(gtr)
      },
      None => named_gtr(prepared.model, prepared.alphabet_name(), log),
    })
    .collect()
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

struct PreparedPartition<'a> {
  model: &'a GtrModel,
  rate: InferredRate,
  data: PreparedData,
}

enum PreparedData {
  Sparse(PartitionFitch),
  DenseFitch(PartitionFitch, NodeMap<NodeSeqInput>),
  Dense(Alphabet, NodeMap<NodeSeqInput>, Inheritance),
}

impl<'a> PreparedPartition<'a> {
  fn new(representation: Representation, input: PartitionInput<'a>, graph: &Graph) -> Result<Self, Report> {
    let PartitionInput {
      model,
      rate,
      inheritance,
      alphabet,
      nodes,
    } = input;
    let data = match representation {
      Representation::Sparse => {
        PreparedData::Sparse(create_inherited_fitch_partition(graph, inheritance, alphabet, nodes)?)
      },
      Representation::Dense if model.is_infer() => PreparedData::DenseFitch(
        create_inherited_fitch_partition(graph, inheritance, alphabet, nodes.clone())?,
        nodes,
      ),
      Representation::Dense => PreparedData::Dense(alphabet, nodes, inheritance),
    };
    Ok(Self { model, rate, data })
  }

  fn inferred_fitch(&self) -> Option<&PartitionFitch> {
    if !self.model.is_infer() {
      return None;
    }
    match &self.data {
      PreparedData::Sparse(fitch) | PreparedData::DenseFitch(fitch, _) => Some(fitch),
      PreparedData::Dense(..) => None,
    }
  }

  fn alphabet_name(&self) -> AlphabetName {
    match &self.data {
      PreparedData::Sparse(fitch) | PreparedData::DenseFitch(fitch, _) => fitch.alphabet.name(),
      PreparedData::Dense(alphabet, ..) => alphabet.name(),
    }
  }

  fn seed(self, model: PartitionModel, graph: &Graph) -> Result<MarginalReconstruction, Report> {
    Ok(match self.data {
      PreparedData::Sparse(fitch) => {
        let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
        MarginalReconstruction::Sparse(SparseReconstruction::seeded(partition, model, node_states))
      },
      PreparedData::DenseFitch(fitch, nodes) => {
        let partition = fitch.into_marginal_dense(graph, &nodes)?;
        MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, model))
      },
      PreparedData::Dense(alphabet, nodes, inheritance) => {
        let partition = PartitionMarginalDense::new(alphabet, graph, &nodes)?.with_inheritance(inheritance);
        MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, model))
      },
    })
  }
}
