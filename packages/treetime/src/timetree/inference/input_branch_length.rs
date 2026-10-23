use crate::branch_lengths::{branch_lengths_or_zero, one_mutation};
use crate::gtr::get_gtr::{GtrModel, JC69Params, log_gtr};
use crate::partition::create::source_gtr;
use crate::partition::inferred_rate::InferredRate;
use crate::partition::resolve::model_alphabet;
use crate::partition::split::{PartitionData, SequenceInput, SplitInput, SplitPartition};
use crate::progress::LogSink;
use crate::timetree::inference::runner::GRID_POINTS;
use crate::timetree::pipeline::{SequenceLength, TimetreeSequences};
use argmin::core::{CostFunction, Error, Executor, TerminationReason, TerminationStatus};
use argmin::solver::brent::BrentRoot;
use eyre::Report;
use log::debug;
use ndarray::Array1;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{Distribution, DistributionFunction, NegLog};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_grid::{BoundaryBehavior, EdgeCellLaw, GridFn, MaxGridPoints, SoftTailLaw};
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_internal_error, make_internal_report, make_report};

const POISSON_MAX_MUTATION_LENGTH: f64 = 0.05;

const MIN_BRANCH_LENGTH_MUTATION_FRACTION: f64 = 1e-3;

const GAUSSIAN_TAIL_OFFSET: f64 = 100.0;

const POISSON_GRID_NEG_LOG_SPAN: f64 = 100.0;

const GRID_END_SOLVE_MAX_ITERS: u64 = 200;

const GRID_POINTS_PER_WIDTH: f64 = 10.0;

pub(crate) fn input_branch_model(
  sequence_length: usize,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  sequences: TimetreeSequences,
  log: &dyn LogSink,
) -> Result<InputBranchModel, Report> {
  let saturated_difference = match sequences {
    TimetreeSequences::Length(SequenceLength { model, alphabet, .. }) => {
      let model = if model.is_infer() {
        GtrModel::JC69(JC69Params::default())
      } else {
        model
      };
      let gtr = model.build(model_alphabet(&model, alphabet))?;
      log_gtr(&gtr, model.name(), log)?;
      saturated_difference(gtr.pi())
    },
    TimetreeSequences::Alignments(input) => partitions_saturated_difference(graph, branch_lengths, input, log)?,
  };
  debug!("Input branch lengths: saturated difference = {saturated_difference:.6}");
  Ok(InputBranchModel::new(sequence_length, saturated_difference))
}

pub(crate) fn saturated_difference(pi: &Array1<f64>) -> f64 {
  1.0 - pi.dot(pi)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct InputBranchModel {
  one_mutation: f64,
  saturated_difference: f64,
}

impl InputBranchModel {
  pub(crate) fn new(sequence_length: usize, saturated_difference: f64) -> Self {
    Self {
      one_mutation: one_mutation(sequence_length),
      saturated_difference,
    }
  }

  pub(crate) fn likelihood(self, mutation_length: f64) -> InputBranchLikelihood {
    let Self {
      one_mutation,
      saturated_difference,
    } = self;
    if mutation_length < POISSON_MAX_MUTATION_LENGTH {
      return InputBranchLikelihood::Poisson {
        mutation_length,
        one_mutation,
      };
    }
    let corrected_length = mutation_length + one_mutation;
    let growth = (corrected_length / saturated_difference).exp();
    let variance =
      saturated_difference * (growth - 1.0) * (growth - saturated_difference * (growth - 1.0)) * one_mutation;
    let tail_scale = (variance + MIN_BRANCH_LENGTH_MUTATION_FRACTION * one_mutation).sqrt();
    InputBranchLikelihood::Gaussian {
      mutation_length,
      variance,
      tail_scale,
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum InputBranchLikelihood {
  Poisson {
    mutation_length: f64,
    one_mutation: f64,
  },
  Gaussian {
    mutation_length: f64,
    variance: f64,
    tail_scale: f64,
  },
}

impl InputBranchLikelihood {
  pub(crate) fn neg_log_lh(self, branch_length: f64) -> f64 {
    match self {
      Self::Poisson {
        mutation_length,
        one_mutation,
      } => {
        let offset = MIN_BRANCH_LENGTH_MUTATION_FRACTION * one_mutation;
        (branch_length - mutation_length * (branch_length + offset).ln()) / one_mutation
      },
      Self::Gaussian {
        mutation_length,
        variance,
        tail_scale,
      } => {
        let deviation = mutation_length - branch_length;
        f64::min(
          0.5 * deviation * deviation / variance,
          GAUSSIAN_TAIL_OFFSET + deviation.abs() / tail_scale,
        )
      },
    }
  }

  pub(crate) fn distribution(
    self,
    branch_rate: f64,
    max_points: MaxGridPoints,
  ) -> Result<Arc<Distribution<NegLog>>, Report> {
    self.distribution_with_points_per_width(branch_rate, GRID_POINTS_PER_WIDTH, max_points)
  }

  pub(crate) fn distribution_with_points_per_width(
    self,
    branch_rate: f64,
    points_per_width: f64,
    max_points: MaxGridPoints,
  ) -> Result<Arc<Distribution<NegLog>>, Report> {
    debug_assert!(branch_rate > 0.0, "branch_rate must be positive, got {branch_rate:.6e}");
    let max_branch_length = self.grid_end()?;
    let time_range = (0.0, max_branch_length / branch_rate);
    let dx = self.width() / points_per_width;
    let n_points = max_points
      .point_count((max_branch_length / dx).ceil() + 1.0, time_range, dx / branch_rate)?
      .max(GRID_POINTS);

    let neg_log_lh =
      Array1::linspace(0.0, max_branch_length, n_points).mapv(|branch_length| self.neg_log_lh(branch_length));
    let peak = neg_log_lh.fold(f64::INFINITY, |min, &value| min.min(value));
    let f = GridFn::from_range_values(time_range, neg_log_lh.mapv(|value| value - peak))?;

    let right_tail = SoftTailLaw {
      slope: self.slope(max_branch_length) * branch_rate,
    };
    let distribution_fn = DistributionFunction::from_grid_fn(f)
      .with_left_extrap(self.start(branch_rate))?
      .with_right_extrap(BoundaryBehavior::Linear(right_tail))?;
    Ok(Arc::new(Distribution::Function(distribution_fn)))
  }

  fn start(self, branch_rate: f64) -> BoundaryBehavior {
    match self {
      Self::Poisson {
        mutation_length,
        one_mutation,
      } => BoundaryBehavior::HardEdgeCell(EdgeCellLaw {
        t_edge: 0.0,
        offset: MIN_BRANCH_LENGTH_MUTATION_FRACTION * one_mutation / branch_rate,
        b: mutation_length / one_mutation,
      }),
      Self::Gaussian { .. } => BoundaryBehavior::Hard,
    }
  }

  fn grid_end(self) -> Result<f64, Report> {
    match self {
      Self::Poisson {
        mutation_length,
        one_mutation,
      } => {
        let peak = self.peak();
        let upper = peak + 2.0 * (mutation_length + POISSON_GRID_NEG_LOG_SPAN * one_mutation);
        let level = PoissonLevel {
          likelihood: self,
          level: self.neg_log_lh(peak) + POISSON_GRID_NEG_LOG_SPAN,
        };
        let result = Executor::new(level, BrentRoot::new(peak, upper, 0.0))
          .configure(|cfg| cfg.max_iters(GRID_END_SOLVE_MAX_ITERS))
          .run()
          .map_err(|e| make_report!("Grid end of the input branch-length likelihood not found: {e}"))?;
        let status = &result.state.termination_status;
        if *status != TerminationStatus::Terminated(TerminationReason::SolverConverged) {
          return make_internal_error!(
            "Grid end solve of the input branch-length likelihood over [{peak}, {upper}] did not converge: {status:?}"
          );
        }
        result.state.best_param.ok_or_else(|| {
          make_internal_report!("Grid end solve of the input branch-length likelihood returned no length")
        })
      },
      Self::Gaussian {
        mutation_length,
        variance,
        tail_scale,
      } => {
        let inverse_scale = tail_scale.recip();
        let transition =
          variance * (inverse_scale + (inverse_scale * inverse_scale + 2.0 * GAUSSIAN_TAIL_OFFSET / variance).sqrt());
        Ok(mutation_length + transition)
      },
    }
  }

  fn slope(self, branch_length: f64) -> f64 {
    match self {
      Self::Poisson {
        mutation_length,
        one_mutation,
      } => {
        let offset = MIN_BRANCH_LENGTH_MUTATION_FRACTION * one_mutation;
        (1.0 - mutation_length / (branch_length + offset)) / one_mutation
      },
      Self::Gaussian { tail_scale, .. } => tail_scale.recip(),
    }
  }

  fn peak(self) -> f64 {
    match self {
      Self::Poisson {
        mutation_length,
        one_mutation,
      } => (mutation_length - MIN_BRANCH_LENGTH_MUTATION_FRACTION * one_mutation).max(0.0),
      Self::Gaussian { mutation_length, .. } => mutation_length.max(0.0),
    }
  }

  fn width(self) -> f64 {
    match self {
      Self::Poisson {
        mutation_length,
        one_mutation,
      } => (mutation_length.max(one_mutation) * one_mutation).sqrt(),
      Self::Gaussian { variance, .. } => variance.sqrt(),
    }
  }
}

fn partitions_saturated_difference(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  input: SequenceInput,
  log: &dyn LogSink,
) -> Result<f64, Report> {
  let SplitInput { partitions, .. } = input.split()?;
  let rate = InferredRate::for_partition_count(partitions.len());
  let branch_lengths = branch_lengths_or_zero(branch_lengths);
  let total_sites: f64 = partitions
    .iter()
    .map(|partition| partition.data.layout.columns.len().to_f64_exact())
    .sum();
  let weighted: f64 = partitions
    .into_iter()
    .map(|SplitPartition { data, model }| -> Result<f64, Report> {
      let PartitionData {
        layout,
        alphabet,
        nodes,
      } = data;
      let gtr = source_gtr(&model, rate, graph, alphabet, nodes, &branch_lengths, log)?;
      Ok(layout.columns.len().to_f64_exact() * saturated_difference(gtr.pi()))
    })
    .sum::<Result<f64, Report>>()?;
  Ok(weighted / total_sites)
}

struct PoissonLevel {
  likelihood: InputBranchLikelihood,
  level: f64,
}

impl CostFunction for PoissonLevel {
  type Param = f64;
  type Output = f64;

  fn cost(&self, branch_length: &f64) -> Result<f64, Error> {
    Ok(self.likelihood.neg_log_lh(*branch_length) - self.level)
  }
}
