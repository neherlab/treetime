use crate::coalescent::mean_field::CoalescentPrior;
use crate::test_utils::reference_grid::{Ends, ExactFunction, GridWindow, convolve_windows, to_isize, to_isize_floor};
use crate::timetree::inference::result::BranchLikelihood;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use ordered_float::OrderedFloat;
use ndarray::Array1;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::{BoundaryBehavior, Distribution, NegLog};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_report};

const WINDOW_FRACTION: f64 = 1e-15;
const WIDE_WINDOW_FRACTION: f64 = 1e-18;
const DOMAIN_TOLERANCE: f64 = 1e-9;
const MIN_OBSERVED_ORDER: f64 = 1.8;
const CONVERGED_DIFFERENCE: f64 = 1e-8;
const TIP_GRID_TOLERANCE: f64 = 1e-6;
const MAX_BRANCH_POINTS: usize = 10_000_000;
const REFINEMENTS: [f64; 3] = [1.0, 0.5, 0.25];

pub(crate) fn reference_log_z(input: &ReferenceInput<'_>) -> Result<ReferenceLogZ, Report> {
  let spacing = input.spacing;
  if !(spacing.is_finite() && spacing > 0.0) {
    return make_error!("Reference grid spacing must be positive and finite, got {spacing}");
  }
  let levels: Vec<ReferenceLevel> = REFINEMENTS
    .iter()
    .map(|refinement| reference_level(input, spacing * refinement, WINDOW_FRACTION))
    .collect::<Result<_, Report>>()?;
  let [coarse, middle, fine] = [levels[0], levels[1], levels[2]];
  let convergence = convergence(coarse.log_z, middle.log_z, fine.log_z)?;
  let extrapolate = |middle: f64, fine: f64| match convergence {
    Convergence::SecondOrder(_) => richardson(middle, fine),
    Convergence::Converged(_) => fine,
  };
  let log_z = extrapolate(middle.log_z, fine.log_z);
  let discretization_error = match convergence {
    Convergence::SecondOrder(_) => (log_z - fine.log_z).abs(),
    Convergence::Converged(difference) => difference,
  };
  let wide = reference_level(input, spacing, WIDE_WINDOW_FRACTION)?;
  let domain_error = (wide.log_z - coarse.log_z).abs();
  if !(domain_error < DOMAIN_TOLERANCE) {
    return make_error!(
      "Reference log Z changes by {domain_error:e} when the windows widen from {WINDOW_FRACTION:e} to {WIDE_WINDOW_FRACTION:e} of their peaks, above the limit {DOMAIN_TOLERANCE:e}"
    );
  }
  Ok(ReferenceLogZ {
    log_z,
    discretization_error,
    domain_error,
    convergence,
    root_time_mean: extrapolate(middle.root_time_mean, fine.root_time_mean),
    root_time_variance: extrapolate(middle.root_time_variance, fine.root_time_variance),
    levels,
  })
}

pub(crate) struct ReferenceInput<'a> {
  pub graph: &'a Graph,
  pub tip_dates: &'a BTreeMap<GraphNodeKey, f64>,
  pub branches: &'a BTreeMap<GraphEdgeKey, BranchTime>,
  pub factors: NodeFactors<'a>,
  pub spacing: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ReferenceLogZ {
  pub log_z: f64,
  pub discretization_error: f64,
  pub domain_error: f64,
  pub convergence: Convergence,
  pub root_time_mean: f64,
  pub root_time_variance: f64,
  pub levels: Vec<ReferenceLevel>,
}

impl ReferenceLogZ {
  pub(crate) fn error(&self) -> f64 {
    self.discretization_error + self.domain_error
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Convergence {
  SecondOrder(f64),
  Converged(f64),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ReferenceLevel {
  pub spacing: f64,
  pub log_z: f64,
  pub root_time_mean: f64,
  pub root_time_variance: f64,
}

#[derive(Clone)]
pub(crate) enum BranchTime {
  Distribution(Arc<Distribution<NegLog>>),
  Analytic(AnalyticBranch),
}

impl BranchTime {
  pub(crate) fn from_likelihoods(
    branches: &BTreeMap<GraphEdgeKey, BranchLikelihood>,
  ) -> Result<BTreeMap<GraphEdgeKey, Self>, Report> {
    branches
      .iter()
      .map(|(&edge, branch)| {
        let distribution = branch
          .distribution
          .clone()
          .ok_or_else(|| make_report!("Edge {edge} has no branch-time distribution"))?;
        Ok((edge, Self::Distribution(distribution)))
      })
      .collect()
  }

  fn support(&self) -> Result<BranchSupport, Report> {
    match self {
      Self::Distribution(distribution) => {
        let (Some((lo, hi)), Some((_, tabulated_end))) = (distribution.support_bounds(), distribution.time_bounds())
        else {
          return make_error!("A branch-time distribution without support cannot enter the reference");
        };
        Ok(BranchSupport {
          lo: lo.max(0.0),
          hi,
          monotone_after: tabulated_end,
        })
      },
      Self::Analytic(analytic) => Ok(BranchSupport {
        lo: 0.0,
        hi: f64::INFINITY,
        monotone_after: analytic.monotone_after,
      }),
    }
  }

  fn knots(&self) -> Vec<f64> {
    match self {
      Self::Distribution(distribution) => distribution.t().to_vec(),
      Self::Analytic(_) => vec![],
    }
  }

  fn inner(&self) -> Result<Option<InnerLaw>, Report> {
    let Self::Distribution(distribution) = self else {
      return Ok(None);
    };
    let Distribution::Function(function) = distribution.as_ref() else {
      return Ok(None);
    };
    let BoundaryBehavior::HardApproach(law) = function.left_extrap() else {
      return Ok(None);
    };
    let knots = distribution.t();
    if knots.len() < 2 || law.t_hard >= knots[0] {
      return Ok(None);
    }
    let (x0, x1) = (knots[0], knots[1]);
    let (y0, y1) = (self.neg_log(x0)?, self.neg_log(x1)?);
    Ok(Some(InnerLaw {
      to: x0,
      neg_log_at_to: y0,
      slope: (y1 - y0) / (x1 - x0),
    }))
  }

  fn neg_log_outer(&self, inner: Option<InnerLaw>, tau: f64) -> Result<f64, Report> {
    match inner {
      Some(inner) if tau < inner.to => Ok(inner.neg_log_at_to + inner.slope * (tau - inner.to)),
      _ => self.neg_log(tau),
    }
  }

  fn neg_log(&self, tau: f64) -> Result<f64, Report> {
    let value = match self {
      Self::Distribution(distribution) => distribution.eval(tau)?,
      Self::Analytic(analytic) => (analytic.neg_log)(tau),
    };
    if value.is_nan() || (value.is_infinite() && value < 0.0) {
      return make_error!("Branch-time function gives the negative log density {value} at time {tau}");
    }
    Ok(value)
  }
}

#[derive(Clone)]
pub(crate) struct AnalyticBranch {
  pub neg_log: Arc<dyn Fn(f64) -> f64 + Send + Sync>,
  pub monotone_after: f64,
}

#[derive(Clone, Copy)]
pub(crate) enum NodeFactors<'a> {
  None,
  Prior(&'a CoalescentPrior),
  Analytic(&'a (dyn Fn(GraphNodeKey, usize, f64) -> f64 + Sync)),
}

impl NodeFactors<'_> {
  fn kinks(self) -> Vec<f64> {
    match self {
      Self::Prior(prior) => prior
        .expectation
        .grid
        .iter()
        .tuple_windows()
        .filter(|(left, right)| left == right)
        .map(|(&time, _)| time)
        .collect(),
      Self::None | Self::Analytic(_) => vec![],
    }
  }

  fn neg_log(self, key: GraphNodeKey, n_children: usize, time: f64) -> Result<f64, Report> {
    match self {
      Self::None => Ok(0.0),
      Self::Prior(prior) => prior.node_contribution(key, n_children, time),
      Self::Analytic(factor) => Ok(factor(key, n_children, time)),
    }
  }
}

#[derive(Clone, Copy)]
struct InnerLaw {
  to: f64,
  neg_log_at_to: f64,
  slope: f64,
}

struct BranchSupport {
  lo: f64,
  hi: f64,
  monotone_after: f64,
}

fn convergence(coarse: f64, middle: f64, fine: f64) -> Result<Convergence, Report> {
  let (first, second) = (coarse - middle, middle - fine);
  let largest = first.abs().max(second.abs());
  if largest <= CONVERGED_DIFFERENCE {
    return Ok(Convergence::Converged(largest));
  }
  let order = (first / second).log2();
  if !(order >= MIN_OBSERVED_ORDER) {
    return make_error!(
      "Reference log Z converges at the observed order {order} below {MIN_OBSERVED_ORDER}: {coarse}, {middle}, {fine} at the spacings h, h/2, h/4"
    );
  }
  Ok(Convergence::SecondOrder(order))
}

fn richardson(middle: f64, fine: f64) -> f64 {
  (4.0 * fine - middle) / 3.0
}

fn reference_level(input: &ReferenceInput<'_>, spacing: f64, fraction: f64) -> Result<ReferenceLevel, Report> {
  let schedule = GraphSchedule::new(input.graph)?;
  schedule.require_tree()?;
  let reference = input.tip_dates.values().copied().fold(f64::NEG_INFINITY, f64::max);
  let time = |index: f64| reference - index * spacing;
  let factor_breaks = sorted(input.factors.kinks().iter().map(|&kink| (reference - kink) / spacing));
  let mut subtrees: BTreeMap<GraphNodeKey, Subtree> = BTreeMap::new();
  for &key in schedule.postorder() {
    let children = schedule.children(key);
    let subtree = if children.is_empty() {
      leaf_subtree(input, key, reference, spacing)?
    } else {
      if input.tip_dates.contains_key(&key) {
        return make_error!("Internal node {key} has a date; the reference dates only tips");
      }
      let messages: Vec<ChildMessage<'_>> = children
        .iter()
        .map(|&(child, edge)| {
          let subtree = subtrees
            .remove(&child)
            .ok_or_else(|| make_report!("Node {child} has no subtree when its parent {key} needs it"))?;
          let branch = input
            .branches
            .get(&edge)
            .ok_or_else(|| make_report!("Edge {edge} has no branch-time function"))?;
          child_message(subtree, branch, spacing, fraction)
            .wrap_err_with(|| format!("When sending the message along edge {edge}"))
        })
        .collect::<Result<_, Report>>()?;
      let factor = |index: f64| input.factors.neg_log(key, children.len(), time(index));
      Subtree::Grid(node_window(key, &messages, factor, &factor_breaks, spacing, fraction)?)
    };
    subtrees.insert(key, subtree);
  }
  let root = input.graph.root_key()?;
  let Some(Subtree::Grid(root_window)) = subtrees.remove(&root) else {
    return make_error!("The reference needs a root with children");
  };
  root_level(&root_window, spacing, time)
}

fn leaf_subtree(
  input: &ReferenceInput<'_>,
  key: GraphNodeKey,
  reference: f64,
  spacing: f64,
) -> Result<Subtree, Report> {
  let date = *input
    .tip_dates
    .get(&key)
    .ok_or_else(|| make_report!("Tip {key} has no date"))?;
  let position = (reference - date) / spacing;
  let index = position.round();
  if (position - index).abs() > TIP_GRID_TOLERANCE {
    return make_error!(
      "Tip date {date} lies {position} grid steps of {spacing} before the youngest tip {reference}; the spacing must divide every difference of tip dates"
    );
  }
  Ok(Subtree::Leaf {
    index: to_isize_floor(index),
    log_factor: -input.factors.neg_log(key, 0, date)?,
  })
}

enum Subtree {
  Leaf { index: isize, log_factor: f64 },
  Grid(GridWindow),
}

struct ChildMessage<'a> {
  window: GridWindow,
  leaf: Option<LeafSource<'a>>,
}

struct LeafSource<'a> {
  index: isize,
  log_factor: f64,
  branch: &'a BranchTime,
  inner: Option<InnerLaw>,
}

impl LeafSource<'_> {
  fn log_values(&self, x: f64, spacing: f64) -> Result<(f64, f64), Report> {
    let tau = (x - self.index.to_f64_exact()) * spacing;
    let true_log = -self.branch.neg_log(tau)? + self.log_factor;
    let outer_log = -self.branch.neg_log_outer(self.inner, tau)? + self.log_factor;
    Ok((true_log, outer_log))
  }

  fn breaks(&self, spacing: f64) -> impl Iterator<Item = f64> {
    let index = self.index.to_f64_exact();
    self.branch.knots().into_iter().map(move |knot| index + knot / spacing)
  }

  fn singular(&self, spacing: f64) -> Result<Option<f64>, Report> {
    self
      .inner
      .map(|_| Ok(self.index.to_f64_exact() + self.branch.support()?.lo / spacing))
      .transpose()
  }
}

fn child_message<'a>(
  subtree: Subtree,
  branch: &'a BranchTime,
  spacing: f64,
  fraction: f64,
) -> Result<ChildMessage<'a>, Report> {
  let window = branch_window(branch, spacing, fraction)?;
  Ok(match subtree {
    Subtree::Leaf { index, log_factor } => ChildMessage {
      window: window.shifted(index, log_factor),
      leaf: Some(LeafSource {
        index,
        log_factor,
        branch,
        inner: branch.inner()?,
      }),
    },
    Subtree::Grid(subtree) => ChildMessage {
      window: convolve_windows(&subtree, &window, spacing, fraction)?,
      leaf: None,
    },
  })
}

fn branch_window(branch: &BranchTime, spacing: f64, fraction: f64) -> Result<GridWindow, Report> {
  let BranchSupport { lo, hi, monotone_after } = branch.support()?;
  let inner = branch.inner()?;
  let first = to_isize_floor((lo / spacing).ceil());
  let cutoff = fraction.ln();
  let mut log_values = vec![];
  let mut peak = f64::NEG_INFINITY;
  let mut tau = first.to_f64_exact() * spacing;
  while tau <= hi {
    let log_value = -branch.neg_log_outer(inner, tau)?;
    peak = peak.max(log_value);
    log_values.push(log_value);
    if tau >= monotone_after && log_value < peak + cutoff {
      break;
    }
    if log_values.len() >= MAX_BRANCH_POINTS {
      return make_error!(
        "Branch-time function stays above {fraction:e} of its peak over {MAX_BRANCH_POINTS} grid points of spacing {spacing}"
      );
    }
    tau = (first + to_isize(log_values.len())).to_f64_exact() * spacing;
  }
  let old = (tau > hi)
    .then(|| -> Result<(f64, f64), Report> { Ok((hi / spacing, -branch.neg_log(hi)?)) })
    .transpose()?;
  let log_value = |x: f64| -> Result<f64, Report> { Ok(-branch.neg_log(x * spacing)?) };
  let breaks = sorted(branch.knots().into_iter().map(|knot| knot / spacing));
  let singular = inner.map(|_| lo / spacing).into_iter().collect_vec();
  let ends = Ends {
    young: Some((lo / spacing, -branch.neg_log_outer(inner, lo)?)),
    old,
    exact: Some(ExactFunction {
      log_value: &log_value,
      breaks: &breaks,
      singular: &singular,
      log_linear: matches!(branch, BranchTime::Distribution(_)),
    }),
  };
  GridWindow::from_log_values(first, &Array1::from_vec(log_values), ends, fraction)
}

fn node_window(
  key: GraphNodeKey,
  messages: &[ChildMessage<'_>],
  factor: impl Fn(f64) -> Result<f64, Report>,
  factor_breaks: &[f64],
  spacing: f64,
  fraction: f64,
) -> Result<GridWindow, Report> {
  let windows = messages.iter().map(|message| &message.window).collect_vec();
  let young = windows
    .iter()
    .map(|window| window.young_end())
    .fold(f64::NEG_INFINITY, f64::max);
  let old = windows
    .iter()
    .map(|window| window.old_end())
    .fold(f64::INFINITY, f64::min);
  let (first, last) = (to_isize_floor(young.ceil()), to_isize_floor(old.floor()));
  if last < first {
    return make_error!("The child messages of node {key} do not overlap");
  }
  let log_values: Array1<f64> = (first..=last)
    .map(|index| {
      let messages: f64 = windows.iter().map(|window| window.log_sample(index)).sum();
      Ok(messages - factor(index.to_f64_exact())?)
    })
    .collect::<Result<Vec<f64>, Report>>()?
    .into();
  let log_values_at = |x: f64| -> Result<(f64, f64), Report> {
    let (true_log, outer_log) = messages.iter().try_fold((0.0, 0.0), |(true_sum, outer_sum), message| {
      let (true_log, outer_log) = match &message.leaf {
        Some(leaf) => leaf.log_values(x, spacing)?,
        None => (message.window.log_at(x), message.window.log_at(x)),
      };
      Ok::<_, Report>((true_sum + true_log, outer_sum + outer_log))
    })?;
    let factor = factor(x)?;
    Ok((true_log - factor, outer_log - factor))
  };
  let young_edged = windows
    .iter()
    .filter(|window| window.young_end().total_cmp(&young).is_eq())
    .all(|window| window.has_young_edge());
  let old_edged = windows
    .iter()
    .filter(|window| window.old_end().total_cmp(&old).is_eq())
    .all(|window| window.has_old_edge());
  let leaves = messages.iter().filter_map(|message| message.leaf.as_ref()).collect_vec();
  let breaks = sorted(
    leaves
      .iter()
      .flat_map(|leaf| leaf.breaks(spacing))
      .chain(factor_breaks.iter().copied()),
  );
  let singular = leaves
    .iter()
    .map(|leaf| leaf.singular(spacing))
    .filter_map(Result::transpose)
    .collect::<Result<Vec<f64>, Report>>()?;
  let log_value = |x: f64| -> Result<f64, Report> { Ok(log_values_at(x)?.0) };
  let ends = Ends {
    young: young_edged
      .then(|| -> Result<(f64, f64), Report> { Ok((young, log_values_at(young)?.1)) })
      .transpose()?,
    old: old_edged
      .then(|| -> Result<(f64, f64), Report> { Ok((old, log_values_at(old)?.0)) })
      .transpose()?,
    exact: Some(ExactFunction {
      log_value: &log_value,
      breaks: &breaks,
      singular: &singular,
      log_linear: false,
    }),
  };
  GridWindow::from_log_values(first, &log_values, ends, fraction)
    .wrap_err_with(|| format!("When forming the subtree function of node {key}"))
}

fn root_level(root: &GridWindow, spacing: f64, time: impl Fn(f64) -> f64) -> Result<ReferenceLevel, Report> {
  let mass = root.integral(|_| 1.0);
  if !(mass.is_finite() && mass > 0.0) {
    return make_error!("The root function has the integral {mass}");
  }
  let mean = root.integral(&time) / mass;
  let variance = root.integral(|index| (time(index) - mean).powi(2)) / mass;
  Ok(ReferenceLevel {
    spacing,
    log_z: root.log_amp() + spacing.ln() + mass.ln(),
    root_time_mean: mean,
    root_time_variance: variance,
  })
}

fn sorted(values: impl Iterator<Item = f64>) -> Vec<f64> {
  values.sorted_by_key(|&value| OrderedFloat(value)).dedup().collect()
}
