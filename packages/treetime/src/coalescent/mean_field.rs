use crate::clock::gaussian::NodeTimeGaussian;
use crate::timetree::inference::result::{NodePosterior, NodeTimes, TimeInference};
use crate::timetree::max_posterior_region::{PosteriorSurvival, SurvivalLevel, posterior_survival};
use eyre::{Report, WrapErr};
use itertools::{Itertools, izip};
use ndarray::Array1;
use ordered_float::OrderedFloat;
use statrs::distribution::{ContinuousCDF, Normal};
use std::collections::BTreeMap;
use std::iter::once;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNodeKey, Node};
use treetime_graph::schedule::GraphSchedule;
use treetime_grid::Side;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_internal_error};

const EXPECTATION_GRID_POINTS: usize = 4001;
const SURVIVAL_TAIL_LEVELS: [f64; 9] = [1e-9, 1e-7, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 2e-2, 3.5e-2];
const INTERIOR_LEVELS: usize = 40;
const MERGER_RATE_STRIDE: usize = 10;
const TRUNCATED_NORMAL_SPAN: f64 = 8.0;

#[derive(Clone, Debug)]
pub(crate) struct CoalescentPrior {
  expectation: Arc<LineageExpectation>,
  tc: Distribution,
  cell_rate: Vec<f64>,
  inverse_tc_integral: Vec<f64>,
  expected_integral: Vec<f64>,
  own_integral: BTreeMap<GraphNodeKey, OwnIntegral>,
  rate_before: f64,
  rate_after: f64,
  log_normalization: f64,
}

impl CoalescentPrior {
  pub(crate) fn new(expectation: &Arc<LineageExpectation>, tc: &Distribution) -> Result<Self, Report> {
    let grid = &expectation.grid;
    let inverse_tc = |t: f64| -> Result<f64, Report> {
      let value = tc
        .eval(t)
        .wrap_err_with(|| format!("When evaluating coalescent Tc at calendar time {t:.6e}"))?;
      if !value.is_finite() || value <= 0.0 {
        return make_error!("Coalescent Tc must be finite and positive at calendar time {t:.6e}, got {value:.6e}");
      }
      Ok(1.0 / value)
    };
    let cell_rate: Vec<f64> = grid
      .iter()
      .tuple_windows()
      .map(|(&t0, &t1)| inverse_tc(f64::midpoint(t0, t1)))
      .try_collect()?;
    let mut inverse_tc_integral = vec![0.0; grid.len()];
    let mut expected_integral = vec![0.0; grid.len()];
    let mut pair_integral = 0.0;
    for i in 1..grid.len() {
      let dt = grid[i] - grid[i - 1];
      let rate = cell_rate[i - 1];
      inverse_tc_integral[i] = inverse_tc_integral[i - 1] + dt * rate;
      expected_integral[i] =
        expected_integral[i - 1] + dt * rate * f64::midpoint(expectation.expected[i - 1], expectation.expected[i]);
      pair_integral += dt * rate * f64::midpoint(expectation.expected_pairs[i - 1], expectation.expected_pairs[i]);
    }
    let rate_before = inverse_tc(grid[0])?;
    let rate_after = inverse_tc(grid[grid.len() - 1])?;
    let rates = CellRates {
      grid,
      inside: &cell_rate,
      before: rate_before,
      after: rate_after,
    };
    let own_integral = expectation
      .nodes
      .iter()
      .map(|(key, node)| (*key, OwnIntegral::new(&node.survival, &rates)))
      .collect();
    let prior = Self {
      expectation: Arc::clone(expectation),
      tc: tc.clone(),
      cell_rate,
      inverse_tc_integral,
      expected_integral,
      own_integral,
      rate_before,
      rate_after,
      log_normalization: pair_integral,
    };
    let own_offset: f64 = expectation
      .nodes
      .values()
      .map(|node| node.weight.powi(2) * prior.inverse_tc_integral_at(node.survival.first()))
      .sum();
    Ok(Self {
      log_normalization: pair_integral + own_offset,
      ..prior
    })
  }

  pub(crate) fn log_normalization(&self) -> f64 {
    self.log_normalization
  }

  pub(crate) fn with_topology(
    &self,
    graph: &Graph,
    merger_times: &BTreeMap<GraphNodeKey, f64>,
  ) -> Result<Self, Report> {
    Self::new(&self.expectation.with_topology(graph, merger_times)?, &self.tc)
  }

  pub(crate) fn time_resolution(&self) -> f64 {
    self.expectation.step
  }

  pub(crate) fn node_contribution(&self, key: GraphNodeKey, n_children: usize, x: f64) -> Result<f64, Report> {
    let weight = 1.0 - n_children.to_f64_exact();
    let n_mergers = n_children.saturating_sub(1).to_f64_exact();
    let own = self.own_integral_at(key, x);
    let own_weight = self.expectation.nodes.get(&key).map_or(0.0, |node| node.weight);
    let merger = if n_mergers > 0.0 {
      let tc = self
        .tc
        .eval(x)
        .wrap_err_with(|| format!("When evaluating coalescent Tc at calendar time {x:.6e}"))?;
      n_mergers * tc.ln()
    } else {
      0.0
    };
    Ok(
      0.5 * weight * (weight - 1.0) * self.inverse_tc_integral_at(x) + weight * self.expected_integral_at(x)
        - weight * own_weight * own
        + merger,
    )
  }

  pub(crate) fn node_contribution_outward_slope(
    &self,
    key: GraphNodeKey,
    n_children: usize,
    side: Side,
    edge: f64,
  ) -> f64 {
    let weight = 1.0 - n_children.to_f64_exact();
    let own_weight = self.expectation.nodes.get(&key).map_or(0.0, |node| node.weight);
    let pairs = 0.5 * weight * (weight - 1.0);
    let grid = &self.expectation.grid;
    let last = grid.len() - 1;
    let inverse_tc = match piece_slope(grid, &self.inverse_tc_integral, side, edge) {
      Piece::Before => self.rate_before,
      Piece::Inside(slope) => slope,
      Piece::After => self.rate_after,
    };
    let expected = match piece_slope(grid, &self.expected_integral, side, edge) {
      Piece::Before => self.rate_before * self.expectation.expected_before,
      Piece::Inside(slope) => slope,
      Piece::After => self.rate_after * self.expectation.expected[last],
    };
    let own = match self.expectation.nodes.get(&key) {
      Some(node) => match node.survival.piece_at(side, edge) {
        Piece::Before => inverse_tc,
        Piece::Inside(i) => {
          let cell = match side {
            Side::Left => grid.partition_point(|&g| g < edge),
            Side::Right => grid.partition_point(|&g| g <= edge),
          };
          self.rates().rate(cell) * node.survival.value_in(i, edge)
        },
        Piece::After => 0.0,
      },
      None => 0.0,
    };
    pairs * inverse_tc + weight * expected - weight * own_weight * own
  }

  fn inverse_tc_integral_at(&self, x: f64) -> f64 {
    let grid = &self.expectation.grid;
    let n = grid.len();
    if x <= grid[0] {
      return (x - grid[0]) * self.rate_before;
    }
    if x >= grid[n - 1] {
      return self.inverse_tc_integral[n - 1] + (x - grid[n - 1]) * self.rate_after;
    }
    interpolate(grid, &self.inverse_tc_integral, x)
  }

  fn expected_integral_at(&self, x: f64) -> f64 {
    let grid = &self.expectation.grid;
    let n = grid.len();
    if x <= grid[0] {
      return (x - grid[0]) * self.rate_before * self.expectation.expected_before;
    }
    if x >= grid[n - 1] {
      return self.expected_integral[n - 1] + (x - grid[n - 1]) * self.rate_after * self.expectation.expected[n - 1];
    }
    interpolate(grid, &self.expected_integral, x)
  }

  pub(crate) fn own_integral_at(&self, key: GraphNodeKey, x: f64) -> f64 {
    let (Some(node), Some(integral)) = (self.expectation.nodes.get(&key), self.own_integral.get(&key)) else {
      return 0.0;
    };
    let survival = &node.survival;
    let first = survival.first();
    if x < first {
      return self.inverse_tc_integral_at(x) - self.inverse_tc_integral_at(first);
    }
    if x >= survival.last() {
      return integral.cumulative[survival.t.len() - 1];
    }
    let i = survival.t.partition_point(|&tk| tk <= x) - 1;
    integral.cumulative[i] + survival_integral_in(survival, i, &self.rates(), x)
  }

  fn rates(&self) -> CellRates<'_> {
    CellRates {
      grid: &self.expectation.grid,
      inside: &self.cell_rate,
      before: self.rate_before,
      after: self.rate_after,
    }
  }
}

#[derive(Clone, Debug)]
struct OwnIntegral {
  cumulative: Vec<f64>,
}

impl OwnIntegral {
  fn new(survival: &NodeSurvival, rates: &CellRates<'_>) -> Self {
    let mut cumulative = Vec::with_capacity(survival.t.len());
    cumulative.push(0.0);
    for i in 0..survival.t.len().saturating_sub(1) {
      cumulative.push(cumulative[i] + survival_integral_in(survival, i, rates, survival.t[i + 1]));
    }
    Self { cumulative }
  }
}

fn survival_integral_in(survival: &NodeSurvival, segment: usize, rates: &CellRates<'_>, x: f64) -> f64 {
  let grid = rates.grid;
  let start = survival.t[segment];
  let mut cell = grid.partition_point(|&g| g <= start);
  let mut from = start;
  let mut total = 0.0;
  for &g in grid[cell..].iter().take_while(|&&g| g < x) {
    total +=
      rates.rate(cell) * (g - from) * f64::midpoint(survival.value_in(segment, from), survival.value_in(segment, g));
    from = g;
    cell += 1;
  }
  total + rates.rate(cell) * (x - from) * f64::midpoint(survival.value_in(segment, from), survival.value_in(segment, x))
}

struct CellRates<'a> {
  grid: &'a [f64],
  inside: &'a [f64],
  before: f64,
  after: f64,
}

impl CellRates<'_> {
  fn rate(&self, cell: usize) -> f64 {
    match cell {
      0 => self.before,
      cell if cell == self.grid.len() => self.after,
      cell => self.inside[cell - 1],
    }
  }
}

#[derive(Clone, Debug)]
pub(crate) struct LineageExpectation {
  grid: Vec<f64>,
  expected: Vec<f64>,
  expected_pairs: Vec<f64>,
  variance: Vec<f64>,
  expected_before: f64,
  step: f64,
  nodes: BTreeMap<GraphNodeKey, WeightedSurvival>,
}

impl LineageExpectation {
  pub(crate) fn of_inference(graph: &Graph, time_inference: &TimeInference) -> Result<Arc<Self>, Report> {
    if !time_inference.posterior_summaries.lineage_survival {
      return make_internal_error!(
        "The time inference holds no survival functions, so it cannot give expected lineage counts"
      );
    }
    Self::from_posterior(graph, &time_inference.posterior, &time_inference.bad_branches)
      .wrap_err("Failed to compute the expected coalescent lineage counts")
  }

  pub(crate) fn from_posterior(
    graph: &Graph,
    posterior: &BTreeMap<GraphNodeKey, NodePosterior>,
    bad_branches: &BTreeMap<GraphNodeKey, bool>,
  ) -> Result<Arc<Self>, Report> {
    Self::from_survivals(graph, |key| {
      if bad_branches.get(&key).copied().unwrap_or(false) {
        return Ok(None);
      }
      let entry = posterior.get(&key);
      Ok(
        match (entry.and_then(|p| p.survival.as_ref()), entry.and_then(|p| p.time)) {
          (Some(survival), _) => Some(Arc::clone(survival)),
          (None, Some(time)) => Some(Arc::new(NodeSurvival::point(time))),
          (None, None) => None,
        },
      )
    })
  }

  pub(crate) fn from_gaussian(
    graph: &Graph,
    marginals: &BTreeMap<GraphNodeKey, NodeTimeGaussian>,
    tip_times: &NodeTimes,
    excluded: &BTreeMap<GraphNodeKey, bool>,
  ) -> Result<Arc<Self>, Report> {
    let tip_time = |key: GraphNodeKey| {
      tip_times
        .get(&key)
        .copied()
        .flatten()
        .filter(|_| !excluded.get(&key).copied().unwrap_or(false))
    };
    let schedule = GraphSchedule::new(graph)?;
    let mut earliest_tip = BTreeMap::new();
    for &key in schedule.postorder() {
      let children = schedule.children(key);
      let earliest = if children.is_empty() {
        tip_time(key).unwrap_or(f64::INFINITY)
      } else {
        children
          .iter()
          .map(|(child, _)| earliest_tip[child])
          .fold(f64::INFINITY, f64::min)
      };
      earliest_tip.insert(key, earliest);
    }
    Self::from_survivals(graph, |key| {
      let is_leaf = graph.get_node(key).is_some_and(Node::is_leaf);
      let survival = if is_leaf {
        tip_time(key).map(NodeSurvival::point)
      } else {
        marginals
          .get(&key)
          .filter(|_| earliest_tip[&key].is_finite())
          .map(|marginal| NodeSurvival::truncated_normal(marginal.mean, marginal.std, earliest_tip[&key]))
      };
      Ok(survival.map(Arc::new))
    })
  }

  fn from_survivals(
    graph: &Graph,
    survival: impl Fn(GraphNodeKey) -> Result<Option<Arc<NodeSurvival>>, Report>,
  ) -> Result<Arc<Self>, Report> {
    let mut always_surviving = 0.0;
    let mut nodes = BTreeMap::new();
    for node in graph.get_nodes() {
      let weight = 1.0 - node.outbound().len().to_f64_exact();
      match survival(node.key())? {
        Some(survival) => {
          nodes.insert(node.key(), WeightedSurvival { weight, survival });
        },
        None => always_surviving += weight,
      }
    }
    Self::aggregate(nodes, always_surviving)
  }

  pub(crate) fn with_topology(
    &self,
    graph: &Graph,
    merger_times: &BTreeMap<GraphNodeKey, f64>,
  ) -> Result<Arc<Self>, Report> {
    Self::from_survivals(graph, |key| {
      Ok(match (self.nodes.get(&key), merger_times.get(&key)) {
        (Some(node), _) => Some(Arc::clone(&node.survival)),
        (None, Some(&time)) => Some(Arc::new(NodeSurvival::point(time))),
        (None, None) => None,
      })
    })
  }

  fn aggregate(nodes: BTreeMap<GraphNodeKey, WeightedSurvival>, always_surviving: f64) -> Result<Arc<Self>, Report> {
    let weighted: Vec<(f64, &NodeSurvival)> = nodes
      .values()
      .map(|node| (node.weight, node.survival.as_ref()))
      .collect();
    if weighted.is_empty() {
      return make_error!("Cannot build the expected lineage count: no node has a time");
    }

    let t_min = weighted.iter().map(|(_, s)| s.first()).fold(f64::INFINITY, f64::min);
    let t_max = weighted.iter().map(|(_, s)| s.last()).fold(f64::NEG_INFINITY, f64::max);
    let span = (t_max - t_min).max(f64::EPSILON);
    let step = span / (EXPECTATION_GRID_POINTS - 1).to_f64_exact();
    let step_times: Vec<f64> = weighted
      .iter()
      .filter(|(_, survival)| survival.is_point())
      .map(|(_, survival)| survival.first())
      .sorted_by_key(|&t| OrderedFloat(t))
      .dedup()
      .collect();
    let grid: Vec<f64> = (0..EXPECTATION_GRID_POINTS)
      .map(|i| t_min + step * i.to_f64_exact())
      .filter(|&g| {
        step_times
          .binary_search_by_key(&OrderedFloat(g), |&t| OrderedFloat(t))
          .is_err()
      })
      .chain(step_times.iter().flat_map(|&t| [t, t]))
      .sorted_by_key(|&t| OrderedFloat(t))
      .collect();

    let mut expected = vec![always_surviving; grid.len()];
    let mut variance = vec![0.0; grid.len()];
    let mut own_squares = vec![0.0; grid.len()];
    let mut surviving_from = vec![0.0; grid.len() + 1];
    let mut squares_from = vec![0.0; grid.len() + 1];
    for &(weight, survival) in &weighted {
      let first = grid.partition_point(|&g| g < survival.first());
      if survival.is_point() {
        surviving_from[first + 1] += weight;
        squares_from[first + 1] += weight * weight;
        continue;
      }
      let last = grid.partition_point(|&g| g < survival.last());
      surviving_from[first] += weight;
      squares_from[first] += weight * weight;
      let values = survival.values_at_sorted(&grid[first..last]);
      for (s, expected, variance, own_square) in izip!(
        values,
        &mut expected[first..last],
        &mut variance[first..last],
        &mut own_squares[first..last]
      ) {
        *expected += weight * s;
        *variance += weight * weight * s * (1.0 - s);
        *own_square += (weight * s).powi(2);
      }
    }
    let (mut running, mut running_squares) = (0.0, 0.0);
    for i in (0..grid.len()).rev() {
      running += surviving_from[i + 1];
      running_squares += squares_from[i + 1];
      expected[i] += running;
      own_squares[i] += running_squares;
    }
    let expected_pairs = izip!(&expected, &own_squares)
      .map(|(k, own_square)| 0.5 * (k * k - own_square))
      .collect();
    let expected_before = always_surviving + weighted.iter().map(|(w, _)| w).sum::<f64>();
    drop(weighted);

    Ok(Arc::new(Self {
      grid,
      expected,
      expected_pairs,
      variance,
      expected_before,
      step,
      nodes,
    }))
  }

  pub(crate) fn span(&self) -> TimeSpan {
    let t_min = self
      .expected
      .iter()
      .position(|&k| k - self.expected_before >= 0.5)
      .map_or(self.grid[0], |i| self.grid[i]);
    TimeSpan {
      t_min,
      t_max: self.grid[self.grid.len() - 1],
    }
  }

  pub(crate) fn merger_rate_schedule(&self, tc: &Distribution) -> Result<PiecewiseConstantFn, Report> {
    let rate = |k: f64, t: f64| -> Result<f64, Report> {
      let tc = tc
        .eval(t)
        .wrap_err_with(|| format!("When evaluating coalescent Tc at calendar time {t:.6e}"))?;
      if !tc.is_finite() || tc <= 0.0 {
        return make_error!("Coalescent Tc must be finite and positive at calendar time {t:.6e}, got {tc:.6e}");
      }
      Ok(0.5 * f64::max(0.5, k - 1.0) / tc)
    };
    let n = self.grid.len();
    let mut edges: Vec<usize> = (0..n).step_by(MERGER_RATE_STRIDE).chain(once(n - 1)).collect();
    edges.dedup_by_key(|&mut i| OrderedFloat(self.grid[i]));
    let mut rates = vec![rate(self.expected_before, self.grid[0])?];
    for (&start, &end) in edges.iter().zip(&edges[1..]) {
      let k = self.expected[usize::midpoint(start, end)];
      rates.push(rate(k, f64::midpoint(self.grid[start], self.grid[end]))?);
    }
    let breakpoints: Vec<f64> = edges[..edges.len() - 1].iter().map(|&i| self.grid[i]).collect();
    Ok(PiecewiseConstantFn::new(
      Array1::from_vec(breakpoints),
      Array1::from_vec(rates),
    ))
  }

  pub(crate) fn segment_statistics(&self, boundaries: &[f64]) -> SegmentStatistics {
    let n_seg = boundaries.len() - 1;
    let segment_of = |t: f64| boundaries.partition_point(|&b| b <= t).saturating_sub(1).min(n_seg - 1);
    let pair_rate = |i: usize| {
      let k = self.expected[i];
      0.5 * (k * k - k + self.variance[i]).max(0.0)
    };
    let mut pair_integral = vec![0.0; n_seg];
    for i in 1..self.grid.len() {
      let (t0, t1) = (self.grid[i - 1], self.grid[i]);
      if t1 <= t0 {
        continue;
      }
      let (r0, r1) = (pair_rate(i - 1), pair_rate(i));
      let rate_at = |t: f64| r0 + (t - t0) / (t1 - t0) * (r1 - r0);
      let cuts = once(t0)
        .chain(boundaries.iter().copied().filter(|&b| b > t0 && b < t1))
        .chain(once(t1));
      for (a, b) in cuts.tuple_windows() {
        pair_integral[segment_of(f64::midpoint(a, b))] += (b - a) * f64::midpoint(rate_at(a), rate_at(b));
      }
    }
    let mut mergers = vec![0.0; n_seg];
    for node in self.nodes.values().filter(|node| node.weight < 0.0) {
      for (i, merger) in mergers.iter_mut().enumerate() {
        let before = if i == 0 { 1.0 } else { node.survival.at(boundaries[i]) };
        let after = if i == n_seg - 1 {
          0.0
        } else {
          node.survival.at(boundaries[i + 1])
        };
        *merger += -node.weight * (before - after);
      }
    }
    SegmentStatistics { pair_integral, mergers }
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct TimeSpan {
  pub t_min: f64,
  pub t_max: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct SegmentStatistics {
  pub pair_integral: Vec<f64>,
  pub mergers: Vec<f64>,
}

#[derive(Clone, Debug)]
struct WeightedSurvival {
  weight: f64,
  survival: Arc<NodeSurvival>,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct NodeSurvival {
  t: Vec<f64>,
  survival: Vec<f64>,
}

impl NodeSurvival {
  pub(crate) fn of_posterior(distribution: &Distribution<NegLog>) -> Result<Option<Self>, Report> {
    let survival = posterior_survival(distribution, &survival_levels())?;
    Ok(survival.map(Self::from))
  }

  fn point(time: f64) -> Self {
    Self {
      t: vec![time],
      survival: vec![0.0],
    }
  }

  fn truncated_normal(mean: f64, std: f64, upper: f64) -> Self {
    let std = std.max(f64::EPSILON * mean.abs().max(1.0));
    let z_upper = (upper - mean) / std;
    let normal = Normal::standard();
    let cdf_upper = normal.cdf(z_upper);
    let tail_rate = (z_upper <= -TRUNCATED_NORMAL_SPAN).then(|| (mean - upper) / std.powi(2));
    let time_below = |mass: f64| match tail_rate {
      Some(rate) => upper + mass.ln() / rate,
      None => mean + std * normal.inverse_cdf(mass * cdf_upper),
    };
    let survival_at = |t: f64| match tail_rate {
      Some(rate) => -(-rate * (upper - t)).exp_m1(),
      None => ((cdf_upper - normal.cdf((t - mean) / std)) / cdf_upper).clamp(0.0, 1.0),
    };
    let mut survival = PosteriorSurvival {
      t: vec![],
      survival: vec![],
    };
    for level in survival_levels() {
      let mass_below = match level {
        SurvivalLevel::Below(mass) => mass,
        SurvivalLevel::Above(mass) => 1.0 - mass,
      };
      let t = time_below(mass_below).min(upper);
      survival.push(t, survival_at(t));
    }
    Self::from(survival.with_unit_ends())
  }

  fn is_point(&self) -> bool {
    self.t.len() == 1
  }

  fn first(&self) -> f64 {
    self.t[0]
  }

  fn last(&self) -> f64 {
    self.t[self.t.len() - 1]
  }

  pub(crate) fn at(&self, x: f64) -> f64 {
    if x < self.first() {
      return 1.0;
    }
    if x >= self.last() {
      return 0.0;
    }
    self.value_in(self.t.partition_point(|&ti| ti <= x) - 1, x)
  }

  fn values_at_sorted<'a>(&'a self, xs: &'a [f64]) -> impl Iterator<Item = f64> + 'a {
    let (first, last) = (self.first(), self.last());
    xs.iter().scan(0, move |segment, &x| {
      if x < first {
        return Some(1.0);
      }
      if x >= last {
        return Some(0.0);
      }
      *segment += self.t[*segment + 1..].partition_point(|&ti| ti <= x);
      Some(self.value_in(*segment, x))
    })
  }

  fn piece_at(&self, side: Side, x: f64) -> Piece<usize> {
    let i = match side {
      Side::Left => self.t.partition_point(|&ti| ti < x),
      Side::Right => self.t.partition_point(|&ti| ti <= x),
    };
    piece_of(i, self.t.len())
  }

  fn value_in(&self, segment: usize, x: f64) -> f64 {
    let (t0, t1) = (self.t[segment], self.t[segment + 1]);
    let (s0, s1) = (self.survival[segment], self.survival[segment + 1]);
    s0 + (x - t0) / (t1 - t0) * (s1 - s0)
  }
}

impl From<PosteriorSurvival> for NodeSurvival {
  fn from(survival: PosteriorSurvival) -> Self {
    Self {
      t: survival.t,
      survival: survival.survival,
    }
  }
}

pub(crate) fn survival_levels() -> Vec<SurvivalLevel> {
  let interior = (1..INTERIOR_LEVELS).map(|i| i.to_f64_exact() / INTERIOR_LEVELS.to_f64_exact());
  let below = SURVIVAL_TAIL_LEVELS
    .iter()
    .copied()
    .chain(interior.filter(|&mass| mass <= 0.5 && mass > SURVIVAL_TAIL_LEVELS[SURVIVAL_TAIL_LEVELS.len() - 1]))
    .map(SurvivalLevel::Below);
  let interior = (1..INTERIOR_LEVELS).map(|i| i.to_f64_exact() / INTERIOR_LEVELS.to_f64_exact());
  let above = interior
    .filter(|&mass| mass < 0.5 && mass > SURVIVAL_TAIL_LEVELS[SURVIVAL_TAIL_LEVELS.len() - 1])
    .rev()
    .chain(SURVIVAL_TAIL_LEVELS.iter().copied().rev())
    .map(SurvivalLevel::Above);
  below.chain(above).collect()
}

fn piece_slope(t: &[f64], values: &[f64], side: Side, x: f64) -> Piece<f64> {
  match piece_at(t, side, x) {
    Piece::Before => Piece::Before,
    Piece::Inside(i) => Piece::Inside((values[i + 1] - values[i]) / (t[i + 1] - t[i])),
    Piece::After => Piece::After,
  }
}

fn piece_at(t: &[f64], side: Side, x: f64) -> Piece<usize> {
  let i = match side {
    Side::Left => t.partition_point(|&ti| ti < x),
    Side::Right => t.partition_point(|&ti| ti <= x),
  };
  piece_of(i, t.len())
}

fn piece_of(i: usize, len: usize) -> Piece<usize> {
  if i == 0 {
    Piece::Before
  } else if i == len {
    Piece::After
  } else {
    Piece::Inside(i - 1)
  }
}

enum Piece<T> {
  Before,
  Inside(T),
  After,
}

fn interpolate(t: &[f64], values: &[f64], x: f64) -> f64 {
  let i = t.partition_point(|&ti| ti <= x).clamp(1, t.len() - 1) - 1;
  let (t0, t1) = (t[i], t[i + 1]);
  values[i] + (x - t0) / (t1 - t0) * (values[i + 1] - values[i])
}
