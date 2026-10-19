use crate::clock::gaussian::NodeTimeGaussian;
use crate::timetree::inference::result::{NodePosterior, NodeTimes, TimeInference};
use crate::timetree::max_posterior_region::posterior_survival;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use statrs::distribution::{ContinuousCDF, Normal};
use std::collections::BTreeMap;
use std::iter::once;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNodeKey, Node};
use treetime_grid::Side;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_internal_error};

const EXPECTATION_GRID_POINTS: usize = 4001;
const SURVIVAL_POINTS: usize = 256;
const SURVIVAL_TRIM: f64 = 1e-9;
const MERGER_RATE_STRIDE: usize = 10;
const TRUNCATED_NORMAL_SPAN: f64 = 8.0;
const TRUNCATED_NORMAL_POINTS: usize = 1024;
const CURVATURE_GRID_REFINEMENT: usize = 16;
const CURVATURE_GRID_UNIFORM_SHARE: f64 = 0.1;

#[derive(Clone, Debug)]
pub(crate) struct CoalescentPrior {
  expectation: Arc<LineageExpectation>,
  tc: Distribution,
  inverse_tc_integral: Vec<f64>,
  expected_integral: Vec<f64>,
  own_integral: BTreeMap<GraphNodeKey, Vec<f64>>,
  rate_before: f64,
  rate_after: f64,
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
    let mut inverse_tc_integral = vec![0.0; grid.len()];
    let mut expected_integral = vec![0.0; grid.len()];
    for i in 1..grid.len() {
      let dt = grid[i] - grid[i - 1];
      let rate = inverse_tc(f64::midpoint(grid[i - 1], grid[i]))?;
      inverse_tc_integral[i] = inverse_tc_integral[i - 1] + dt * rate;
      expected_integral[i] =
        expected_integral[i - 1] + dt * rate * f64::midpoint(expectation.expected[i - 1], expectation.expected[i]);
    }
    let own_integral = expectation
      .nodes
      .iter()
      .map(|(key, node)| {
        let survival = &node.survival;
        let mut integral = vec![0.0; survival.t.len()];
        for i in 1..survival.t.len() {
          let rate = inverse_tc(f64::midpoint(survival.t[i - 1], survival.t[i]))?;
          integral[i] = integral[i - 1]
            + (survival.t[i] - survival.t[i - 1])
              * rate
              * f64::midpoint(survival.survival[i - 1], survival.survival[i]);
        }
        Ok((*key, integral))
      })
      .collect::<Result<_, Report>>()?;
    let rate_before = inverse_tc(grid[0])?;
    let rate_after = inverse_tc(grid[grid.len() - 1])?;
    Ok(Self {
      expectation: Arc::clone(expectation),
      tc: tc.clone(),
      inverse_tc_integral,
      expected_integral,
      own_integral,
      rate_before,
      rate_after,
    })
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
      PieceSlope::Before => self.rate_before,
      PieceSlope::Inside(slope) => slope,
      PieceSlope::After => self.rate_after,
    };
    let expected = match piece_slope(grid, &self.expected_integral, side, edge) {
      PieceSlope::Before => self.rate_before * self.expectation.expected_before,
      PieceSlope::Inside(slope) => slope,
      PieceSlope::After => self.rate_after * self.expectation.expected[last],
    };
    let own = match (self.expectation.nodes.get(&key), self.own_integral.get(&key)) {
      (Some(node), Some(integral)) => match piece_slope(&node.survival.t, integral, side, edge) {
        PieceSlope::Before => inverse_tc,
        PieceSlope::Inside(slope) => slope,
        PieceSlope::After => 0.0,
      },
      _ => 0.0,
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

  fn own_integral_at(&self, key: GraphNodeKey, x: f64) -> f64 {
    let (Some(node), Some(integral)) = (self.expectation.nodes.get(&key), self.own_integral.get(&key)) else {
      return 0.0;
    };
    let t = &node.survival.t;
    let n = t.len();
    if x < t[0] {
      return self.inverse_tc_integral_at(x) - self.inverse_tc_integral_at(t[0]);
    }
    if x >= t[n - 1] {
      return integral[n - 1];
    }
    interpolate(t, integral, x)
  }
}

#[derive(Clone, Debug)]
pub(crate) struct LineageExpectation {
  grid: Vec<f64>,
  expected: Vec<f64>,
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
    let mut earliest_tip = BTreeMap::new();
    graph.iter_depth_first_postorder_forward(|node| {
      let earliest = if node.is_leaf {
        tip_time(node.key).unwrap_or(f64::INFINITY)
      } else {
        node
          .child_keys
          .iter()
          .map(|(child, _)| earliest_tip[child])
          .fold(f64::INFINITY, f64::min)
      };
      earliest_tip.insert(node.key, earliest);
      Ok(())
    })?;
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

    let t_min = weighted.iter().map(|(_, s)| s.t[0]).fold(f64::INFINITY, f64::min);
    let t_max = weighted
      .iter()
      .map(|(_, s)| s.t[s.t.len() - 1])
      .fold(f64::NEG_INFINITY, f64::max);
    let span = (t_max - t_min).max(f64::EPSILON);
    let step = span / (EXPECTATION_GRID_POINTS - 1).to_f64_exact();
    let step_times: Vec<f64> = weighted
      .iter()
      .filter(|(_, survival)| survival.is_point())
      .map(|(_, survival)| survival.t[0])
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
    let mut surviving_from = vec![0.0; grid.len() + 1];
    for &(weight, survival) in &weighted {
      let first = grid.partition_point(|&g| g < survival.t[0]);
      if survival.is_point() {
        surviving_from[first + 1] += weight;
        continue;
      }
      let last = grid.partition_point(|&g| g < survival.t[survival.t.len() - 1]);
      surviving_from[first] += weight;
      for i in first..last {
        let s = survival.at(grid[i]);
        expected[i] += weight * s;
        variance[i] += weight * weight * s * (1.0 - s);
      }
    }
    let mut running = 0.0;
    for i in (0..grid.len()).rev() {
      running += surviving_from[i + 1];
      expected[i] += running;
    }
    let expected_before = always_surviving + weighted.iter().map(|(w, _)| w).sum::<f64>();
    drop(weighted);

    Ok(Arc::new(Self {
      grid,
      expected,
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
    Ok(posterior_survival(distribution, SURVIVAL_TRIM)?.map(|survival| Self::compact(survival.t, survival.survival)))
  }

  fn compact(t: Vec<f64>, survival: Vec<f64>) -> Self {
    let n = t.len();
    if n <= 2 {
      return Self { t, survival };
    }
    let first = survival
      .iter()
      .position(|&s| s < 1.0 - SURVIVAL_TRIM)
      .unwrap_or(0)
      .saturating_sub(1);
    let last = (survival.iter().rposition(|&s| s > SURVIVAL_TRIM).unwrap_or(n - 1) + 1).min(n - 1);
    let full = Self {
      t: t[first..=last].to_vec(),
      survival: survival[first..=last].to_vec(),
    };
    if full.t.len() <= SURVIVAL_POINTS {
      return full;
    }
    let (start, end) = (full.t[0], full.t[full.t.len() - 1]);
    let step = (end - start) / (SURVIVAL_POINTS - 1).to_f64_exact();
    let t: Vec<f64> = (0..SURVIVAL_POINTS).map(|i| start + step * i.to_f64_exact()).collect();
    let survival = t.iter().map(|&x| full.at(x)).collect();
    Self { t, survival }
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
    let (start, end) = match tail_rate {
      Some(rate) => (upper + SURVIVAL_TRIM.ln() / rate, upper),
      None => (
        mean + std * normal.inverse_cdf(SURVIVAL_TRIM * cdf_upper),
        upper.min(mean + std * normal.inverse_cdf((1.0 - SURVIVAL_TRIM) * cdf_upper)),
      ),
    };
    let survival = |t: f64| match tail_rate {
      Some(rate) => 1.0 - (-rate * (upper - t)).exp(),
      None => ((cdf_upper - normal.cdf((t - mean) / std)) / cdf_upper).clamp(0.0, 1.0),
    };
    let t = curvature_grid(start, end, TRUNCATED_NORMAL_POINTS, survival);
    let mut survival: Vec<f64> = t.iter().map(|&x| survival(x)).collect();
    survival[0] = 1.0;
    survival[TRUNCATED_NORMAL_POINTS - 1] = 0.0;
    Self { t, survival }
  }

  fn is_point(&self) -> bool {
    self.t.len() == 1
  }

  fn at(&self, x: f64) -> f64 {
    let n = self.t.len();
    if x < self.t[0] {
      return 1.0;
    }
    if x >= self.t[n - 1] {
      return 0.0;
    }
    let i = self.t.partition_point(|&ti| ti <= x) - 1;
    let (t0, t1) = (self.t[i], self.t[i + 1]);
    self.survival[i] + (x - t0) / (t1 - t0) * (self.survival[i + 1] - self.survival[i])
  }
}

fn curvature_grid(start: f64, end: f64, n_points: usize, survival: impl Fn(f64) -> f64) -> Vec<f64> {
  let fine = n_points * CURVATURE_GRID_REFINEMENT;
  let spacing = (end - start) / (fine - 1).to_f64_exact();
  let times: Vec<f64> = (0..fine).map(|i| start + spacing * i.to_f64_exact()).collect();
  let values: Vec<f64> = times.iter().map(|&time| survival(time)).collect();
  let curvature: Vec<f64> = (0..fine)
    .map(|i| {
      let center = i.clamp(1, fine - 2);
      let second_difference = values[center - 1] - 2.0 * values[center] + values[center + 1];
      (second_difference / (spacing * spacing)).abs().sqrt()
    })
    .collect();
  let mean_curvature = curvature.iter().sum::<f64>() / fine.to_f64_exact();
  let density: Vec<f64> = curvature
    .iter()
    .map(|&local| local + CURVATURE_GRID_UNIFORM_SHARE * mean_curvature.max(f64::MIN_POSITIVE))
    .collect();
  let cumulative: Vec<f64> = once(0.0)
    .chain(density.iter().tuple_windows().scan(0.0, |acc, (left, right)| {
      *acc += spacing * f64::midpoint(*left, *right);
      Some(*acc)
    }))
    .collect();
  let total = cumulative[fine - 1];
  (0..n_points)
    .map(|index| {
      if index == 0 {
        return start;
      }
      if index == n_points - 1 {
        return end;
      }
      let target = total * index.to_f64_exact() / (n_points - 1).to_f64_exact();
      interpolate(&cumulative, &times, target)
    })
    .collect()
}

fn piece_slope(t: &[f64], values: &[f64], side: Side, x: f64) -> PieceSlope {
  let i = match side {
    Side::Left => t.partition_point(|&ti| ti < x),
    Side::Right => t.partition_point(|&ti| ti <= x),
  };
  if i == 0 {
    PieceSlope::Before
  } else if i == t.len() {
    PieceSlope::After
  } else {
    PieceSlope::Inside((values[i] - values[i - 1]) / (t[i] - t[i - 1]))
  }
}

enum PieceSlope {
  Before,
  Inside(f64),
  After,
}

fn interpolate(t: &[f64], values: &[f64], x: f64) -> f64 {
  let i = t.partition_point(|&ti| ti <= x).clamp(1, t.len() - 1) - 1;
  let (t0, t1) = (t[i], t[i + 1]);
  values[i] + (x - t0) / (t1 - t0) * (values[i + 1] - values[i])
}
