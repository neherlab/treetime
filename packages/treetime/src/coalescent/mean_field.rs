//! Mean-field coalescent prior on node times.
//!
//! With merger indicators `I_u(t) = [t_u > t]` and node weights `w_u = 1 - n_children(u)` (+1 for a tip, -m for a node
//! of m mergers), the lineage count is `k(t) = sum_u w_u I_u(t)`. For a fixed topology the Kingman log-density is
//!
//! ```text
//! log p = -int k(t) (k(t) - 1) / (2 Tc(t)) dt - sum_v m_v ln Tc(t_v)
//! ```
//!
//! which couples all node times through `k(t)`. Linearized in the indicators of the other nodes around their
//! expectations `S_u(t) = P(t_u > t)`, the terms of node v are local:
//!
//! ```text
//! -log f_v(x) = int^x [w_v (w_v - 1) / 2 + w_v kbar_{-v}(t)] / Tc(t) dt + m_v ln Tc(x)
//! ```
//!
//! with `kbar_{-v}(t) = kbar(t) - w_v S_v(t)` and `kbar(t) = sum_u w_u S_u(t)`. For a binary internal node the slope
//! is `(kbar_{-v} - 1) / Tc` and the merger factor `1 / Tc` carries no step, unlike a lineage count frozen at point
//! times. The expectations come from the node posteriors of an earlier pass.

use crate::timetree::inference::result::NodePosterior;
use crate::timetree::max_posterior_region::posterior_survival;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::make_error;

/// Number of points of the global grid that carries `kbar(t)` and its integrals
const EXPECTATION_GRID_POINTS: usize = 4001;
/// Largest number of points kept of one node's survival function
const SURVIVAL_POINTS: usize = 256;
/// Survival below this, or above one minus it, counts as 0 or 1 when a survival function is trimmed
const SURVIVAL_TRIM: f64 = 1e-9;

/// Survival function `S(t) = P(T > t)` of one node on its own grid: 1 before the first time, 0 after the last
#[derive(Clone, Debug)]
struct NodeSurvival {
  t: Vec<f64>,
  survival: Vec<f64>,
  weight: f64,
}

impl NodeSurvival {
  /// Trimmed to where the survival moves away from 1 and 0, and resampled to at most `SURVIVAL_POINTS` uniform times
  fn compact(t: Vec<f64>, survival: Vec<f64>) -> Self {
    let n = t.len();
    if n <= 2 {
      return Self {
        t,
        survival,
        weight: 0.0,
      };
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
      weight: 0.0,
    };
    if full.t.len() <= SURVIVAL_POINTS {
      return full;
    }
    let (start, end) = (full.t[0], full.t[full.t.len() - 1]);
    let step = (end - start) / (SURVIVAL_POINTS - 1).to_f64_exact();
    let t: Vec<f64> = (0..SURVIVAL_POINTS).map(|i| start + step * i.to_f64_exact()).collect();
    let survival = t.iter().map(|&x| full.at(x)).collect();
    Self {
      t,
      survival,
      weight: 0.0,
    }
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

/// Expected lineage count `kbar(t)` and the survival function of every node, from the node posteriors of one pass.
/// Nodes on bad branches are counted as surviving at all times, as `compute_lineage_counts` counts them at the
/// present.
#[derive(Clone, Debug)]
pub(crate) struct LineageExpectation {
  grid: Vec<f64>,
  expected: Vec<f64>,
  /// `V(t) = sum_u w_u^2 S_u(t) (1 - S_u(t))`, the variance of the lineage count for independent node times
  variance: Vec<f64>,
  /// `kbar` before the first grid time
  expected_before: f64,
  nodes: BTreeMap<GraphNodeKey, NodeSurvival>,
}

impl LineageExpectation {
  pub(crate) fn from_posterior(
    graph: &Graph,
    posterior: &BTreeMap<GraphNodeKey, NodePosterior>,
    bad_branches: &BTreeMap<GraphNodeKey, bool>,
  ) -> Result<Arc<Self>, Report> {
    let mut always_surviving = 0.0;
    let mut nodes = BTreeMap::new();
    for node in graph.get_nodes() {
      let key = node.key();
      let weight = 1.0 - node.outbound().len().to_f64_exact();
      if bad_branches.get(&key).copied().unwrap_or(false) {
        always_surviving += weight;
        continue;
      }
      let entry = posterior.get(&key);
      let survival = match entry.and_then(|p| p.distribution.as_deref()) {
        Some(distribution) => posterior_survival(distribution)
          .wrap_err_with(|| format!("When computing the survival function of node {key}"))?,
        None => None,
      };
      let survival = match (survival, entry.and_then(|p| p.time)) {
        (Some((t, survival)), _) => NodeSurvival {
          weight,
          ..NodeSurvival::compact(t, survival)
        },
        (None, Some(time)) => NodeSurvival {
          t: vec![time],
          survival: vec![0.0],
          weight,
        },
        (None, None) => {
          always_surviving += weight;
          continue;
        },
      };
      nodes.insert(key, survival);
    }
    let weighted: Vec<(f64, &NodeSurvival)> = nodes.values().map(|survival| (survival.weight, survival)).collect();
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
    let grid: Vec<f64> = (0..EXPECTATION_GRID_POINTS)
      .map(|i| t_min + step * i.to_f64_exact())
      .collect();

    let mut expected = vec![always_surviving; grid.len()];
    let mut variance = vec![0.0; grid.len()];
    let mut surviving_from = vec![0.0; grid.len() + 1];
    for &(weight, survival) in &weighted {
      let first = grid.partition_point(|&g| g < survival.t[0]);
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
      nodes,
    }))
  }

  /// Expected coalescent statistics of the segments between `boundaries` (the first and last segments extend to the
  /// past and the future): the pair-rate integral `I_i = int_i E[k (k - 1) / 2] dt` with
  /// `E[k (k - 1)] = kbar^2 - kbar + V` for independent node times, and the expected number of mergers
  /// `M_i = sum_v m_v P(t_v in segment i)`. With point times they are the statistics of the lineage counts.
  pub(crate) fn segment_statistics(&self, boundaries: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n_seg = boundaries.len() - 1;
    let segment_of = |t: f64| boundaries.partition_point(|&b| b <= t).saturating_sub(1).min(n_seg - 1);
    let pair_rate = |i: usize| {
      let k = self.expected[i];
      0.5 * (k * k - k + self.variance[i]).max(0.0)
    };
    let mut pair_integral = vec![0.0; n_seg];
    for i in 1..self.grid.len() {
      let mid = f64::midpoint(self.grid[i - 1], self.grid[i]);
      pair_integral[segment_of(mid)] +=
        (self.grid[i] - self.grid[i - 1]) * f64::midpoint(pair_rate(i - 1), pair_rate(i));
    }
    let mut mergers = vec![0.0; n_seg];
    for survival in self.nodes.values().filter(|survival| survival.weight < 0.0) {
      for (i, merger) in mergers.iter_mut().enumerate() {
        let before = if i == 0 { 1.0 } else { survival.at(boundaries[i]) };
        let after = if i == n_seg - 1 {
          0.0
        } else {
          survival.at(boundaries[i + 1])
        };
        *merger += -survival.weight * (before - after);
      }
    }
    (pair_integral, mergers)
  }
}

/// Mean-field coalescent prior for one coalescent time scale; see the module documentation
#[derive(Clone, Debug)]
pub(crate) struct CoalescentPrior {
  expectation: Arc<LineageExpectation>,
  tc: Distribution,
  /// `Q(g) = int_{g_0}^{g} dt / Tc(t)` on the grid of the expectation
  inverse_tc_integral: Vec<f64>,
  /// `B(g) = int_{g_0}^{g} kbar(t) / Tc(t) dt`
  expected_integral: Vec<f64>,
  /// `C_v(t) = int_{t_0}^{t} S_v(s) / Tc(s) ds` on the grid of each node's survival function
  own_integral: BTreeMap<GraphNodeKey, Vec<f64>>,
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
      .map(|(key, survival)| {
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
    Ok(Self {
      expectation: Arc::clone(expectation),
      tc: tc.clone(),
      inverse_tc_integral,
      expected_integral,
      own_integral,
    })
  }

  /// Negative log prior of node `key` with `n_children` children at time `x` (the module's `-log f_v`, up to a
  /// constant per node). Tips (no children) take the same form with weight 1 and no merger factor.
  pub(crate) fn node_contribution(&self, key: GraphNodeKey, n_children: usize, x: f64) -> Result<f64, Report> {
    let weight = 1.0 - n_children.to_f64_exact();
    let n_mergers = n_children.saturating_sub(1).to_f64_exact();
    let own = self.own_integral_at(key, x);
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
        - weight * weight * own
        + merger,
    )
  }

  fn rate_outside(&self, t: f64) -> f64 {
    self.tc.eval(t).map_or(0.0, |tc| if tc > 0.0 { 1.0 / tc } else { 0.0 })
  }

  fn inverse_tc_integral_at(&self, x: f64) -> f64 {
    let grid = &self.expectation.grid;
    let n = grid.len();
    if x <= grid[0] {
      return (x - grid[0]) * self.rate_outside(grid[0]);
    }
    if x >= grid[n - 1] {
      return self.inverse_tc_integral[n - 1] + (x - grid[n - 1]) * self.rate_outside(grid[n - 1]);
    }
    interpolate(grid, &self.inverse_tc_integral, x)
  }

  fn expected_integral_at(&self, x: f64) -> f64 {
    let grid = &self.expectation.grid;
    let n = grid.len();
    if x <= grid[0] {
      return (x - grid[0]) * self.rate_outside(grid[0]) * self.expectation.expected_before;
    }
    if x >= grid[n - 1] {
      return self.expected_integral[n - 1]
        + (x - grid[n - 1]) * self.rate_outside(grid[n - 1]) * self.expectation.expected[n - 1];
    }
    interpolate(grid, &self.expected_integral, x)
  }

  /// `C_v(x)`; 0 for a node without an earlier posterior, which then sees `kbar` unchanged
  fn own_integral_at(&self, key: GraphNodeKey, x: f64) -> f64 {
    let (Some(survival), Some(integral)) = (self.expectation.nodes.get(&key), self.own_integral.get(&key)) else {
      return 0.0;
    };
    let t = &survival.t;
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

fn interpolate(t: &[f64], values: &[f64], x: f64) -> f64 {
  let i = t.partition_point(|&ti| ti <= x).clamp(1, t.len() - 1) - 1;
  let (t0, t1) = (t[i], t[i + 1]);
  values[i] + (x - t0) / (t1 - t0) * (values[i + 1] - values[i])
}
