use eyre::Report;
use itertools::Itertools;
use ndarray::{Array1, ArrayView1, s};
use ordered_float::OrderedFloat;
use std::iter::once;
use treetime_distribution::{Distribution, NegLog};
use treetime_utils::array::ndarray::partition_point;
use treetime_utils::{make_error, make_internal_error};

pub fn max_posterior_region(
  distribution: &Distribution<NegLog>,
  fraction: f64,
) -> Result<Option<MaxPosteriorRegion>, Report> {
  let normalized = distribution.normalize()?;
  let region = match &normalized {
    Distribution::Empty => return Ok(None),
    Distribution::Point(_) => return Ok(Some(MaxPosteriorRegion::Point)),
    Distribution::Range(_) => grid_max_posterior_region(normalized.t().view(), normalized.y()?.view(), fraction)?,
    Distribution::Function(function) => {
      grid_max_posterior_region(normalized.t().view(), function.y().view(), fraction)?
    },
    Distribution::Formula(_) => {
      return make_internal_error!("Normalization left a formula distribution undiscretized");
    },
  };
  Ok(Some(MaxPosteriorRegion::Grid(region)))
}

pub(crate) fn grid_max_posterior_region(
  t: ArrayView1<f64>,
  neglog: ArrayView1<f64>,
  fraction: f64,
) -> Result<PosteriorRegion, Report> {
  let n_points = t.len();
  if n_points != neglog.len() {
    return make_error!(
      "Time grid and posterior ordinates differ in length: {n_points} and {}",
      neglog.len()
    );
  }
  if n_points < 2 {
    return make_error!("Maximum posterior region needs at least two grid points, got {n_points}");
  }
  if !(fraction > 0.0 && fraction < 1.0) {
    return make_error!("Posterior mass fraction must be in (0, 1), got {fraction}");
  }
  if t.iter().tuple_windows().any(|(a, b)| a >= b) {
    return make_error!("Time grid of a posterior must be strictly increasing");
  }
  if neglog.iter().any(|y| y.is_nan()) {
    return make_error!(
      "Posterior on [{}, {}] has NaN negative log-likelihood values",
      t[0],
      t[n_points - 1]
    );
  }

  let peak = latest_peak_index(neglog);
  let y_peak = neglog[peak];
  if !y_peak.is_finite() {
    return make_error!("Posterior on [{}, {}] has no probability mass", t[0], t[n_points - 1]);
  }

  let cdf = GridCdf::new(t, neglog.iter().map(|y| y - y_peak));

  let region = if peak == n_points - 1 {
    PosteriorRegion {
      lower: cdf.quantile(1.0 - fraction),
      upper: t[n_points - 1],
      peak: PeakLocation::Latest,
    }
  } else if peak == 0 {
    PosteriorRegion {
      lower: t[0],
      upper: cdf.quantile(fraction),
      peak: PeakLocation::Earliest,
    }
  } else {
    let (lower, upper) = highest_density_interval(t, neglog, y_peak, peak, &cdf, fraction);
    PosteriorRegion {
      lower,
      upper,
      peak: PeakLocation::Interior,
    }
  };
  Ok(region)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MaxPosteriorRegion {
  Point,
  Grid(PosteriorRegion),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PosteriorRegion {
  pub lower: f64,
  pub upper: f64,
  pub peak: PeakLocation,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeakLocation {
  Earliest,
  Interior,
  Latest,
}

fn latest_peak_index(neglog: ArrayView1<f64>) -> usize {
  neglog
    .iter()
    .enumerate()
    .rev()
    .min_by_key(|(_, y)| OrderedFloat(**y))
    .map_or(0, |(i, _)| i)
}

fn highest_density_interval(
  t: ArrayView1<f64>,
  neglog: ArrayView1<f64>,
  y_peak: f64,
  peak: usize,
  cdf: &GridCdf<'_>,
  fraction: f64,
) -> (f64, f64) {
  let height = |y: &f64| y - y_peak;
  let earlier = LevelSide::new(
    t.slice_move(s![..=peak;-1]),
    neglog.slice(s![..=peak;-1]).iter().map(height),
  );
  let later = LevelSide::new(t.slice_move(s![peak..]), neglog.slice(s![peak..]).iter().map(height));

  let mass_at = |level: f64, cells: (usize, usize)| {
    cdf.at(later.position(cells.1, level)) - cdf.at(earlier.position(cells.0, level))
  };
  let cells_at = |level: f64| (earlier.cell(level), later.cell(level));

  let levels = earlier
    .envelope
    .iter()
    .merge(&later.envelope)
    .copied()
    .filter(|level| level.is_finite())
    .dedup();

  let reached = levels
    .tuple_windows()
    .find(|&(_, above)| mass_at(above, cells_at(above)) >= fraction);

  let Some((below, above)) = reached else {
    // TODO: this region can hold less than `fraction`; kb/issues/M-timetree-max-posterior-region-fallback-holds-less-than-fraction.md
    return (earlier.support_end(), later.support_end());
  };
  let cells = cells_at(above);
  let mass_below = mass_at(below, cells);
  let mass_above = mass_at(above, cells);
  let level = if mass_below >= fraction {
    below
  } else {
    below + (fraction - mass_below) / (mass_above - mass_below) * (above - below)
  };

  (earlier.position(cells.0, level), later.position(cells.1, level))
}

struct GridCdf<'a> {
  t: ArrayView1<'a, f64>,
  cumulative: Array1<f64>,
}

impl<'a> GridCdf<'a> {
  fn new(t: ArrayView1<'a, f64>, heights: impl Iterator<Item = f64>) -> Self {
    let cells = t
      .iter()
      .zip(heights.map(|h| (-h).exp()))
      .tuple_windows()
      .map(|((t0, p0), (t1, p1))| (t1 - t0) * f64::midpoint(p0, p1));
    let mut cumulative: Array1<f64> = once(0.0)
      .chain(cells.scan(0.0, |acc, mass| {
        *acc += mass;
        Some(*acc)
      }))
      .collect();
    let total = cumulative[cumulative.len() - 1];
    cumulative.mapv_inplace(|c| c / total);
    Self { t, cumulative }
  }

  fn at(&self, x: f64) -> f64 {
    let n_points = self.t.len();
    let i = partition_point(&self.t, |&ti| ti <= x).clamp(1, n_points - 1) - 1;
    let (t0, t1) = (self.t[i], self.t[i + 1]);
    let (c0, c1) = (self.cumulative[i], self.cumulative[i + 1]);
    c0 + (x - t0) / (t1 - t0) * (c1 - c0)
  }

  fn quantile(&self, q: f64) -> f64 {
    let i = self
      .cumulative
      .iter()
      .position(|&c| c >= q)
      .unwrap_or(self.cumulative.len() - 1);
    if i == 0 {
      return self.t[0];
    }
    let (t0, t1) = (self.t[i - 1], self.t[i]);
    let (c0, c1) = (self.cumulative[i - 1], self.cumulative[i]);
    t0 + (q - c0) / (c1 - c0) * (t1 - t0)
  }
}

struct LevelSide<'a> {
  t: ArrayView1<'a, f64>,
  envelope: Array1<f64>,
}

impl<'a> LevelSide<'a> {
  fn new(t: ArrayView1<'a, f64>, heights: impl Iterator<Item = f64>) -> Self {
    let envelope = heights
      .scan(f64::NEG_INFINITY, |running_max, h| {
        *running_max = running_max.max(h);
        Some(*running_max)
      })
      .collect();
    Self { t, envelope }
  }

  fn support_end(&self) -> f64 {
    let finite = partition_point(&self.envelope, |h| h.is_finite());
    self.t[finite - 1]
  }

  fn cell(&self, level: f64) -> usize {
    partition_point(&self.envelope, |&h| h < level)
  }

  fn position(&self, cell: usize, level: f64) -> f64 {
    let n_points = self.t.len();
    if cell == 0 {
      return self.t[0];
    }
    if cell >= n_points {
      return self.t[n_points - 1];
    }
    let (h0, h1) = (self.envelope[cell - 1], self.envelope[cell]);
    let (x0, x1) = (self.t[cell - 1], self.t[cell]);
    if h1.is_infinite() {
      return x0;
    }
    x0 + (level - h0) / (h1 - h0) * (x1 - x0)
  }
}
