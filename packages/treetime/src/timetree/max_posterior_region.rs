use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use std::iter::once;
use treetime_distribution::{Distribution, NegLog};
use treetime_utils::{make_error, make_internal_error};

pub fn max_posterior_region(
  distribution: &Distribution<NegLog>,
  fraction: f64,
) -> Result<Option<MaxPosteriorRegion>, Report> {
  let normalized = distribution.normalize()?;
  match &normalized {
    Distribution::Empty => Ok(None),
    Distribution::Point(_) => Ok(Some(MaxPosteriorRegion::Point)),
    Distribution::Range(_) | Distribution::Function(_) => {
      let region = grid_max_posterior_region(&normalized.t(), &normalized.y()?, fraction)?;
      Ok(Some(MaxPosteriorRegion::Grid(region)))
    },
    Distribution::Formula(_) => make_internal_error!("Normalization left a formula distribution undiscretized"),
  }
}

pub(crate) fn grid_max_posterior_region(
  t: &Array1<f64>,
  neglog: &Array1<f64>,
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

  let heights = neglog.mapv(|y| y - y_peak);
  let cdf = GridCdf::new(t, &heights);

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
    let (lower, upper) = highest_density_interval(t, &heights, peak, &cdf, fraction);
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

fn latest_peak_index(neglog: &Array1<f64>) -> usize {
  neglog
    .iter()
    .enumerate()
    .rev()
    .min_by_key(|(_, y)| OrderedFloat(**y))
    .map_or(0, |(i, _)| i)
}

fn highest_density_interval(
  t: &Array1<f64>,
  heights: &Array1<f64>,
  peak: usize,
  cdf: &GridCdf,
  fraction: f64,
) -> (f64, f64) {
  let earlier = LevelSide::new((0..=peak).rev().map(|i| (t[i], heights[i])));
  let later = LevelSide::new((peak..t.len()).map(|i| (t[i], heights[i])));

  let mass_at = |level: f64, cells: (usize, usize)| {
    cdf.at(later.position(cells.1, level)) - cdf.at(earlier.position(cells.0, level))
  };
  let cells_at = |level: f64| (earlier.cell(level), later.cell(level));

  let levels = earlier
    .envelope
    .iter()
    .chain(&later.envelope)
    .copied()
    .filter(|level| level.is_finite())
    .sorted_by_key(|level| OrderedFloat(*level))
    .dedup()
    .collect_vec();

  let reached = levels
    .iter()
    .copied()
    .tuple_windows()
    .find(|&(_, above)| mass_at(above, cells_at(above)) >= fraction);

  let Some((below, above)) = reached else {
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

struct GridCdf {
  t: Vec<f64>,
  cumulative: Vec<f64>,
}

impl GridCdf {
  fn new(t: &Array1<f64>, heights: &Array1<f64>) -> Self {
    let density = heights.mapv(|h| (-h).exp());
    let cells = t
      .iter()
      .zip(&density)
      .tuple_windows()
      .map(|((t0, p0), (t1, p1))| (t1 - t0) * f64::midpoint(*p0, *p1));
    let unnormalized = once(0.0)
      .chain(cells.scan(0.0, |acc, mass| {
        *acc += mass;
        Some(*acc)
      }))
      .collect_vec();
    let total = unnormalized[unnormalized.len() - 1];
    let cumulative = unnormalized.iter().map(|c| c / total).collect_vec();
    Self {
      t: t.to_vec(),
      cumulative,
    }
  }

  fn at(&self, x: f64) -> f64 {
    let n_points = self.t.len();
    let i = self.t.partition_point(|&ti| ti <= x).clamp(1, n_points - 1) - 1;
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

struct LevelSide {
  t: Vec<f64>,
  envelope: Vec<f64>,
}

impl LevelSide {
  fn new(points: impl Iterator<Item = (f64, f64)>) -> Self {
    let (t, heights): (Vec<f64>, Vec<f64>) = points.unzip();
    let envelope = heights
      .iter()
      .scan(f64::NEG_INFINITY, |running_max, &h| {
        *running_max = running_max.max(h);
        Some(*running_max)
      })
      .collect_vec();
    Self { t, envelope }
  }

  fn support_end(&self) -> f64 {
    let finite = self.envelope.partition_point(|h| h.is_finite());
    self.t[finite - 1]
  }

  fn cell(&self, level: f64) -> usize {
    self.envelope.partition_point(|&h| h < level)
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
