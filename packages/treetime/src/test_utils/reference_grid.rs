use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::{Array1, s};
use num_traits::ToPrimitive;
use std::iter::once;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_internal_error};

const GRADING_LEVELS: i32 = 60;
const LOG_CHANGE_PER_PANEL: f64 = 4.0;
const MAX_PANELS: usize = 64;
const SINGULAR_TOLERANCE: f64 = 1e-9;
const EXP_SERIES_RANGE: f64 = 0.5;
const EXP_SERIES_TERMS: u32 = 24;
const GAUSS_LEGENDRE_8: [(f64, f64); 8] = [
  (-0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
  (-0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
  (-0.525_532_409_916_329, 0.313_706_645_877_887_3),
  (-0.183_434_642_495_649_8, 0.362_683_783_378_362),
  (0.183_434_642_495_649_8, 0.362_683_783_378_362),
  (0.525_532_409_916_329, 0.313_706_645_877_887_3),
  (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
  (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
];

#[derive(Clone, Debug)]
pub(crate) struct GridWindow {
  start: isize,
  log_amp: f64,
  values: Array1<f64>,
  young: Option<HardEnd>,
  old: Option<HardEnd>,
  corrections: Option<CellMoments>,
}

impl GridWindow {
  pub(crate) fn from_log_values(
    start: isize,
    log_values: &Array1<f64>,
    ends: Ends<'_>,
    fraction: f64,
  ) -> Result<Self, Report> {
    let peak = log_values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !peak.is_finite() {
      return make_error!(
        "The reference function vanishes on all {} grid points from index {start}",
        log_values.len()
      );
    }
    let values = log_values.mapv(|value| (value - peak).exp());
    let last = start + to_isize(values.len()) - 1;
    let relative = |(at, log_value): (f64, f64), from: isize| HardEnd {
      delta: (at - from.to_f64_exact()).abs(),
      value: (log_value - peak).exp(),
    };
    let window = Self::truncated(
      Self {
        start,
        log_amp: peak,
        values,
        young: ends.young.map(|end| relative(end, start)),
        old: ends.old.map(|end| relative(end, last)),
        corrections: None,
      },
      fraction,
    )?;
    let corrections = ends.exact.map(|exact| CellMoments::of(&window, &exact)).transpose()?;
    Ok(Self { corrections, ..window })
  }

  pub(crate) fn from_values(start: isize, log_amp: f64, values: &Array1<f64>, fraction: f64) -> Result<Self, Report> {
    let values = values.mapv(|value| value.max(0.0));
    let peak = values.iter().copied().fold(0.0, f64::max);
    if !(peak.is_finite() && peak > 0.0) {
      return make_error!(
        "The reference function vanishes on all {} grid points from index {start}",
        values.len()
      );
    }
    Self::truncated(
      Self {
        start,
        log_amp: log_amp + peak.ln(),
        values: values.mapv(|value| value / peak),
        young: None,
        old: None,
        corrections: None,
      },
      fraction,
    )
  }

  pub(crate) fn shifted(&self, offset: isize, log_factor: f64) -> Self {
    Self {
      start: self.start + offset,
      log_amp: self.log_amp + log_factor,
      corrections: self.corrections.as_ref().map(|corrections| corrections.shifted(offset)),
      ..self.clone()
    }
  }

  pub(crate) fn end(&self) -> isize {
    self.start + to_isize(self.values.len()) - 1
  }

  pub(crate) fn log_amp(&self) -> f64 {
    self.log_amp
  }

  pub(crate) fn young_end(&self) -> f64 {
    self.start.to_f64_exact() - self.young.map_or(0.0, |end| end.delta)
  }

  pub(crate) fn old_end(&self) -> f64 {
    self.end().to_f64_exact() + self.old.map_or(0.0, |end| end.delta)
  }

  pub(crate) fn has_young_edge(&self) -> bool {
    self.young.is_some()
  }

  pub(crate) fn has_old_edge(&self) -> bool {
    self.old.is_some()
  }

  pub(crate) fn log_at(&self, x: f64) -> f64 {
    let below = to_isize_floor(x);
    let (lower, upper) = (self.sample(below), self.sample(below + 1));
    let inside = below >= self.start && below < self.end();
    let log_value = if inside && lower > 0.0 && upper > 0.0 {
      lerp(lower.ln(), upper.ln(), x - below.to_f64_exact())
    } else {
      self.at(x).ln()
    };
    log_value + self.log_amp
  }

  pub(crate) fn log_sample(&self, index: isize) -> f64 {
    self.sample(index).ln() + self.log_amp
  }

  pub(crate) fn integral(&self, weight: impl Fn(f64) -> f64) -> f64 {
    let (young, old) = (self.young_end(), self.old_end());
    let product = |x: f64| self.at(x) * weight(x);
    let interior: f64 = once(young)
      .chain(((to_isize_floor(young.floor()) + 1)..to_isize_floor(old.ceil())).map(ToF64Exact::to_f64_exact))
      .chain(once(old))
      .tuple_windows()
      .filter(|(left, right)| right > left)
      .map(|(left, right)| {
        (right - left) / 6.0 * (product(left) + 4.0 * product(f64::midpoint(left, right)) + product(right))
      })
      .sum();
    let corrections = self
      .corrections
      .as_ref()
      .map_or(0.0, |corrections| corrections.integral(&weight));
    interior + corrections
  }

  fn truncated(window: Self, fraction: f64) -> Result<Self, Report> {
    let values = &window.values;
    let (Some(first), Some(last)) = (
      values.iter().position(|&value| value >= fraction),
      values.iter().rposition(|&value| value >= fraction),
    ) else {
      return make_internal_error!("A normalized reference function has no value above {fraction}");
    };
    let first = first.saturating_sub(1);
    let last = (last + 1).min(values.len() - 1);
    Ok(Self {
      start: window.start + to_isize(first),
      values: values.slice(s![first..=last]).to_owned(),
      young: window.young.filter(|_| first == 0),
      old: window.old.filter(|_| last + 1 == values.len()),
      ..window
    })
  }

  fn sample(&self, index: isize) -> f64 {
    usize::try_from(index - self.start)
      .ok()
      .and_then(|offset| self.values.get(offset).copied())
      .unwrap_or(0.0)
  }

  fn at(&self, x: f64) -> f64 {
    if x < self.young_end() || x > self.old_end() {
      return 0.0;
    }
    let (first, last) = (self.start.to_f64_exact(), self.end().to_f64_exact());
    if x < first {
      let end = self.young.unwrap_or(HardEnd { delta: 1.0, value: 0.0 });
      return lerp(end.value, self.values[0], (x - (first - end.delta)) / end.delta);
    }
    if x > last {
      let end = self.old.unwrap_or(HardEnd { delta: 1.0, value: 0.0 });
      return lerp(self.values[self.values.len() - 1], end.value, (x - last) / end.delta);
    }
    let below = x.floor();
    let lower = self.sample(to_isize_floor(below));
    if x > below {
      lerp(lower, self.sample(to_isize_floor(below) + 1), x - below)
    } else {
      lower
    }
  }

  fn stencils(&self) -> Stencils {
    let (young, old) = (self.young_end(), self.old_end());
    let first = to_isize_floor(young.floor());
    let last = to_isize_floor(old.ceil()) + 1;
    let (mean, slope): (Vec<f64>, Vec<f64>) = (first..=last)
      .map(|upper| {
        let upper = upper.to_f64_exact();
        let (from, to) = ((upper - 1.0).max(young), upper.min(old));
        if to < from {
          return (0.0, 0.0);
        }
        let (left, right) = (self.at(from), self.at(to));
        let inside = to - from >= 1.0;
        (
          0.5 * (to - from) * (left + right),
          if inside { right - left } else { 0.0 },
        )
      })
      .unzip();
    Stencils {
      mean: GridSamples {
        start: first,
        values: Array1::from_vec(mean),
      },
      slope: GridSamples {
        start: first,
        values: Array1::from_vec(slope),
      },
    }
  }
}

#[derive(Default)]
pub(crate) struct Ends<'a> {
  pub young: Option<(f64, f64)>,
  pub old: Option<(f64, f64)>,
  pub exact: Option<ExactFunction<'a>>,
}

pub(crate) struct ExactFunction<'a> {
  pub log_value: &'a dyn Fn(f64) -> Result<f64, Report>,
  pub breaks: &'a [f64],
  pub singular: &'a [f64],
  pub log_linear: bool,
}

#[derive(Clone, Copy, Debug)]
struct HardEnd {
  delta: f64,
  value: f64,
}

#[derive(Clone, Debug)]
struct CellMoments {
  start: isize,
  d0: Array1<f64>,
  d1: Array1<f64>,
}

impl CellMoments {
  fn of(window: &GridWindow, exact: &ExactFunction<'_>) -> Result<Self, Report> {
    let (young, old) = (window.young_end(), window.old_end());
    let first = to_isize_floor(young.floor());
    let defect = |x: f64| -> Result<f64, Report> { Ok(((exact.log_value)(x)? - window.log_amp).exp() - window.at(x)) };
    let (d0, d1): (Vec<f64>, Vec<f64>) = (first..to_isize_floor(old.ceil()))
      .map(|cell| {
        let left = cell.to_f64_exact();
        let (from, to) = (left.max(young), (left + 1.0).min(old));
        if to <= from {
          return Ok((0.0, 0.0));
        }
        cell_moments(from, to, left + 0.5, window, exact, &defect)
      })
      .collect::<Result<Vec<_>, Report>>()?
      .into_iter()
      .unzip();
    Ok(Self {
      start: first,
      d0: Array1::from_vec(d0),
      d1: Array1::from_vec(d1),
    })
  }

  fn shifted(&self, offset: isize) -> Self {
    Self {
      start: self.start + offset,
      ..self.clone()
    }
  }

  fn integral(&self, weight: impl Fn(f64) -> f64) -> f64 {
    izip!(self.d0.iter(), self.d1.iter())
      .enumerate()
      .map(|(offset, (&d0, &d1))| {
        let left = (self.start + to_isize(offset)).to_f64_exact();
        let (w0, w1) = (weight(left), weight(left + 1.0));
        d0 * 0.5 * (w0 + w1) + d1 * (w1 - w0)
      })
      .sum()
  }

  fn point_masses(&self) -> GridSamples {
    let mut values = Array1::zeros(self.d0.len() + 1);
    for (offset, (&d0, &d1)) in izip!(self.d0.iter(), self.d1.iter()).enumerate() {
      values[offset] += 0.5 * d0 - d1;
      values[offset + 1] += 0.5 * d0 + d1;
    }
    GridSamples {
      start: self.start,
      values,
    }
  }

  fn as_samples(&self) -> (GridSamples, GridSamples) {
    let samples = |values: &Array1<f64>| GridSamples {
      start: self.start,
      values: values.clone(),
    };
    (samples(&self.d0), samples(&self.d1))
  }
}

struct Stencils {
  mean: GridSamples,
  slope: GridSamples,
}

struct CorrectionTerms {
  terms: Vec<(GridSamples, GridSamples, f64)>,
}

impl CorrectionTerms {
  fn new(subtree: &GridWindow, branch: &GridWindow) -> Self {
    let against = |moments: &CellMoments, other: &GridWindow| {
      let ((d0, d1), stencils) = (moments.as_samples(), other.stencils());
      [(d0, stencils.mean, 1.0), (d1, stencils.slope, -1.0)]
    };
    let cross = branch
      .corrections
      .as_ref()
      .zip(subtree.corrections.as_ref())
      .map(|(left, right)| (left.point_masses(), right.point_masses(), 1.0));
    let terms = branch
      .corrections
      .as_ref()
      .map(|moments| against(moments, subtree))
      .into_iter()
      .chain(subtree.corrections.as_ref().map(|moments| against(moments, branch)))
      .flatten()
      .chain(cross)
      .collect();
    Self { terms }
  }

  fn at(&self, parent: isize) -> f64 {
    self
      .terms
      .iter()
      .map(|(left, right, sign)| sign * left.convolved_at(right, parent))
      .sum()
  }
}

#[derive(Clone, Debug)]
struct GridSamples {
  start: isize,
  values: Array1<f64>,
}

impl GridSamples {
  fn convolved_at(&self, other: &Self, parent: isize) -> f64 {
    let lo = self.start.max(parent - other.end());
    let hi = self.end().min(parent - other.start);
    if hi < lo {
      return 0.0;
    }
    self
      .values
      .slice(s![to_usize(lo - self.start)..=to_usize(hi - self.start)])
      .dot(
        &other
          .values
          .slice(s![to_usize(parent - hi - other.start)..=to_usize(parent - lo - other.start);-1]),
      )
  }

  fn end(&self) -> isize {
    self.start + to_isize(self.values.len()) - 1
  }
}

pub(crate) fn convolve_windows(
  subtree: &GridWindow,
  branch: &GridWindow,
  spacing: f64,
  fraction: f64,
) -> Result<GridWindow, Report> {
  let first = to_isize_floor((subtree.young_end() + branch.young_end()).ceil());
  let last = to_isize_floor((subtree.old_end() + branch.old_end()).floor());
  let reversed = subtree.values.slice(s![..;-1]).to_owned();
  let corrections = CorrectionTerms::new(subtree, branch);
  let values: Array1<f64> = (first..=last)
    .map(|parent| message_value(subtree, branch, &reversed, parent) + corrections.at(parent))
    .collect();
  GridWindow::from_values(
    first,
    subtree.log_amp + branch.log_amp + spacing.ln(),
    &values,
    fraction,
  )
}

fn to_isize(value: usize) -> isize {
  isize::try_from(value).expect("grid index must fit into isize")
}

pub(crate) fn to_isize_floor(value: f64) -> isize {
  value
    .floor()
    .to_isize()
    .expect("grid index must be finite and fit into isize")
}

fn message_value(subtree: &GridWindow, branch: &GridWindow, reversed: &Array1<f64>, parent: isize) -> f64 {
  let parent_x = parent.to_f64_exact();
  let k_lo = branch.young_end().max(parent_x - subtree.old_end());
  let k_hi = branch.old_end().min(parent_x - subtree.young_end());
  if k_hi <= k_lo {
    return 0.0;
  }
  let product = |k: f64| branch.at(k) * subtree.at(parent_x - k);
  let k1 = to_isize_floor(k_lo.ceil());
  let k2 = to_isize_floor(k_hi);
  if k1 > k2 {
    return 0.5 * (k_hi - k_lo) * (product(k_lo) + product(k_hi));
  }
  let sample = |k: isize| branch.sample(k) * subtree.sample(parent - k);
  let reversed_first = to_usize(to_isize(subtree.values.len()) - 1 - parent + k1 + subtree.start);
  let sum = branch
    .values
    .slice(s![to_usize(k1 - branch.start)..=to_usize(k2 - branch.start)])
    .dot(&reversed.slice(s![reversed_first..=reversed_first + to_usize(k2 - k1)]));
  let (g1, g2) = (sample(k1), sample(k2));
  sum - f64::midpoint(g1, g2)
    + 0.5 * (k1.to_f64_exact() - k_lo) * (product(k_lo) + g1)
    + 0.5 * (k_hi - k2.to_f64_exact()) * (g2 + product(k_hi))
}

fn cell_moments(
  from: f64,
  to: f64,
  center: f64,
  window: &GridWindow,
  exact: &ExactFunction<'_>,
  defect: &dyn Fn(f64) -> Result<f64, Report>,
) -> Result<(f64, f64), Report> {
  let (lo, hi) = (
    exact.breaks.partition_point(|&x| x <= from),
    exact.breaks.partition_point(|&x| x < to),
  );
  once(from)
    .chain(exact.breaks[lo..hi.max(lo)].iter().copied())
    .chain(once(to))
    .tuple_windows()
    .map(|(left, right)| {
      if let Some(moments) = log_linear_moments(left, right, center, window, exact)? {
        return Ok(moments);
      }
      let points = piece_points(left, right, exact)?;
      points
        .iter()
        .tuple_windows()
        .try_fold((0.0, 0.0), |(d0, d1), (&a, &b)| {
          let (middle, half) = (f64::midpoint(a, b), 0.5 * (b - a));
          GAUSS_LEGENDRE_8.iter().try_fold((d0, d1), |(d0, d1), &(node, weight)| {
            let x = middle + half * node;
            let value = half * weight * defect(x)?;
            Ok::<_, Report>((d0 + value, d1 + value * (x - center)))
          })
        })
    })
    .try_fold((0.0, 0.0), |(d0, d1), moments: Result<(f64, f64), Report>| {
      let (piece_d0, piece_d1) = moments?;
      Ok((d0 + piece_d0, d1 + piece_d1))
    })
}

fn log_linear_moments(
  from: f64,
  to: f64,
  center: f64,
  window: &GridWindow,
  exact: &ExactFunction<'_>,
) -> Result<Option<(f64, f64)>, Report> {
  if !exact.log_linear || starts_at_singular(from, exact) {
    return Ok(None);
  }
  let (left, right) = (
    (exact.log_value)(from)? - window.log_amp,
    (exact.log_value)(to)? - window.log_amp,
  );
  if !(left.is_finite() && right.is_finite()) {
    return Ok(None);
  }
  let (length, change, start) = (to - from, right - left, left.exp());
  let mass = start * length * exp_moment_0(change);
  let first = start * length * length * exp_moment_1(change) + (from - center) * mass;
  let (lower, upper) = (window.at(from), window.at(to));
  let linear_mass = 0.5 * length * (lower + upper);
  let linear_first = linear_mass * (f64::midpoint(from, to) - center) + (upper - lower) * length * length / 12.0;
  Ok(Some((mass - linear_mass, first - linear_first)))
}

fn exp_moment_0(change: f64) -> f64 {
  if change == 0.0 { 1.0 } else { change.exp_m1() / change }
}

fn exp_moment_1(change: f64) -> f64 {
  if change.abs() < EXP_SERIES_RANGE {
    (0..EXP_SERIES_TERMS)
      .scan(1.0, |term, n| {
        let value = *term / (n + 2).to_f64_exact();
        *term *= change / (n + 1).to_f64_exact();
        Some(value)
      })
      .sum()
  } else {
    (change.exp() * (change - 1.0) + 1.0) / (change * change)
  }
}

fn starts_at_singular(from: f64, exact: &ExactFunction<'_>) -> bool {
  exact
    .singular
    .iter()
    .any(|&point| (point - from).abs() <= SINGULAR_TOLERANCE * (1.0 + from.abs()))
}

fn piece_points(from: f64, to: f64, exact: &ExactFunction<'_>) -> Result<Vec<f64>, Report> {
  let length = to - from;
  if starts_at_singular(from, exact) {
    return Ok(
      once(from)
        .chain(
          (1..GRADING_LEVELS)
            .rev()
            .map(|level| from + length * 0.5_f64.powi(level)),
        )
        .chain(once(to))
        .collect(),
    );
  }
  let (left, right) = ((exact.log_value)(from)?, (exact.log_value)(to)?);
  let change = (right - left).abs();
  let panels = if change.is_finite() {
    (change / LOG_CHANGE_PER_PANEL)
      .floor()
      .to_usize()
      .map_or(MAX_PANELS, |panels| (panels + 1).min(MAX_PANELS))
  } else if left.is_finite() || right.is_finite() {
    MAX_PANELS
  } else {
    1
  };
  Ok(
    (0..=panels)
      .map(|panel| from + length * panel.to_f64_exact() / panels.to_f64_exact())
      .collect(),
  )
}

fn to_usize(value: isize) -> usize {
  usize::try_from(value).expect("grid offset must be non-negative")
}

fn lerp(from: f64, to: f64, fraction: f64) -> f64 {
  from + (to - from) * fraction
}
