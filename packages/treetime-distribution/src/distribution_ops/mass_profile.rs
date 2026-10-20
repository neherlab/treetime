use crate::DistributionFunction;
use crate::policy::NegLog;
use eyre::Report;
use itertools::Itertools;
use ndarray::{Array1, Axis, s};
use std::iter::{once, repeat_n};
use std::ops::Add;
use treetime_grid::{BoundaryBehavior, Side};
use treetime_utils::array::ndarray::partition_point;
use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};
use treetime_utils::make_error;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TailMass {
  Required,
  ZeroWhereUndefined,
}

#[derive(Clone, Debug)]
pub struct MassProfile {
  y_peak: f64,
  x_min: f64,
  dx: f64,
  r: Array1<f64>,
  hard_cells: Vec<HardCellPower>,
  prefix: Array1<f64>,
  suffix: Array1<f64>,
  cells: CellSums,
  left_extrap: BoundaryBehavior,
  right_extrap: BoundaryBehavior,
  left: Option<TailLaw>,
  right: Option<TailLaw>,
}

impl MassProfile {
  pub fn new(f: &DistributionFunction<f64, NegLog>) -> Result<Self, Report> {
    let n = f.len();
    if n < 2 {
      return make_error!("Mass profile needs at least two grid points, got {n}");
    }
    let y_peak = f.grid_fn().y_min();
    if !y_peak.is_finite() {
      return make_error!("Mass profile requires a finite peak ordinate, got {y_peak}");
    }
    if f.y().iter().any(|y| y.is_nan()) {
      return make_error!("Mass profile requires ordinates without NaN");
    }
    let dx = f.dx();
    let r = f.y().mapv(|y| y - y_peak);
    let (left_extrap, right_extrap) = (f.left_extrap(), f.right_extrap());
    let hard_cells: Vec<HardCellPower> = (0..n - 1)
      .filter(|&cell| r[cell].is_infinite() != r[cell + 1].is_infinite())
      .filter_map(|cell| {
        f.grid_fn()
          .hard_cell(cell, left_extrap, right_extrap)
          .map(|hard| hard.map(|hard| HardCellPower::new(cell, hard.law.b)))
          .transpose()
      })
      .try_collect()?;
    let cells = (0..n - 1).map(|j| full_cell_mass(&r, &hard_cells, dx, j)).collect_vec();
    let masses = cells.iter().map(|cell| cell.linear()).collect_vec();
    let prefix = running_sums(masses.iter().copied());
    let mut suffix = running_sums(masses.iter().rev().copied());
    suffix.invert_axis(Axis(0));
    let left = TailLaw::new(left_extrap, Side::Left, f.x_min(), r[0]);
    let right = TailLaw::new(right_extrap, Side::Right, f.x_max(), r[n - 1]);
    Ok(Self {
      y_peak,
      x_min: f.x_min(),
      dx,
      r,
      hard_cells,
      prefix,
      suffix,
      cells: CellSums::new(cells.into_iter()),
      left_extrap,
      right_extrap,
      left,
      right,
    })
  }

  pub fn y_peak(&self) -> f64 {
    self.y_peak
  }

  pub fn neg_log_mass(&self, relative_mass: f64) -> f64 {
    self.y_peak - relative_mass.ln()
  }

  pub fn mass_between(&self, a: f64, b: f64) -> Result<f64, Report> {
    Ok(self.interval_mass(self.left, self.right, a, b)?.linear())
  }

  pub fn neg_log_mass_between(&self, a: f64, b: f64) -> Result<f64, Report> {
    Ok(self.interval_mass(self.left, self.right, a, b)?.neg_log())
  }

  pub fn fractions(&self, policy: TailMass) -> Option<MassFractions<'_>> {
    let left = resolve_tail(self.left, policy, Side::Left, self.x_min)?;
    let right = resolve_tail(self.right, policy, Side::Right, self.x_max())?;
    let left_total = left.total()?;
    let right_total = right.total()?;
    let grid_total = self.prefix[self.prefix.len() - 1];
    let total = left_total + grid_total + right_total;
    (total.is_finite() && total > 0.0).then_some(MassFractions {
      profile: self,
      left,
      right,
      left_total,
      grid_total,
      grid_total_from_right: self.suffix[0],
      right_total,
      total,
    })
  }

  pub fn x_min(&self) -> f64 {
    self.x_min
  }

  pub fn x_max(&self) -> f64 {
    self.x_at(self.r.len() - 1)
  }

  pub(crate) fn extrap(&self, side: Side) -> BoundaryBehavior {
    match side {
      Side::Left => self.left_extrap,
      Side::Right => self.right_extrap,
    }
  }

  fn x_at(&self, idx: usize) -> f64 {
    self.x_min + self.dx * idx.to_f64_exact()
  }

  fn cell(&self, j: usize) -> Cell {
    cell_at(&self.r, &self.hard_cells, j)
  }

  fn interval_mass(&self, left: Option<TailLaw>, right: Option<TailLaw>, a: f64, b: f64) -> Result<ScaledMass, Report> {
    if a.is_nan() || b.is_nan() || a > b {
      return make_error!("Mass between {a} and {b} requires ordered bounds");
    }
    let x_min = self.x_min;
    let x_max = self.x_max();
    let left_mass = if a < x_min {
      let Some(law) = left else {
        return make_error!(
          "Mass between {a} and {b} reaches below the grid edge {x_min}, where the left side has no tail law"
        );
      };
      law.mass_within(a, b)
    } else {
      ScaledMass::ZERO
    };
    let right_mass = if b > x_max {
      let Some(law) = right else {
        return make_error!(
          "Mass between {a} and {b} reaches above the grid edge {x_max}, where the right side has no tail law"
        );
      };
      law.mass_within(a, b)
    } else {
      ScaledMass::ZERO
    };
    let lo = a.max(x_min);
    let hi = b.min(x_max);
    let grid_mass = if lo < hi {
      self.grid_mass_between(lo, hi)
    } else {
      ScaledMass::ZERO
    };
    Ok(left_mass + grid_mass + right_mass)
  }

  fn locate(&self, x: f64) -> CellPosition {
    let last_cell = self.r.len() - 2;
    let j = f64_to_usize_saturating(((x - self.x_min) / self.dx).floor()).min(last_cell);
    let j = if j < last_cell && x >= self.x_at(j + 1) {
      j + 1
    } else if j > 0 && x < self.x_at(j) {
      j - 1
    } else {
      j
    };
    let offset = if x <= self.x_at(j) {
      Offset::start(self.dx)
    } else if x >= self.x_at(j + 1) {
      Offset::end(self.dx)
    } else {
      let from_left = (x - self.x_at(j)).min(self.dx);
      let from_right = (self.x_at(j + 1) - x).min(self.dx);
      if from_left <= from_right {
        Offset {
          from_left,
          from_right: self.dx - from_left,
        }
      } else {
        Offset {
          from_left: self.dx - from_right,
          from_right,
        }
      }
    };
    CellPosition { j, offset }
  }

  fn grid_mass_below(&self, x: f64) -> f64 {
    let at = self.locate(x);
    let cell = self.cell(at.j).mass_between(self.dx, Offset::start(self.dx), at.offset);
    (self.prefix[at.j] + cell.linear()).min(self.prefix[at.j + 1])
  }

  fn grid_mass_above(&self, x: f64) -> f64 {
    let at = self.locate(x);
    let cell = self.cell(at.j).mass_between(self.dx, at.offset, Offset::end(self.dx));
    (self.suffix[at.j + 1] + cell.linear()).min(self.suffix[at.j])
  }

  fn grid_mass_between(&self, a: f64, b: f64) -> ScaledMass {
    let lo = self.locate(a);
    let hi = self.locate(b);
    if lo.j == hi.j {
      return self.cell(lo.j).mass_between(self.dx, lo.offset, hi.offset);
    }
    let head = self.cell(lo.j).mass_between(self.dx, lo.offset, Offset::end(self.dx));
    let tail = self.cell(hi.j).mass_between(self.dx, Offset::start(self.dx), hi.offset);
    head + self.cells.sum(lo.j + 1, hi.j) + tail
  }

  fn grid_position_below(&self, m: f64) -> f64 {
    let j = partition_point(&self.prefix.slice(s![1..]), |&p| p < m).min(self.r.len() - 2);
    let u = self.cell(j).offset_from(self.dx, Side::Left, m - self.prefix[j]);
    self.position_in_cell(j, u, self.dx - u)
  }

  fn grid_position_above(&self, m: f64) -> f64 {
    let j = partition_point(&self.suffix, |&s| s >= m)
      .saturating_sub(1)
      .min(self.r.len() - 2);
    let u = self.cell(j).offset_from(self.dx, Side::Right, m - self.suffix[j + 1]);
    self.position_in_cell(j, self.dx - u, u)
  }

  fn position_in_cell(&self, j: usize, from_left: f64, from_right: f64) -> f64 {
    if from_left <= from_right {
      self.x_at(j) + from_left
    } else {
      self.x_at(j + 1) - from_right
    }
  }
}

#[derive(Clone, Copy, Debug)]
pub struct MassFractions<'a> {
  profile: &'a MassProfile,
  left: TailLaw,
  right: TailLaw,
  left_total: f64,
  grid_total: f64,
  grid_total_from_right: f64,
  right_total: f64,
  total: f64,
}

impl MassFractions<'_> {
  pub fn total(&self) -> f64 {
    self.total
  }

  pub fn mass_between(&self, a: f64, b: f64) -> Result<f64, Report> {
    Ok(
      self
        .profile
        .interval_mass(Some(self.left), Some(self.right), a, b)?
        .linear(),
    )
  }

  pub fn cdf(&self, x: f64) -> f64 {
    self.mass_below(x) / self.total
  }

  pub fn survival(&self, x: f64) -> f64 {
    self.mass_above(x) / self.total
  }

  pub fn quantile(&self, q: f64) -> f64 {
    self.position_below(q * self.total)
  }

  pub fn survival_inverse(&self, q: f64) -> f64 {
    self.position_above(q * self.total)
  }

  pub fn quantile_from(&self, lower: f64, q: f64) -> f64 {
    self.position_below(self.mass_below(lower) + q * self.total)
  }

  pub fn survival_inverse_from(&self, upper: f64, q: f64) -> f64 {
    self.position_above(self.mass_above(upper) + q * self.total)
  }

  fn mass_below(&self, x: f64) -> f64 {
    let x_min = self.profile.x_min;
    let x_max = self.profile.x_max();
    if x <= x_min {
      self.left.mass_within(f64::NEG_INFINITY, x).linear()
    } else if x >= x_max {
      self.left_total + self.grid_total + self.right.mass_within(x_max, x).linear()
    } else {
      self.left_total + self.profile.grid_mass_below(x)
    }
  }

  fn mass_above(&self, x: f64) -> f64 {
    let x_min = self.profile.x_min;
    let x_max = self.profile.x_max();
    if x >= x_max {
      self.right.mass_within(x, f64::INFINITY).linear()
    } else if x <= x_min {
      self.right_total + self.grid_total_from_right + self.left.mass_within(x, x_min).linear()
    } else {
      self.right_total + self.profile.grid_mass_above(x)
    }
  }

  fn position_below(&self, m: f64) -> f64 {
    if m < self.left_total {
      return self.left.position_at_outer_mass(m);
    }
    let m = m - self.left_total;
    if m <= self.grid_total {
      return self.profile.grid_position_below(m);
    }
    self.right.position_at_inner_mass(m - self.grid_total)
  }

  fn position_above(&self, m: f64) -> f64 {
    let right_total = self.right_total;
    if m < right_total {
      return self.right.position_at_outer_mass(m);
    }
    let m = m - right_total;
    if m <= self.grid_total_from_right {
      return self.profile.grid_position_above(m);
    }
    self.left.position_at_inner_mass(m - self.grid_total_from_right)
  }
}

#[derive(Clone, Copy, Debug)]
struct HardCellPower {
  cell: usize,
  power: f64,
}

impl HardCellPower {
  fn new(cell: usize, b: f64) -> Self {
    Self { cell, power: b + 1.0 }
  }
}

#[derive(Clone, Copy, Debug)]
enum Cell {
  Trapezoid { ra: f64, rb: f64 },
  PowerLaw { r_anchor: f64, power: f64, hard_end: Side },
  Empty,
}

impl Cell {
  fn mass_between(self, h: f64, lo: Offset, hi: Offset) -> ScaledMass {
    let width = if lo.from_left + hi.from_left <= lo.from_right + hi.from_right {
      hi.from_left - lo.from_left
    } else {
      lo.from_right - hi.from_right
    };
    if width <= 0.0 {
      return ScaledMass::ZERO;
    }
    match self {
      Cell::Empty => ScaledMass::ZERO,
      Cell::Trapezoid { ra, rb } => {
        let r_min = ra.min(rb);
        let (qa, qb) = ((r_min - ra).exp(), (r_min - rb).exp());
        let density = |at: Offset| {
          if at.from_left <= at.from_right {
            qa + (qb - qa) * (at.from_left / h)
          } else {
            qb + (qa - qb) * (at.from_right / h)
          }
        };
        ScaledMass::new(r_min, width * f64::midpoint(density(lo), density(hi)))
      },
      Cell::PowerLaw {
        r_anchor,
        power,
        hard_end,
      } => {
        let distance = |at: Offset| match hard_end {
          Side::Left => at.from_left,
          Side::Right => at.from_right,
        };
        power_law_mass(r_anchor, h, power, distance(lo).max(distance(hi)), width)
      },
    }
  }

  fn offset_from(self, h: f64, end: Side, m: f64) -> f64 {
    if m <= 0.0 {
      return 0.0;
    }
    match self {
      Cell::Empty => 0.0,
      Cell::Trapezoid { ra, rb } => {
        let (r_e, r_f) = match end {
          Side::Left => (ra, rb),
          Side::Right => (rb, ra),
        };
        let r_min = r_e.min(r_f);
        let (q_e, q_f) = ((r_min - r_e).exp(), (r_min - r_f).exp());
        let target = scaled_by_density(m, r_min);
        if target >= h * f64::midpoint(q_e, q_f) {
          return h;
        }
        let slope = (q_f - q_e) / h;
        let root = (q_e * q_e + 2.0 * slope * target).max(0.0).sqrt();
        (2.0 * target / (q_e + root)).clamp(0.0, h)
      },
      Cell::PowerLaw {
        r_anchor,
        power,
        hard_end,
      } => {
        let fraction = (scaled_by_density(m, r_anchor) * power / h).min(1.0);
        if end == hard_end {
          power_law_distance_from_hard_point(fraction, h, power)
        } else {
          power_law_distance_from_anchor(fraction, h, power)
        }
      },
    }
  }
}

#[derive(Clone, Copy, Debug)]
struct Offset {
  from_left: f64,
  from_right: f64,
}

impl Offset {
  fn start(h: f64) -> Self {
    Self {
      from_left: 0.0,
      from_right: h,
    }
  }

  fn end(h: f64) -> Self {
    Self {
      from_left: h,
      from_right: 0.0,
    }
  }
}

struct CellPosition {
  j: usize,
  offset: Offset,
}

#[derive(Clone, Copy, Debug)]
struct TailLaw {
  side: Side,
  x_edge: f64,
  r_edge: f64,
  kind: TailKind,
}

impl TailLaw {
  fn new(behavior: BoundaryBehavior, side: Side, x_edge: f64, r_edge: f64) -> Option<Self> {
    let kind = match behavior {
      BoundaryBehavior::Error => return None,
      BoundaryBehavior::Hard => TailKind::Zero,
      BoundaryBehavior::Linear(law) => TailKind::Exponential {
        rate: match side {
          Side::Left => -law.slope,
          Side::Right => law.slope,
        },
      },
      BoundaryBehavior::HardApproach(law) => TailKind::PowerLaw {
        t_hard: law.t_hard,
        delta_edge: (x_edge - law.t_hard).abs(),
        power: law.b + 1.0,
      },
    };
    Some(Self {
      side,
      x_edge,
      r_edge,
      kind,
    })
  }

  fn zero(side: Side, x_edge: f64) -> Self {
    Self {
      side,
      x_edge,
      r_edge: f64::INFINITY,
      kind: TailKind::Zero,
    }
  }

  fn total(self) -> Option<f64> {
    let total = self.total_mass();
    total.is_finite().then_some(total)
  }

  fn total_mass(self) -> f64 {
    self.mass_within(f64::NEG_INFINITY, f64::INFINITY).linear()
  }

  fn mass_within(self, a: f64, b: f64) -> ScaledMass {
    let (near, far) = match self.side {
      Side::Left => (b.min(self.x_edge), a),
      Side::Right => (a.max(self.x_edge), b),
    };
    let width = self.outward(near, far);
    if self.r_edge.is_infinite() || width <= 0.0 {
      return ScaledMass::ZERO;
    }
    match self.kind {
      TailKind::Zero => ScaledMass::ZERO,
      TailKind::Exponential { rate } => {
        let r = self.r_edge + rate * self.outward(self.x_edge, near);
        if rate == 0.0 {
          return ScaledMass::new(r, width);
        }
        if rate < 0.0 {
          return growing_exponential_mass(r, -rate, width);
        }
        let decay = -(-rate * width).exp_m1();
        let factor = decay / rate;
        if factor.is_finite() {
          ScaledMass::new(r, factor)
        } else {
          ScaledMass::new(r - (decay.ln() - rate.ln()), 1.0)
        }
      },
      TailKind::PowerLaw {
        t_hard,
        delta_edge,
        power,
      } => {
        let delta_near = self.outward(near, t_hard).clamp(0.0, delta_edge);
        let delta_far = self.outward(far, t_hard).clamp(0.0, delta_edge);
        power_law_mass(self.r_edge, delta_edge, power, delta_near, delta_near - delta_far)
      },
    }
  }

  fn position_at_outer_mass(self, m: f64) -> f64 {
    match self.kind {
      TailKind::Zero => self.x_edge,
      TailKind::Exponential { rate } => self.at_distance_beyond_edge((-self.r_edge - rate.ln() - m.ln()) / rate),
      TailKind::PowerLaw {
        t_hard,
        delta_edge,
        power,
      } => {
        let fraction = (m / self.total_mass()).min(1.0);
        self.at_distance_from_hard_point(t_hard, power_law_distance_from_hard_point(fraction, delta_edge, power))
      },
    }
  }

  fn position_at_inner_mass(self, m: f64) -> f64 {
    match self.kind {
      TailKind::Zero => self.x_edge,
      TailKind::Exponential { rate } => {
        let scaled = m * rate * self.r_edge.exp();
        let scaled = if scaled.is_finite() {
          scaled
        } else {
          (m.ln() + rate.ln() + self.r_edge).exp()
        };
        self.at_distance_beyond_edge(-(-scaled).max(-1.0).ln_1p() / rate)
      },
      TailKind::PowerLaw { delta_edge, power, .. } => {
        let fraction = (m / self.total_mass()).min(1.0);
        self.at_distance_beyond_edge(power_law_distance_from_anchor(fraction, delta_edge, power))
      },
    }
  }

  fn outward(self, from: f64, to: f64) -> f64 {
    match self.side {
      Side::Left => from - to,
      Side::Right => to - from,
    }
  }

  fn at_distance_beyond_edge(self, distance: f64) -> f64 {
    match self.side {
      Side::Left => self.x_edge - distance,
      Side::Right => self.x_edge + distance,
    }
  }

  fn at_distance_from_hard_point(self, t_hard: f64, delta: f64) -> f64 {
    match self.side {
      Side::Left => t_hard + delta,
      Side::Right => t_hard - delta,
    }
  }
}

#[derive(Clone, Copy, Debug)]
enum TailKind {
  Zero,
  Exponential { rate: f64 },
  PowerLaw { t_hard: f64, delta_edge: f64, power: f64 },
}

#[derive(Clone, Copy, Debug)]
struct ScaledMass {
  r: f64,
  factor: f64,
}

impl ScaledMass {
  const ZERO: Self = Self {
    r: f64::INFINITY,
    factor: 0.0,
  };

  const INFINITE: Self = Self {
    r: f64::NEG_INFINITY,
    factor: 1.0,
  };

  fn new(r: f64, factor: f64) -> Self {
    if factor > 0.0 && r < f64::INFINITY {
      Self { r, factor }
    } else {
      Self::ZERO
    }
  }

  fn linear(self) -> f64 {
    if self.factor == 0.0 {
      return 0.0;
    }
    let scale = (-self.r).exp();
    if scale.is_normal() {
      self.factor * scale
    } else {
      (self.factor.ln() - self.r).exp()
    }
  }

  fn neg_log(self) -> f64 {
    self.r - self.factor.ln()
  }

  fn ln_mass(self) -> f64 {
    -self.neg_log()
  }
}

impl Add for ScaledMass {
  type Output = Self;

  fn add(self, other: Self) -> Self {
    if self.factor == 0.0 || other.ln_mass() == f64::INFINITY {
      return other;
    }
    if other.factor == 0.0 || self.ln_mass() == f64::INFINITY {
      return self;
    }
    let r = self.r.min(other.r);
    let (scale, other_scale) = ((r - self.r).exp(), (r - other.r).exp());
    if scale.is_normal() && other_scale.is_normal() {
      return Self {
        r,
        factor: self.factor * scale + other.factor * other_scale,
      };
    }
    let (large, small) = if self.ln_mass() >= other.ln_mass() {
      (self, other)
    } else {
      (other, self)
    };
    Self {
      r: large.r,
      factor: large.factor * (1.0 + (small.ln_mass() - large.ln_mass()).exp()),
    }
  }
}

#[derive(Clone, Debug)]
struct CellSums {
  nodes: Vec<ScaledMass>,
}

impl CellSums {
  fn new(cells: impl ExactSizeIterator<Item = ScaledMass>) -> Self {
    let n = cells.len();
    let mut nodes = repeat_n(ScaledMass::ZERO, n).chain(cells).collect_vec();
    for i in (1..n).rev() {
      nodes[i] = nodes[2 * i] + nodes[2 * i + 1];
    }
    Self { nodes }
  }

  fn sum(&self, from: usize, to: usize) -> ScaledMass {
    let n = self.nodes.len() >> 1;
    let (mut lo, mut hi) = (from + n, to + n);
    let mut sum = ScaledMass::ZERO;
    while lo < hi {
      if lo & 1 == 1 {
        sum = sum + self.nodes[lo];
        lo += 1;
      }
      if hi & 1 == 1 {
        hi -= 1;
        sum = sum + self.nodes[hi];
      }
      lo >>= 1;
      hi >>= 1;
    }
    sum
  }
}

fn cell_at(r: &Array1<f64>, hard_cells: &[HardCellPower], j: usize) -> Cell {
  let (ra, rb) = (r[j], r[j + 1]);
  let power = || {
    hard_cells
      .binary_search_by_key(&j, |hard| hard.cell)
      .map_or(1.0, |i| hard_cells[i].power)
  };
  match (ra.is_infinite(), rb.is_infinite()) {
    (false, false) => Cell::Trapezoid { ra, rb },
    (true, false) => Cell::PowerLaw {
      r_anchor: rb,
      power: power(),
      hard_end: Side::Left,
    },
    (false, true) => Cell::PowerLaw {
      r_anchor: ra,
      power: power(),
      hard_end: Side::Right,
    },
    (true, true) => Cell::Empty,
  }
}

fn full_cell_mass(r: &Array1<f64>, hard_cells: &[HardCellPower], dx: f64, j: usize) -> ScaledMass {
  cell_at(r, hard_cells, j).mass_between(dx, Offset::start(dx), Offset::end(dx))
}

fn resolve_tail(law: Option<TailLaw>, policy: TailMass, side: Side, x_edge: f64) -> Option<TailLaw> {
  match (law.filter(|law| law.total().is_some()), policy) {
    (Some(law), _) => Some(law),
    (None, TailMass::Required) => None,
    (None, TailMass::ZeroWhereUndefined) => Some(TailLaw::zero(side, x_edge)),
  }
}

fn running_sums(masses: impl Iterator<Item = f64>) -> Array1<f64> {
  once(0.0)
    .chain(masses.scan((0.0_f64, 0.0_f64), |(sum, compensation), m| {
      let next = *sum + m;
      *compensation += if sum.abs() >= m.abs() {
        (*sum - next) + m
      } else {
        (m - next) + *sum
      };
      *sum = next;
      Some(*sum + *compensation)
    }))
    .collect()
}

fn growing_exponential_mass(r_near: f64, growth: f64, width: f64) -> ScaledMass {
  if width == f64::INFINITY {
    return ScaledMass::INFINITE;
  }
  let factor = (growth * width).exp_m1() / growth;
  if factor.is_finite() {
    ScaledMass::new(r_near, factor)
  } else {
    let ln_factor = growth * width + (-(-growth * width).exp()).ln_1p() - growth.ln();
    ScaledMass::new(r_near - ln_factor, 1.0)
  }
}

fn power_law_mass(r_edge: f64, delta_edge: f64, power: f64, delta_near: f64, width: f64) -> ScaledMass {
  if delta_near <= 0.0 || width <= 0.0 {
    return ScaledMass::ZERO;
  }
  let log_ratio = (-width.min(delta_near) / delta_near).ln_1p();
  let shape = delta_edge / power * -(power * log_ratio).exp_m1();
  let ratio = delta_near / delta_edge;
  let scale = ratio.powf(power);
  if scale.is_normal() {
    ScaledMass::new(r_edge, scale * shape)
  } else {
    ScaledMass::new(r_edge - power * ratio.ln(), shape)
  }
}

fn power_law_distance_from_hard_point(fraction: f64, delta_edge: f64, power: f64) -> f64 {
  (delta_edge * fraction.powf(power.recip())).min(delta_edge)
}

fn power_law_distance_from_anchor(fraction: f64, delta_edge: f64, power: f64) -> f64 {
  (delta_edge * -((-fraction).ln_1p() / power).exp_m1()).clamp(0.0, delta_edge)
}

fn scaled_by_density(m: f64, r: f64) -> f64 {
  let scaled = m * r.exp();
  if scaled.is_finite() { scaled } else { (m.ln() + r).exp() }
}
