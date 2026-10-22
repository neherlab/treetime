use crate::{EdgeCellLaw, GridEdge};
use eyre::Report;
use std::iter::successors;
use treetime_utils::make_error;

const TRUNCATION_DROP: f64 = 40.0;

const PANEL_MAX_WIDTH: f64 = 0.5;

const BISECTION_STEPS: u32 = 200;

const NEWTON_STEPS: u32 = 100;

const GAUSS_LEGENDRE_12: [(f64, f64); 6] = [
  (0.125_233_408_511_468_9, 0.249_147_045_813_402_77),
  (0.367_831_498_998_180_2, 0.233_492_536_538_354_8),
  (0.587_317_954_286_617_5, 0.203_167_426_723_065_92),
  (0.769_902_674_194_304_7, 0.160_078_328_543_346_22),
  (0.904_117_256_370_474_9, 0.106_939_325_995_318_43),
  (0.981_560_634_246_719_2, 0.047_175_336_386_511_83),
];

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeCell {
  pub law: EdgeCellLaw,
  pub edge: GridEdge,
  pub anchor: GridEdge,
}

impl EdgeCell {
  pub fn new(law: EdgeCellLaw, edge: GridEdge, anchor: GridEdge) -> Result<Self, Report> {
    if !(law.b.is_finite() && law.b >= 0.0) {
      return make_error!(
        "The law of an edge cell needs a finite non-negative exponent, got {}",
        law.b
      );
    }
    if !(edge.y.is_finite() && anchor.y.is_finite()) {
      return make_error!(
        "The law of an edge cell needs finite values at both cell ends, got {} at {} and {} at {}",
        edge.y,
        edge.t,
        anchor.y,
        anchor.t
      );
    }
    if !(law.offset.is_finite() && law.offset > 0.0 && (anchor.t - edge.t).abs() > 0.0) {
      return make_error!(
        "The law of an edge cell needs a finite positive offset of its origin beyond the cell edge {} and a cell \
         of positive width, got the offset {} and the other cell end {}",
        edge.t,
        law.offset,
        anchor.t
      );
    }
    Ok(Self { law, edge, anchor })
  }

  pub fn eval(&self, t: f64) -> f64 {
    self.eval_at_distance(self.distance(t))
  }

  pub fn distance(&self, t: f64) -> f64 {
    (t - self.edge.t).abs() + self.law.offset
  }

  pub fn position(&self, distance: f64) -> f64 {
    let inward = (self.anchor.t - self.edge.t).signum();
    self.edge.t + inward * (distance - self.law.offset)
  }

  pub fn edge_distance(&self) -> f64 {
    self.law.offset
  }

  pub fn anchor_distance(&self) -> f64 {
    self.distance(self.anchor.t)
  }

  pub fn slope(&self) -> f64 {
    let (d_edge, d_anchor) = (self.edge_distance(), self.anchor_distance());
    (self.edge.y - self.anchor.y + self.law.b * (d_edge / d_anchor).ln()) / (d_edge - d_anchor)
  }

  pub fn ln_relative_density(&self, distance: f64) -> f64 {
    self.anchor.y - self.eval_at_distance(distance)
  }

  pub fn mode(&self) -> Option<GridEdge> {
    let slope = self.slope();
    if !(slope > 0.0 && self.law.b > 0.0) {
      return None;
    }
    let distance = self.law.b / slope;
    (distance > self.edge_distance() && distance < self.anchor_distance()).then(|| GridEdge {
      t: self.position(distance),
      y: self.eval_at_distance(distance),
    })
  }

  pub fn ln_relative_mass(&self, near: f64, far: f64) -> f64 {
    if far <= near {
      return f64::NEG_INFINITY;
    }
    let slope = self.slope();
    let densest = if slope > 0.0 { self.law.b / slope } else { far }.clamp(near, far);
    let peak = self.ln_relative_density(densest);
    let mass = self.monotone_mass(densest, near, slope) + self.monotone_mass(densest, far, slope);
    peak + densest.ln() + mass.ln()
  }

  pub fn distance_at_mass(&self, from: CellEnd, ln_mass: f64) -> f64 {
    let (d_edge, d_anchor) = (self.edge_distance(), self.anchor_distance());
    let (start, end) = match from {
      CellEnd::Edge => (d_edge, d_anchor),
      CellEnd::Anchor => (d_anchor, d_edge),
    };
    if ln_mass == f64::NEG_INFINITY {
      return start;
    }
    if ln_mass >= self.ln_relative_mass(d_edge, d_anchor) {
      return end;
    }
    let excess = |d: f64| match from {
      CellEnd::Edge => self.ln_relative_mass(d_edge, d) - ln_mass,
      CellEnd::Anchor => self.ln_relative_mass(d, d_anchor) - ln_mass,
    };
    let (mut lower, mut upper) = (d_edge, d_anchor);
    let mut d = f64::midpoint(lower, upper);
    for _ in 0..NEWTON_STEPS {
      let gap = excess(d);
      let rises_past = match from {
        CellEnd::Edge => gap > 0.0,
        CellEnd::Anchor => gap < 0.0,
      };
      if rises_past {
        upper = d;
      } else {
        lower = d;
      }
      let rate = (self.ln_relative_density(d) - (gap + ln_mass)).exp();
      let step = match from {
        CellEnd::Edge => gap / rate,
        CellEnd::Anchor => -gap / rate,
      };
      if gap.abs() <= 2.0 * f64::EPSILON || upper - lower <= f64::EPSILON * upper {
        break;
      }
      let next = d - step;
      d = if next > lower && next < upper {
        next
      } else {
        f64::midpoint(lower, upper)
      };
    }
    d
  }

  pub fn position_between(&self, level: f64, from: f64, to: f64) -> f64 {
    let (d_from, d_to) = (self.distance(from), self.distance(to));
    let (y_from, y_to) = (self.eval_at_distance(d_from), self.eval_at_distance(d_to));
    if (level - y_from) * (level - y_to) >= 0.0 {
      return if (level - y_to).abs() <= (level - y_from).abs() {
        to
      } else {
        from
      };
    }
    let distance = bisect(|d| self.eval_at_distance(d) - level, d_from, d_to);
    self.position(distance)
  }

  fn monotone_mass(&self, densest: f64, end: f64, slope: f64) -> f64 {
    let b = self.law.b;
    let scaled_slope = slope * densest;
    let ln_density = |z: f64| b * z - scaled_slope * z.exp_m1();
    let z_end = ((end - densest) / densest).ln_1p();
    let z_stop = if ln_density(z_end) >= -TRUNCATION_DROP {
      z_end
    } else {
      bisect(|z| ln_density(z) + TRUNCATION_DROP, 0.0, z_end)
    };
    let direction = z_stop.signum();
    let rate = |z: f64| (b + 1.0 - scaled_slope * z.exp()).abs();
    let bend = |z: f64| scaled_slope.abs() * z.exp();
    let narrowest = 4.0 * f64::EPSILON * z_stop.abs();
    let panel_end = |z: f64| {
      let fits = |width: f64| {
        let next = z + direction * width;
        width <= narrowest || (width * rate(z).max(rate(next)) <= 1.0 && width * width * bend(z).max(bend(next)) <= 1.0)
      };
      let widest = (z_stop - z).abs().min(PANEL_MAX_WIDTH);
      let width = successors(Some(widest), |width| Some(0.5 * width))
        .find(|&width| fits(width))
        .unwrap_or(widest);
      if width >= (z_stop - z).abs() {
        z_stop
      } else {
        z + direction * width.max(narrowest)
      }
    };
    let ends: Vec<f64> = successors(Some(0.0_f64), |&z| {
      ((z_stop - z) * direction > 0.0).then(|| panel_end(z))
    })
    .collect();
    ends
      .iter()
      .zip(ends.iter().skip(1))
      .map(|(&start, &next)| gauss_legendre(|z| (ln_density(z) + z).exp(), start.min(next), start.max(next)))
      .sum()
  }

  fn eval_at_distance(&self, distance: f64) -> f64 {
    let d_anchor = self.anchor_distance();
    self.anchor.y - self.law.b * (distance / d_anchor).ln() + self.slope() * (distance - d_anchor)
  }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CellEnd {
  Edge,
  Anchor,
}

fn gauss_legendre(f: impl Fn(f64) -> f64, lower: f64, upper: f64) -> f64 {
  let (middle, half_width) = (f64::midpoint(lower, upper), 0.5 * (upper - lower));
  let sum: f64 = GAUSS_LEGENDRE_12
    .iter()
    .map(|&(node, weight)| weight * (f(middle - half_width * node) + f(middle + half_width * node)))
    .sum();
  half_width * sum
}

fn bisect(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
  let (mut lower, mut upper) = (a.min(b), a.max(b));
  let rises = f(upper) > f(lower);
  for _ in 0..BISECTION_STEPS {
    let middle = f64::midpoint(lower, upper);
    if middle <= lower || middle >= upper {
      break;
    }
    if (f(middle) < 0.0) == rises {
      lower = middle;
    } else {
      upper = middle;
    }
  }
  f64::midpoint(lower, upper)
}
