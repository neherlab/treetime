use eyre::Report;
use ndarray::{Array1, Array2, aview1};
use rand::Rng;
use rand::distr::weighted::WeightedIndex;
use rand_distr::{Distribution, Open01};
use treetime_primitives::AsciiChar;
use treetime_utils::array::ndarray::max_or;
use treetime_utils::{make_error, make_report};

pub const NUCLEOTIDES: [AsciiChar; 4] = [
  AsciiChar::from_byte_unchecked(b'A'),
  AsciiChar::from_byte_unchecked(b'C'),
  AsciiChar::from_byte_unchecked(b'G'),
  AsciiChar::from_byte_unchecked(b'T'),
];

#[derive(Clone, Debug, PartialEq)]
pub struct SubstitutionModel {
  pi: Array1<f64>,
  w: Array2<f64>,
  mu: f64,
  root_states: WeightedIndex<f64>,
  jumps: Vec<Option<Jump>>,
}

impl SubstitutionModel {
  pub fn new(pi: Array1<f64>, w: Array2<f64>, mu: f64) -> Result<Self, Report> {
    let n = NUCLEOTIDES.len();
    if pi.len() != n {
      return make_error!(
        "Equilibrium frequencies (pi) need {n} entries for the states A, C, G, T, got {}",
        pi.len()
      );
    }
    if !pi.iter().all(|value| value.is_finite() && *value > 0.0) {
      return make_error!("Equilibrium frequencies (pi) must be finite and positive, got {pi}");
    }
    if w.shape() != [n, n] {
      return make_error!(
        "Exchangeabilities (w) must be a {n}x{n} matrix, got shape {:?}",
        w.shape()
      );
    }
    if !w.iter().all(|value| value.is_finite() && *value >= 0.0) {
      return make_error!("Exchangeabilities (w) must be finite and non-negative, got {w}");
    }
    if w.diag() != aview1(&[0.0; 4]) {
      return make_error!("Exchangeabilities (w) must have a zero diagonal, got {w}");
    }
    if w != w.t() {
      return make_error!("Exchangeabilities (w) must be symmetric, got {w}");
    }
    if !(mu.is_finite() && mu >= 0.0) {
      return make_error!("Substitution rate (mu) must be finite and non-negative, got {mu}");
    }

    let pi_max = max_or(&pi, 0.0);
    let pi = pi / pi_max;
    let pi_sum = pi.sum();
    let pi = pi / pi_sum;
    if !pi.iter().all(|value| value.is_finite() && *value > 0.0) {
      return make_error!("Equilibrium frequencies (pi) must be finite and positive, got {pi}");
    }
    let average_rate = pi.dot(&w.dot(&pi));
    if !(average_rate.is_finite() && average_rate > 0.0) {
      return make_error!("Average rate pi^T w pi of the exchangeabilities (w) must be positive, got {average_rate}");
    }
    let w = w / average_rate;

    let root_states =
      WeightedIndex::new(pi.iter()).map_err(|error| make_report!("Equilibrium frequencies (pi) {pi}: {error}"))?;
    let jumps = (0..n)
      .map(|state| Jump::new(state, &pi, &w, mu))
      .collect::<Result<_, _>>()?;
    Ok(Self {
      pi,
      w,
      mu,
      root_states,
      jumps,
    })
  }

  pub fn random(mu: f64, rng: &mut impl Rng) -> Result<Self, Report> {
    let n = NUCLEOTIDES.len();
    let pi = Array1::from_shape_fn(n, |_| Open01.sample(rng));
    let w: Array2<f64> = Array2::from_shape_fn((n, n), |_| Open01.sample(rng));
    let mut w = &w + &w.t();
    w.diag_mut().fill(0.0);
    Self::new(pi, w, mu)
  }

  pub fn jc69(mu: f64) -> Result<Self, Report> {
    Self::f81(Array1::from_elem(NUCLEOTIDES.len(), 0.25), mu)
  }

  pub fn f81(pi: Array1<f64>, mu: f64) -> Result<Self, Report> {
    let n = NUCLEOTIDES.len();
    let mut w = Array2::ones((n, n));
    w.diag_mut().fill(0.0);
    Self::new(pi, w, mu)
  }

  pub fn pi(&self) -> &Array1<f64> {
    &self.pi
  }

  pub fn w(&self) -> &Array2<f64> {
    &self.w
  }

  pub fn mu(&self) -> f64 {
    self.mu
  }

  pub fn sample_state(&self, rng: &mut impl Rng) -> usize {
    self.root_states.sample(rng)
  }

  pub(crate) fn jump(&self, state: usize) -> Option<&Jump> {
    self.jumps[state].as_ref()
  }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Jump {
  pub(crate) exit_rate: f64,
  pub(crate) targets: WeightedIndex<f64>,
}

impl Jump {
  fn new(state: usize, pi: &Array1<f64>, w: &Array2<f64>, mu: f64) -> Result<Option<Self>, Report> {
    let rates = mu * &w.row(state) * pi;
    let exit_rate = rates.sum();
    if exit_rate <= 0.0 {
      return Ok(None);
    }
    let targets = WeightedIndex::new(rates.iter())
      .map_err(|error| make_report!("Substitution rates {rates} from state {state}: {error}"))?;
    Ok(Some(Self { exit_rate, targets }))
  }
}
