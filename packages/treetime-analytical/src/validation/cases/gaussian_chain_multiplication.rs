use crate::gaussian::GaussianParams;
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct GaussianChainMultiplicationTestCase {
  pub name: &'static str,
  pub description: &'static str,
  pub stress_type: &'static str,
  pub analytical_caution: &'static str,
  pub slowness: f64,
  pub factors: &'static [GaussianParams],
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
}

const TWO_IDENTICAL_FACTORS: &[GaussianParams] = &[
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
];

const THREE_SHIFTED_FACTORS: &[GaussianParams] = &[
  GaussianParams {
    mu: -1.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 1.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
];

const FIVE_IDENTICAL_FACTORS: &[GaussianParams] = &[
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
];

const TEN_IDENTICAL_FACTORS: &[GaussianParams] = &[
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
];

const TWENTY_IDENTICAL_FACTORS: &[GaussianParams] = &[
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
];

const FIFTY_IDENTICAL_FACTORS: &[GaussianParams] = &[
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
];

const FIVE_MIXED_FACTORS: &[GaussianParams] = &[
  GaussianParams {
    mu: -2.0,
    sigma: 1.5,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: -0.5,
    sigma: 0.8,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.5,
    sigma: 1.2,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 2.0,
    sigma: 0.9,
    amplitude: 1.0,
  },
];

const TEN_SMALL_AMPLITUDE_FACTORS: &[GaussianParams] = &[
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 0.1,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 0.1,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 0.1,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 0.1,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 0.1,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 0.1,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 0.1,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 0.1,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 0.1,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 0.1,
  },
];

const HUNDRED_IDENTICAL_FACTORS: &[GaussianParams] = &[
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
  GaussianParams {
    mu: 0.0,
    sigma: 1.0,
    amplitude: 1.0,
  },
];

const SINGLE_FACTOR: &[GaussianParams] = &[GaussianParams {
  mu: 0.0,
  sigma: 1.0,
  amplitude: 1.0,
}];

const EMPTY_FACTORS: &[GaussianParams] = &[];

pub const GAUSSIAN_CHAIN_MULTIPLICATION_CASES: &[GaussianChainMultiplicationTestCase] = &[
  GaussianChainMultiplicationTestCase {
    name: "two_identical",
    description: "Two identical Gaussians (baseline for chain).",
    stress_type: "baseline",
    analytical_caution: "none",
    slowness: 0.0,
    factors: TWO_IDENTICAL_FACTORS,
    input_grid_domain: (-5.0, 5.0),
    input_grid_n_points: 201,
  },
  GaussianChainMultiplicationTestCase {
    name: "three_shifted",
    description: "Three Gaussians with shifted means.",
    stress_type: "translation",
    analytical_caution: "none",
    slowness: 0.1,
    factors: THREE_SHIFTED_FACTORS,
    input_grid_domain: (-6.0, 6.0),
    input_grid_n_points: 241,
  },
  GaussianChainMultiplicationTestCase {
    name: "five_identical",
    description: "Five identical Gaussians testing moderate chain.",
    stress_type: "moderate chain",
    analytical_caution: "none",
    slowness: 0.15,
    factors: FIVE_IDENTICAL_FACTORS,
    input_grid_domain: (-5.0, 5.0),
    input_grid_n_points: 201,
  },
  GaussianChainMultiplicationTestCase {
    name: "ten_identical",
    description: "Ten identical Gaussians testing underflow resistance.",
    stress_type: "underflow risk",
    analytical_caution: "product amplitude decreases rapidly",
    slowness: 0.25,
    factors: TEN_IDENTICAL_FACTORS,
    input_grid_domain: (-5.0, 5.0),
    input_grid_n_points: 201,
  },
  GaussianChainMultiplicationTestCase {
    name: "twenty_identical",
    description: "Twenty identical Gaussians for severe underflow stress.",
    stress_type: "severe underflow",
    analytical_caution: "requires log-scale handling",
    slowness: 0.4,
    factors: TWENTY_IDENTICAL_FACTORS,
    input_grid_domain: (-5.0, 5.0),
    input_grid_n_points: 201,
  },
  GaussianChainMultiplicationTestCase {
    name: "fifty_identical",
    description: "Fifty identical Gaussians for extreme underflow stress.",
    stress_type: "extreme underflow",
    analytical_caution: "plain multiplication would underflow to zero",
    slowness: 0.7,
    factors: FIFTY_IDENTICAL_FACTORS,
    input_grid_domain: (-5.0, 5.0),
    input_grid_n_points: 201,
  },
  GaussianChainMultiplicationTestCase {
    name: "five_mixed",
    description: "Five Gaussians with varied parameters.",
    stress_type: "varied params",
    analytical_caution: "none",
    slowness: 0.2,
    factors: FIVE_MIXED_FACTORS,
    input_grid_domain: (-8.0, 8.0),
    input_grid_n_points: 321,
  },
  GaussianChainMultiplicationTestCase {
    name: "small_amplitudes_chain",
    description: "Chain of Gaussians with small amplitudes.",
    stress_type: "amplitude underflow",
    analytical_caution: "requires log-scale handling",
    slowness: 0.3,
    factors: TEN_SMALL_AMPLITUDE_FACTORS,
    input_grid_domain: (-5.0, 5.0),
    input_grid_n_points: 201,
  },
  GaussianChainMultiplicationTestCase {
    name: "hundred_identical",
    description: "One hundred identical Gaussians for stress testing.",
    stress_type: "stress test",
    analytical_caution: "extreme precision requirements",
    slowness: 0.9,
    factors: HUNDRED_IDENTICAL_FACTORS,
    input_grid_domain: (-5.0, 5.0),
    input_grid_n_points: 201,
  },
  GaussianChainMultiplicationTestCase {
    name: "single_factor",
    description: "Single Gaussian (trivial case, result equals input).",
    stress_type: "edge case",
    analytical_caution: "none",
    slowness: 0.0,
    factors: SINGLE_FACTOR,
    input_grid_domain: (-5.0, 5.0),
    input_grid_n_points: 201,
  },
  GaussianChainMultiplicationTestCase {
    name: "empty_factors",
    description: "Empty factors list (edge case, returns identity).",
    stress_type: "edge case",
    analytical_caution: "trivial case",
    slowness: 0.0,
    factors: EMPTY_FACTORS,
    input_grid_domain: (-5.0, 5.0),
    input_grid_n_points: 201,
  },
];
