use crate::traits::ConvolveAlgo;
use eyre::{Report, WrapErr};
use itertools::{Itertools, izip};
use ndarray::{Array1, Zip, s};
use ndarray_conv::{ConvExt, ConvMode, PaddingMode};
use parking_lot::RwLock;
use realfft::num_complex::Complex;
use realfft::{ComplexToReal, RealFftPlanner, RealToComplex};
use std::collections::BTreeMap;
use std::f64::consts::{SQRT_2, TAU};
use std::ops::RangeInclusive;
use std::sync::{Arc, OnceLock};
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_internal_report, make_report};

const SINGLE_PAIR_WEIGHT: f64 = 0.25;

const NEGLIGIBLE_TERM_DEPTH: f64 = 50.0;

pub const UNIT_ROUND_OFF: f64 = f64::EPSILON / 2.0;

#[derive(Debug)]
pub struct DirectConvolve;

impl ConvolveAlgo for DirectConvolve {
  fn name(&self) -> &'static str {
    "direct"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve_direct(dx, f_values, g_values)
  }
}

fn convolve_direct(dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
  check_operand_lengths(f_values.len(), g_values.len())?;
  let mut overlap_sums = Array1::zeros(f_values.len() + g_values.len() - 1);
  for (i, &f_val) in f_values.iter().enumerate() {
    for (j, &g_val) in g_values.iter().enumerate() {
      overlap_sums[i + j] += f_val * g_val;
    }
  }
  Ok(trapezoid_from_overlap_sums(overlap_sums, dx, f_values, g_values))
}

#[derive(Debug)]
pub struct NdarrayConvolve;

impl ConvolveAlgo for NdarrayConvolve {
  fn name(&self) -> &'static str {
    "ndarray-conv"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve(dx, f_values, g_values)
  }
}

fn convolve(dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
  check_operand_lengths(f_values.len(), g_values.len())?;
  let overlap_sums = f_values.conv(g_values, ConvMode::Full, PaddingMode::Zeros)?;
  Ok(trapezoid_from_overlap_sums(overlap_sums, dx, f_values, g_values))
}

#[derive(Debug)]
pub struct FftConvolve;

impl ConvolveAlgo for FftConvolve {
  fn name(&self) -> &'static str {
    "fft"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve_fft(dx, f_values, g_values)
  }
}

pub(crate) fn convolve_fft(dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
  FftConvolver::new(f_values.len(), g_values.len(), &FftPlans::new())?.convolve(dx, f_values, g_values)
}

pub struct FftConvolver {
  f_len: usize,
  g_len: usize,
  fft_len: usize,
  plans: FftPlanPair,
  f_buffer: Array1<f64>,
  g_buffer: Array1<f64>,
  f_spectrum: Array1<Complex<f64>>,
  g_spectrum: Array1<Complex<f64>>,
  scratch: Array1<Complex<f64>>,
}

impl FftConvolver {
  pub fn new(f_len: usize, g_len: usize, plans: &FftPlans) -> Result<Self, Report> {
    check_operand_lengths(f_len, g_len)?;
    let fft_len = power_of_two_len(f_len, g_len)?;
    let plans = plans.get(fft_len);
    let spectrum_len = plans.forward.complex_len();
    let scratch_len = plans.forward.get_scratch_len().max(plans.inverse.get_scratch_len());
    Ok(Self {
      f_len,
      g_len,
      fft_len,
      plans,
      f_buffer: Array1::zeros(fft_len),
      g_buffer: Array1::zeros(fft_len),
      f_spectrum: Array1::zeros(spectrum_len),
      g_spectrum: Array1::zeros(spectrum_len),
      scratch: Array1::zeros(scratch_len),
    })
  }

  pub fn convolve(&mut self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let overlap_sums = self.overlap_sums(f_values, g_values)?;
    Ok(trapezoid_from_overlap_sums(overlap_sums, dx, f_values, g_values))
  }

  pub fn overlap_sums(&mut self, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    if (f_values.len(), g_values.len()) != (self.f_len, self.g_len) {
      return make_error!(
        "FFT convolver for operands of {} and {} samples got operands of {} and {} samples",
        self.f_len,
        self.g_len,
        f_values.len(),
        g_values.len()
      );
    }
    self.f_buffer.fill(0.0);
    self.f_buffer.slice_mut(s![..self.f_len]).assign(f_values);
    self.g_buffer.fill(0.0);
    self.g_buffer.slice_mut(s![..self.g_len]).assign(g_values);
    self
      .plans
      .forward
      .process_with_scratch(
        contiguous(&mut self.f_buffer)?,
        contiguous(&mut self.f_spectrum)?,
        contiguous(&mut self.scratch)?,
      )
      .wrap_err("When computing the forward FFT of the first convolution operand")?;
    self
      .plans
      .forward
      .process_with_scratch(
        contiguous(&mut self.g_buffer)?,
        contiguous(&mut self.g_spectrum)?,
        contiguous(&mut self.scratch)?,
      )
      .wrap_err("When computing the forward FFT of the second convolution operand")?;
    self.f_spectrum.zip_mut_with(&self.g_spectrum, |f, g| *f *= *g);
    let last = self.f_spectrum.len() - 1;
    self.f_spectrum[0].im = 0.0;
    self.f_spectrum[last].im = 0.0;
    self
      .plans
      .inverse
      .process_with_scratch(
        contiguous(&mut self.f_spectrum)?,
        contiguous(&mut self.f_buffer)?,
        contiguous(&mut self.scratch)?,
      )
      .wrap_err("When computing the inverse FFT of the convolution product")?;
    let fft_len = self.fft_len.to_f64_exact();
    Ok(
      self
        .f_buffer
        .slice(s![..self.f_len + self.g_len - 1])
        .mapv(|sum| sum / fft_len),
    )
  }

  pub fn overlap_sums_round_off(&self, f_values: &Array1<f64>, g_values: &Array1<f64>) -> f64 {
    let fft_len = self.fft_len.to_f64_exact();
    let twiddle = TAU * rounding_gamma(2) + SQRT_2 * UNIT_ROUND_OFF;
    let level = twiddle + rounding_gamma(4) * (SQRT_2 + twiddle);
    let levels = fft_len.log2() - 1.0;
    let complex_transform = levels * level / (1.0 - levels * level);
    let packing = 10_f64.sqrt() * (rounding_gamma(4) + (1.0 + rounding_gamma(4)) * twiddle);
    let real_transform = complex_transform + packing + complex_transform * packing;
    let product = SQRT_2 * rounding_gamma(2);
    let norm_1 = |x: &Array1<f64>| x.iter().map(|v| v.abs()).sum::<f64>();
    let norm_2 = |x: &Array1<f64>| x.dot(x).sqrt();
    let (f_1, f_2, g_1, g_2) = (norm_1(f_values), norm_2(f_values), norm_1(g_values), norm_2(g_values));
    let spectrum_error = |norm_2: f64| real_transform * fft_len.sqrt() * norm_2;
    let forward_and_product =
      (2.0 * real_transform + real_transform * real_transform + product * (1.0 + real_transform).powi(2)) * f_2 * g_2;
    let inverse = real_transform
      * (1.0 + real_transform)
      * (1.0 + product)
      * (f_2 * (g_1 + spectrum_error(g_2))).min((f_1 + spectrum_error(f_2)) * g_2);
    let underflow = fft_len * fft_len * f64::MIN_POSITIVE * (1.0 + f_2) * (1.0 + g_2);
    (forward_and_product + inverse + underflow) / (1.0 - rounding_gamma(2 * self.fft_len + 8))
  }
}

#[derive(Default)]
pub struct FftPlans {
  by_len: RwLock<BTreeMap<usize, Arc<OnceLock<FftPlanPair>>>>,
}

impl FftPlans {
  pub fn new() -> Self {
    Self::default()
  }

  pub(crate) fn get(&self, fft_len: usize) -> FftPlanPair {
    self.slot(fft_len).get_or_init(|| FftPlanPair::plan(fft_len)).clone()
  }

  fn slot(&self, fft_len: usize) -> Arc<OnceLock<FftPlanPair>> {
    if let Some(slot) = self.by_len.read().get(&fft_len) {
      return Arc::clone(slot);
    }
    Arc::clone(self.by_len.write().entry(fft_len).or_default())
  }
}

#[derive(Clone)]
pub(crate) struct FftPlanPair {
  pub(crate) forward: Arc<dyn RealToComplex<f64>>,
  pub(crate) inverse: Arc<dyn ComplexToReal<f64>>,
}

impl FftPlanPair {
  fn plan(fft_len: usize) -> Self {
    let mut planner = RealFftPlanner::<f64>::new();
    let forward = planner.plan_fft_forward(fft_len);
    let inverse = planner.plan_fft_inverse(fft_len);
    Self { forward, inverse }
  }
}

pub fn convolve_direct_at(dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>, n: usize) -> Result<f64, Report> {
  check_operand_lengths(f_values.len(), g_values.len())?;
  let (first_i, last_i) = overlap_of_sample(f_values.len(), g_values.len(), n)?;
  if first_i == last_i {
    return Ok(SINGLE_PAIR_WEIGHT * f_values[first_i] * g_values[n - first_i] * dx);
  }
  let sum: f64 = (first_i..=last_i).map(|i| f_values[i] * g_values[n - i]).sum();
  let ends = f64::midpoint(
    f_values[first_i] * g_values[n - first_i],
    f_values[last_i] * g_values[n - last_i],
  );
  Ok((sum - ends) * dx)
}

pub fn convolve_direct_neg_log_at(
  dx: f64,
  f_neg_log: &Array1<f64>,
  g_neg_log: &Array1<f64>,
  n: usize,
  pairs: RangeInclusive<usize>,
) -> Result<f64, Report> {
  check_operand_lengths(f_neg_log.len(), g_neg_log.len())?;
  let (full_first, full_last) = overlap_of_sample(f_neg_log.len(), g_neg_log.len(), n)?;
  let (first_i, last_i) = (full_first.max(*pairs.start()), full_last.min(*pairs.end()));
  if first_i > last_i {
    return Ok(f64::INFINITY);
  }
  let neg_ln_weight = |i: usize| -trapezoid_weight(f_neg_log.len(), g_neg_log.len(), n, i).ln();
  let (first_penalty, last_penalty) = (neg_ln_weight(full_first), neg_ln_weight(full_last));
  let penalty = |q: usize| match first_i + q {
    i if i == full_first => first_penalty,
    i if i == full_last => last_penalty,
    _ => 0.0,
  };
  let f_pairs = f_neg_log.slice(s![first_i..=last_i]);
  let g_pairs = g_neg_log.slice(s![n - last_i..=n - first_i;-1]);
  if let Some((q, _)) = izip!(&f_pairs, &g_pairs).find_position(|(f, g)| f.is_nan() || g.is_nan()) {
    let i = first_i + q;
    return make_error!(
      "Convolution sample {n} has a NaN operand value in its pair ({i}, {})",
      n - i
    );
  }
  let pairs = || Zip::indexed(&f_pairs).and(&g_pairs);
  let min_cost = pairs().fold(f64::INFINITY, |min, q, f, g| min.min(f + g + penalty(q)));
  if !min_cost.is_finite() {
    return Ok(f64::INFINITY);
  }
  let sum = pairs().fold(0.0, |sum, q, f, g| {
    let depth = f + g + penalty(q) - min_cost;
    if depth <= NEGLIGIBLE_TERM_DEPTH {
      sum + (-depth).exp()
    } else {
      sum
    }
  });
  Ok(min_cost - (sum * dx).ln())
}

pub fn trapezoid_weight(f_len: usize, g_len: usize, n: usize, i: usize) -> f64 {
  let (first_i, last_i) = (n.saturating_sub(g_len - 1), n.min(f_len - 1));
  if first_i == last_i {
    SINGLE_PAIR_WEIGHT
  } else if i == first_i || i == last_i {
    0.5
  } else {
    1.0
  }
}

fn overlap_of_sample(f_len: usize, g_len: usize, n: usize) -> Result<(usize, usize), Report> {
  let last_n = f_len + g_len - 2;
  if n > last_n {
    return make_error!(
      "Convolution sample {n} is outside the full convolution of {} samples",
      last_n + 1
    );
  }
  Ok((n.saturating_sub(g_len - 1), n.min(f_len - 1)))
}

fn check_operand_lengths(f_len: usize, g_len: usize) -> Result<(), Report> {
  if f_len < 2 || g_len < 2 {
    return make_error!("Trapezoid convolution needs at least 2 samples per operand, got {f_len} and {g_len}");
  }
  Ok(())
}

pub(crate) fn trapezoid_from_overlap_sums(
  mut overlap_sums: Array1<f64>,
  dx: f64,
  f_values: &Array1<f64>,
  g_values: &Array1<f64>,
) -> Array1<f64> {
  let (last_f, last_g) = (f_values.len() - 1, g_values.len() - 1);
  let last_n = last_f + last_g;
  for (n, sum) in overlap_sums.iter_mut().enumerate().take(last_n).skip(1) {
    let (first_i, last_i) = (n.saturating_sub(last_g), n.min(last_f));
    *sum -= f64::midpoint(
      f_values[first_i] * g_values[n - first_i],
      f_values[last_i] * g_values[n - last_i],
    );
  }
  overlap_sums[0] = SINGLE_PAIR_WEIGHT * f_values[0] * g_values[0];
  overlap_sums[last_n] = SINGLE_PAIR_WEIGHT * f_values[last_f] * g_values[last_g];
  overlap_sums *= dx;
  overlap_sums
}

pub(crate) fn power_of_two_len(f_len: usize, g_len: usize) -> Result<usize, Report> {
  f_len
    .checked_add(g_len - 1)
    .and_then(usize::checked_next_power_of_two)
    .ok_or_else(|| {
      make_report!("FFT convolution of operands of {f_len} and {g_len} samples needs a transform length beyond usize")
    })
}

pub fn rounding_gamma(n: usize) -> f64 {
  let nu = n.to_f64_exact() * UNIT_ROUND_OFF;
  nu / (1.0 - nu)
}

fn contiguous<T>(array: &mut Array1<T>) -> Result<&mut [T], Report> {
  array
    .as_slice_mut()
    .ok_or_else(|| make_internal_report!("FFT buffer is not contiguous in memory"))
}
