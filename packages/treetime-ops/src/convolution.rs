use crate::traits::ConvolveAlgo;
use eyre::{Report, WrapErr};
use ndarray::{Array1, s};
use ndarray_conv::{ConvExt, ConvMode, PaddingMode};
use parking_lot::RwLock;
use realfft::num_complex::Complex;
use realfft::{ComplexToReal, RealFftPlanner, RealToComplex};
use std::collections::BTreeMap;
use std::iter::successors;
use std::sync::{Arc, OnceLock};
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_internal_report};

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
    let fft_len = smooth_even_len(f_len + g_len - 1);
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
    let overlap_sums = self
      .f_buffer
      .slice(s![..self.f_len + self.g_len - 1])
      .mapv(|sum| sum / fft_len);
    Ok(trapezoid_from_overlap_sums(overlap_sums, dx, f_values, g_values))
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
  let (last_f, last_g) = (f_values.len() - 1, g_values.len() - 1);
  let last_n = last_f + last_g;
  if n > last_n {
    return make_error!(
      "Convolution sample {n} is outside the full convolution of {} samples",
      last_n + 1
    );
  }
  if n == 0 || n == last_n {
    let (i, j) = if n == 0 { (0, 0) } else { (last_f, last_g) };
    return Ok(0.25 * f_values[i] * g_values[j] * dx);
  }
  let (first_i, last_i) = (n.saturating_sub(last_g), n.min(last_f));
  let sum: f64 = (first_i..=last_i).map(|i| f_values[i] * g_values[n - i]).sum();
  let ends = f64::midpoint(
    f_values[first_i] * g_values[n - first_i],
    f_values[last_i] * g_values[n - last_i],
  );
  Ok((sum - ends) * dx)
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
  overlap_sums[0] = 0.25 * f_values[0] * g_values[0];
  overlap_sums[last_n] = 0.25 * f_values[last_f] * g_values[last_g];
  overlap_sums *= dx;
  overlap_sums
}

pub(crate) fn smooth_even_len(min_len: usize) -> usize {
  let min_len = min_len.max(1);
  let powers_up_to_min = |base: usize| {
    successors(Some(1_usize), move |&power| power.checked_mul(base)).take_while(move |&power| power <= min_len)
  };
  powers_up_to_min(3)
    .flat_map(|power_of_3| powers_up_to_min(5).filter_map(move |power_of_5| power_of_3.checked_mul(power_of_5)))
    .filter(|&odd| odd <= min_len)
    .map(|odd| 2 * odd * min_len.div_ceil(2 * odd).next_power_of_two())
    .fold(usize::MAX, usize::min)
}

fn contiguous<T>(array: &mut Array1<T>) -> Result<&mut [T], Report> {
  array
    .as_slice_mut()
    .ok_or_else(|| make_internal_report!("FFT buffer is not contiguous in memory"))
}
