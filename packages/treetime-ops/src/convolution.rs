use crate::traits::ConvolveAlgo;
use eyre::{Report, WrapErr};
use ndarray::{Array1, s};
use ndarray_conv::{ConvExt, ConvFFTExt, ConvMode, PaddingMode, get_fft_processor};
use std::any::Any;
use std::cell::RefCell;
use parking_lot::RwLock;
use realfft::num_complex::Complex;
use realfft::{ComplexToReal, RealFftPlanner, RealToComplex};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::sync::{Arc, OnceLock};
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_internal_report};

const FFT_PROCESSORS_PER_THREAD: usize = 8;
const MAX_CACHED_FFT_LEN: usize = 1 << 14;

thread_local! {
  static FFT_PROCESSORS_BY_FFT_LEN: RefCell<Vec<(usize, Box<dyn Any>)>> = const { RefCell::new(Vec::new()) };
}

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
    "ndarray-conv-fft"
  }

  fn convolve(&self, dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
    convolve_fft(dx, f_values, g_values)
  }
}

pub fn convolve_fft(dx: f64, f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<Array1<f64>, Report> {
  check_operand_lengths(f_values.len(), g_values.len())?;
  let padded_len = f_values.len() + 2 * g_values.len().saturating_sub(1);
  let fft_len = fft_len(padded_len.max(g_values.len()));
  let overlap_sums = if fft_len > MAX_CACHED_FFT_LEN {
    f_values.conv_fft(g_values, ConvMode::Full, PaddingMode::Zeros)?
  } else {
    FFT_PROCESSORS_BY_FFT_LEN.with_borrow_mut(|processors| {
      let processor = processor_for_fft_len(processors, fft_len, get_fft_processor::<f64, f64>)?;
      Ok::<_, Report>(f_values.conv_fft_with_processor(g_values, ConvMode::Full, PaddingMode::Zeros, processor)?)
    })?
  };
  Ok(trapezoid_from_overlap_sums(overlap_sums, dx, f_values, g_values))
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
    let fft_len = fft_len((f_len + 2 * (g_len - 1)).max(g_len));
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
    let offset = self.g_len - 1;
    self.f_buffer.fill(0.0);
    self.f_buffer.slice_mut(s![offset..offset + self.f_len]).assign(f_values);
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
    self.f_spectrum[0].im = 0.0;
    if self.fft_len.is_multiple_of(2) {
      let last = self.f_spectrum.len() - 1;
      self.f_spectrum[last].im = 0.0;
    }
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
      .slice(s![offset..offset + self.f_len + self.g_len - 1])
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

#[expect(
  clippy::integer_division,
  reason = "the transform length of ndarray-conv steps down from a power of two by integer ratios"
)]
fn fft_len(padded_len: usize) -> usize {
  let mut len = padded_len.next_power_of_two();
  for (numerator, denominator) in [(3, 4), (5, 6)] {
    loop {
      let smaller = len / denominator * numerator;
      match smaller.cmp(&padded_len) {
        Ordering::Less => break,
        Ordering::Equal => return padded_len,
        Ordering::Greater => len = smaller,
      }
    }
  }
  len
}

fn processor_for_fft_len<P: Any>(
  processors: &mut Vec<(usize, Box<dyn Any>)>,
  fft_len: usize,
  new_processor: impl FnOnce() -> P,
) -> Result<&mut P, Report> {
  let entry = if let Some(index) = processors.iter().position(|(len, _)| *len == fft_len) {
    processors.remove(index)
  } else {
    if processors.len() >= FFT_PROCESSORS_PER_THREAD {
      processors.remove(0);
    }
    let processor: Box<dyn Any> = Box::new(new_processor());
    (fft_len, processor)
  };
  processors.push(entry);
  let index = processors.len() - 1;
  processors[index]
    .1
    .downcast_mut::<P>()
    .ok_or_else(|| make_internal_report!("FFT processor cache holds a processor of another type"))
}

fn contiguous<T>(array: &mut Array1<T>) -> Result<&mut [T], Report> {
  array
    .as_slice_mut()
    .ok_or_else(|| make_internal_report!("FFT buffer is not contiguous in memory"))
}
