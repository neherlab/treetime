use crate::traits::ConvolveAlgo;
use eyre::Report;
use ndarray::Array1;
use ndarray_conv::{ConvExt, ConvFFTExt, ConvMode, PaddingMode, get_fft_processor};
use std::any::Any;
use std::cell::RefCell;
use std::cmp::Ordering;
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
  check_operand_lengths(f_values, g_values)?;
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
  check_operand_lengths(f_values, g_values)?;
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
  check_operand_lengths(f_values, g_values)?;
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

fn check_operand_lengths(f_values: &Array1<f64>, g_values: &Array1<f64>) -> Result<(), Report> {
  if f_values.len() < 2 || g_values.len() < 2 {
    return make_error!(
      "Trapezoid convolution needs at least 2 samples per operand, got {} and {}",
      f_values.len(),
      g_values.len()
    );
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
