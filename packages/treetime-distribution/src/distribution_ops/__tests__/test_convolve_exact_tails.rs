#[cfg(test)]
mod tests {
  use self::helpers::{hard_operand, mixture_closed_form, traced_convolution};
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::cast::ToF64Exact;

  #[test]
  fn test_convolution_of_hard_operands_with_steep_corners_is_exact_on_the_whole_support() -> Result<(), Report> {
    let a = hard_operand(0.0, 0.05, 201, |i, t| {
      2.0 * (t - 3.0).powi(2) + 70.0 * i.saturating_sub(192).to_f64_exact()
    })?;
    let b = hard_operand(0.0, 0.05, 121, |i, t| {
      3.0 * t + 65.0 * 6_usize.saturating_sub(i).to_f64_exact()
    })?;
    let traced = traced_convolution(&a, &b)?;
    traced.assert_covers_the_whole_support();
    assert_eq!(Vec::<usize>::new(), traced.infinite_samples());
    traced.assert_matches_the_quadrature_oracle();
    Ok(())
  }

  #[test]
  fn test_convolution_of_hard_operands_whose_tail_flattens_beyond_the_trusted_range_is_exact_on_the_whole_support()
  -> Result<(), Report> {
    let a = hard_operand(0.0, 0.05, 241, |_, t| (8.0 * (t - 1.0).powi(2)).min(70.0 + 0.2 * t))?;
    let b = hard_operand(0.0, 0.05, 61, |_, t| 4.0 * (t - 1.0).powi(2))?;
    let traced = traced_convolution(&a, &b)?;
    traced.assert_covers_the_whole_support();
    assert_eq!(Vec::<usize>::new(), traced.infinite_samples());
    traced.assert_matches_the_quadrature_oracle();
    Ok(())
  }

  #[test]
  fn test_convolution_of_truncated_gaussians_with_a_far_mode_matches_the_closed_form_on_the_whole_support()
  -> Result<(), Report> {
    let dx = 0.05;
    let a = hard_operand(-6.0, dx, 741, |_, s| {
      let (main, far) = (s * s / 2.0, 200.0 + (s - 25.0).powi(2) / 2.0);
      main.min(far) - (-(main - far).abs()).exp().ln_1p()
    })?;
    let b = hard_operand(-6.0, dx, 241, |_, s| s * s / 2.0)?;
    let traced = traced_convolution(&a, &b)?;
    traced.assert_covers_the_whole_support();
    traced.assert_matches_the_quadrature_oracle();

    let last = traced.result_len() - 1;
    let failures = traced
      .samples()
      .filter(|&(k, _)| k != 0 && k != last)
      .filter_map(|(k, _)| {
        let t = -12.0 + k.to_f64_exact() * dx;
        let exact = mixture_closed_form(t, dx);
        let oracle = traced.oracle(k);
        let tolerance = exact.discretization_bound + 1e-6;
        ((oracle - exact.neg_log).abs() > tolerance).then_some((k, oracle, exact.neg_log, tolerance))
      })
      .collect_vec();
    assert_eq!(Vec::<(usize, f64, f64, f64)>::new(), failures);
    Ok(())
  }

  #[test]
  fn test_convolution_of_a_subtree_message_with_a_plateau_tail_and_a_branch_is_exact_on_the_whole_support()
  -> Result<(), Report> {
    let (len_a, len_b) = (1163, 4321);
    let a = hard_operand(0.0, 0.026, len_a, |i, _| {
      let d = i.to_f64_exact() / len_a.to_f64_exact() - 0.264;
      if d < 0.0 {
        700.0 * d * d
      } else {
        1050.0 * d * d + 8000.0 * (d - 0.5).max(0.0).powi(2)
      }
    })?;
    let b = hard_operand(0.0, 0.026, len_b, |j, _| {
      let v = j.to_f64_exact() / len_b.to_f64_exact() - 0.2815;
      if v < 0.0 {
        115.0 * (v / 0.2815).powi(2)
      } else {
        -139.0 * (-(v / 0.27).powi(2)).exp_m1()
      }
    })?;
    let traced = traced_convolution(&a, &b)?;
    traced.assert_covers_the_whole_support();
    assert_eq!(Vec::<usize>::new(), traced.infinite_samples());
    traced.assert_matches_the_quadrature_oracle();
    Ok(())
  }

  #[test]
  fn test_convolution_window_of_a_gap_keeps_a_pair_of_narrow_components_that_its_ends_and_middle_do_not_reach()
  -> Result<(), Report> {
    let (len_a, len_b) = (1163, 4321);
    let dip = |index: usize, center: usize| 190.0 * (1.0 - index.abs_diff(center).to_f64_exact() / 3.0).max(0.0);
    let a = hard_operand(0.0, 0.026, len_a, |i, _| {
      let d = i.to_f64_exact() / len_a.to_f64_exact() - 0.264;
      let smooth = if d < 0.0 {
        700.0 * d * d
      } else {
        1050.0 * d * d + 8000.0 * (d - 0.5).max(0.0).powi(2)
      };
      smooth - dip(i, 1085)
    })?;
    let b = hard_operand(0.0, 0.026, len_b, |j, _| {
      let v = j.to_f64_exact() / len_b.to_f64_exact() - 0.2815;
      let smooth = if v < 0.0 {
        115.0 * (v / 0.2815).powi(2)
      } else {
        -139.0 * (-(v / 0.27).powi(2)).exp_m1()
      };
      smooth - dip(j, 4200)
    })?;
    let traced = traced_convolution(&a, &b)?;
    traced.assert_covers_the_whole_support();
    assert_eq!(Vec::<usize>::new(), traced.infinite_samples());
    assert_ne!(0, traced.windowed_fft_samples());
    traced.assert_matches_the_quadrature_oracle();
    Ok(())
  }

  #[test]
  fn test_convolution_with_wide_overlaps_at_a_hard_end_computes_its_tail_with_fewer_quadrature_terms_than_one_fft()
  -> Result<(), Report> {
    let (len_a, len_b) = (2642, 1476);
    let a = hard_operand(0.0, 0.01, len_a, |i, _| {
      let x = i.to_f64_exact() / (len_a - 1).to_f64_exact();
      match x {
        x if x <= 0.35 => 3700.0 * (x - 0.09).powi(2),
        x if x <= 0.55 => 250.12 + 140.0 * (x - 0.35),
        x if x <= 0.85 => 278.12 + 300.0 * (x - 0.55) + 4300.0 * (x - 0.55).powi(2),
        x if x <= 0.925 => 755.12 + 380.0 * (x - 0.85),
        x => 783.62 + 380.0 * (x - 0.925) + 8e5 * (x - 0.925).powi(3),
      }
    })?;
    let b = hard_operand(0.0, 0.01, len_b, |j, _| {
      let y = j.to_f64_exact() / (len_b - 1).to_f64_exact();
      if y <= 0.6295 {
        1150.0 * (y - 0.3165).powi(2)
      } else {
        112.66 + 48.0 * (y - 0.6295)
      }
    })?;
    let traced = traced_convolution(&a, &b)?;
    traced.assert_covers_the_whole_support();
    assert_eq!(Vec::<usize>::new(), traced.infinite_samples());
    assert_ne!(0, traced.windowed_fft_samples());
    traced.assert_matches_the_quadrature_oracle();
    let points = (len_a + len_b).to_f64_exact();
    let one_fft_of_the_operands = points * points.log2();
    assert!(traced.quadrature_terms().to_f64_exact() <= one_fft_of_the_operands);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::support_gap(  f64::INFINITY, (61..=159).collect_vec())]
  #[case::deep_valley(  1500.0,        vec![]                  )]
  #[trace]
  fn test_convolution_valley_quadrature_keeps_exact_zeros_and_depths_below_the_smallest_double(
    #[case] gap_neg_log: f64,
    #[case] expected_infinite: Vec<usize>,
  ) -> Result<(), Report> {
    let a = hard_operand(0.0, 0.05, 201, |_, t| match t {
      t if t <= 2.0 + 1e-9 => (t - 1.0).powi(2),
      t if t >= 8.0 - 1e-9 => (t - 9.0).powi(2),
      _ => gap_neg_log,
    })?;
    let b = hard_operand(0.0, 0.05, 21, |_, t| 4.0 * (t - 0.5).powi(2))?;
    let traced = traced_convolution(&a, &b)?;
    traced.assert_covers_the_whole_support();
    assert_eq!(expected_infinite, traced.infinite_samples());
    traced.assert_matches_the_quadrature_oracle();
    Ok(())
  }

  mod helpers {
    use crate::Distribution;
    use crate::distribution_core::function::DistributionFunction;
    use crate::distribution_ops::__tests__::test_edge_convolution::tests::helpers::ln_gauss_integral;
    use crate::distribution_ops::convolve::{SampleTrace, convolution_function_function_fine_traced};
    use crate::policy::NegLog;
    use eyre::Report;
    use itertools::Itertools;
    use ndarray::{Array1, Zip};
    use pretty_assertions::assert_eq;
    use std::collections::BTreeMap;
    use treetime_grid::{BoundaryBehavior, MaxGridPoints};
    use treetime_ops::{FftConvolver, FftPlans, UNIT_ROUND_OFF, rounding_gamma};
    use treetime_utils::array::ndarray::{max_or, min_or};
    use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};
    use treetime_utils::make_error;

    pub(super) fn hard_operand(
      x_min: f64,
      dx: f64,
      n: usize,
      neg_log: impl Fn(usize, f64) -> f64,
    ) -> Result<DistributionFunction<f64, NegLog>, Report> {
      let y: Array1<f64> = (0..n).map(|i| neg_log(i, x_min + i.to_f64_exact() * dx)).collect();
      DistributionFunction::from_start_dx_values(x_min, dx, y)?
        .with_left_extrap(BoundaryBehavior::Hard)?
        .with_right_extrap(BoundaryBehavior::Hard)
    }

    pub(super) fn traced_convolution(
      a: &DistributionFunction<f64, NegLog>,
      b: &DistributionFunction<f64, NegLog>,
    ) -> Result<TracedConvolution, Report> {
      let mut sources = SourceRecorder::default();
      let result =
        convolution_function_function_fine_traced(a, b, MaxGridPoints::default(), &FftPlans::new(), &mut sources)?;
      let result_x_start = a.x_min() + b.x_min();
      let Distribution::Function(result) = result else {
        return make_error!("the convolution of two grid functions must be a grid function, got {result}");
      };
      let (peak_a, peak_b) = (min_or(a.y(), f64::INFINITY), min_or(b.y(), f64::INFINITY));
      let (a, b, dx) = (a.y().mapv(|y| y - peak_a), b.y().mapv(|y| y - peak_b), a.dx());
      let round_offs = sources
        .sources
        .values()
        .filter_map(|source| match *source {
          SampleSource::Fft { theta } => Some((theta.to_bits(), theta)),
          SampleSource::Quadrature | SampleSource::WindowedFft => None,
        })
        .collect::<BTreeMap<_, _>>()
        .into_iter()
        .map(|(key, theta)| Ok((key, FftRoundOff::of(&a, &b, theta, dx)?)))
        .collect::<Result<_, Report>>()?;
      Ok(TracedConvolution {
        a,
        b,
        peaks: peak_a + peak_b,
        dx,
        x_start: result_x_start,
        result,
        sources: sources.sources,
        quadrature_terms: sources.quadrature_terms,
        round_offs,
      })
    }

    pub(super) struct TracedConvolution {
      a: Array1<f64>,
      b: Array1<f64>,
      peaks: f64,
      dx: f64,
      x_start: f64,
      result: DistributionFunction<f64, NegLog>,
      sources: BTreeMap<usize, SampleSource>,
      quadrature_terms: usize,
      round_offs: BTreeMap<u64, FftRoundOff>,
    }

    impl TracedConvolution {
      pub(super) fn result_len(&self) -> usize {
        self.a.len() + self.b.len() - 1
      }

      pub(super) fn samples(&self) -> impl Iterator<Item = (usize, f64)> + '_ {
        let first = f64_to_usize_saturating(((self.result.x_min() - self.x_start) / self.dx).round());
        self.result.y().iter().enumerate().map(move |(i, &y)| (first + i, y))
      }

      pub(super) fn assert_covers_the_whole_support(&self) {
        let covered = self.samples().map(|(k, _)| k).collect_vec();
        assert_eq!((0..self.result_len()).collect_vec(), covered);
      }

      pub(super) fn quadrature_terms(&self) -> usize {
        self.quadrature_terms
      }

      pub(super) fn windowed_fft_samples(&self) -> usize {
        self
          .sources
          .values()
          .filter(|source| matches!(source, SampleSource::WindowedFft))
          .count()
      }

      pub(super) fn infinite_samples(&self) -> Vec<usize> {
        self.samples().filter(|(_, y)| !y.is_finite()).map(|(k, _)| k).collect()
      }

      pub(super) fn assert_matches_the_quadrature_oracle(&self) {
        let failures = self
          .samples()
          .filter_map(|(k, y)| self.oracle_mismatch(k, y))
          .collect_vec();
        assert_eq!(Vec::<(usize, f64, f64, f64)>::new(), failures);
      }

      fn oracle_mismatch(&self, k: usize, y: f64) -> Option<(usize, f64, f64, f64)> {
        let Some(&source) = self.sources.get(&k) else {
          panic!("convolution sample {k} has no recorded source");
        };
        if let SampleSource::Fft { theta } = source {
          let fft = &self.round_offs[&theta.to_bits()];
          let tilt = theta * (k.to_f64_exact() * self.dx);
          let computed = (fft.min_a + fft.min_b + tilt - (y - self.peaks)).exp();
          if computed < 1e-7 * fft.peak {
            let expected = fft.trapezoid_sum(k);
            let tolerance = fft.absolute_bound(k, computed, [y, self.peaks, tilt]);
            return ((computed - expected).abs() > tolerance).then_some((k, computed, expected, tolerance));
          }
        }
        let oracle = self.oracle(k);
        let same_infinity =
          y.is_infinite() && oracle.is_infinite() && y.is_sign_positive() == oracle.is_sign_positive();
        (!(same_infinity || (y - oracle).abs() <= 1e-6)).then_some((k, y, oracle, 1e-6))
      }

      pub(super) fn oracle(&self, k: usize) -> f64 {
        let (last_a, last_b) = (self.a.len() - 1, self.b.len() - 1);
        let (lo, hi) = (k.saturating_sub(last_b), k.min(last_a));
        let costs = (lo..=hi)
          .map(|i| {
            let weight: f64 = match i {
              _ if lo == hi => 0.25,
              i if i == lo || i == hi => 0.5,
              _ => 1.0,
            };
            self.a[i] + self.b[k - i] - weight.ln()
          })
          .collect_vec();
        let min = costs.iter().copied().fold(f64::INFINITY, f64::min);
        if !min.is_finite() {
          return f64::INFINITY;
        }
        let sum: f64 = costs.iter().map(|cost| (min - cost).exp()).sum();
        self.peaks + min - (sum * self.dx).ln()
      }
    }

    #[derive(Default)]
    struct SourceRecorder {
      sources: BTreeMap<usize, SampleSource>,
      quadrature_terms: usize,
    }

    impl SampleTrace for SourceRecorder {
      fn fft_sample(&mut self, sample: usize, theta: f64) {
        self.sources.insert(sample, SampleSource::Fft { theta });
      }

      fn windowed_fft_sample(&mut self, sample: usize) {
        self.sources.insert(sample, SampleSource::WindowedFft);
      }

      fn quadrature_sample(&mut self, sample: usize, terms: usize) {
        self.sources.insert(sample, SampleSource::Quadrature);
        self.quadrature_terms += terms;
      }
    }

    #[derive(Clone, Copy, Debug)]
    enum SampleSource {
      Quadrature,
      WindowedFft,
      Fft { theta: f64 },
    }

    struct FftRoundOff {
      min_a: f64,
      min_b: f64,
      plain_a: Array1<f64>,
      plain_b: Array1<f64>,
      dx: f64,
      peak: f64,
      overlap_sums_bound: f64,
    }

    impl FftRoundOff {
      fn of(a: &Array1<f64>, b: &Array1<f64>, theta: f64, dx: f64) -> Result<Self, Report> {
        let (min_a, plain_a) = tilted_plain(a, theta, dx);
        let (min_b, plain_b) = tilted_plain(b, theta, dx);
        let peak = dx * max_or(&direct_overlap_sums(&plain_a, &plain_b), 0.0);
        let overlap_sums_bound =
          FftConvolver::new(plain_a.len(), plain_b.len(), &FftPlans::new())?.overlap_sums_round_off(&plain_a, &plain_b);
        Ok(Self {
          min_a,
          min_b,
          plain_a,
          plain_b,
          dx,
          peak,
          overlap_sums_bound,
        })
      }

      fn trapezoid_sum(&self, k: usize) -> f64 {
        let (lo, hi) = (k.saturating_sub(self.plain_b.len() - 1), k.min(self.plain_a.len() - 1));
        let products = (lo..=hi).map(|i| self.plain_a[i] * self.plain_b[k - i]);
        let sum: f64 = if lo == hi {
          0.25 * self.plain_a[lo] * self.plain_b[k - lo]
        } else {
          products.sum::<f64>()
            - f64::midpoint(
              self.plain_a[lo] * self.plain_b[k - lo],
              self.plain_a[hi] * self.plain_b[k - hi],
            )
        };
        sum * self.dx
      }

      fn absolute_bound(&self, k: usize, computed: f64, assembled_terms: [f64; 3]) -> f64 {
        let pairs = k.min(self.plain_a.len() - 1) + 1 - k.saturating_sub(self.plain_b.len() - 1);
        let exponents = [self.min_a, self.min_b].into_iter().chain(assembled_terms);
        let assembly = 8.0 * UNIT_ROUND_OFF * exponents.map(f64::abs).sum::<f64>();
        let fft = self.dx * self.overlap_sums_bound * (1.0 + UNIT_ROUND_OFF) + rounding_gamma(5) * self.peak;
        let underflow = pairs.to_f64_exact() * f64::MIN_POSITIVE * self.dx;
        fft + computed * assembly.exp_m1() + rounding_gamma(pairs + 4) * self.trapezoid_sum(k) + underflow
      }
    }

    fn tilted_plain(neg_log: &Array1<f64>, theta: f64, dx: f64) -> (f64, Array1<f64>) {
      let tilted = Zip::indexed(neg_log).map_collect(|i, &value| value - theta * (i.to_f64_exact() * dx));
      let min = min_or(&tilted, f64::INFINITY);
      (min, tilted.mapv(|value| (min - value).exp()))
    }

    fn direct_overlap_sums(f: &Array1<f64>, g: &Array1<f64>) -> Array1<f64> {
      (0..f.len() + g.len() - 1)
        .map(|n| {
          let (lo, hi) = (n.saturating_sub(g.len() - 1), n.min(f.len() - 1));
          (lo..=hi).map(|i| f[i] * g[n - i]).sum()
        })
        .collect()
    }

    pub(super) struct ClosedForm {
      pub(super) neg_log: f64,
      pub(super) discretization_bound: f64,
    }

    pub(super) fn mixture_closed_form(t: f64, dx: f64) -> ClosedForm {
      let (lo, hi) = ((t - 6.0).max(-6.0), (t + 6.0).min(31.0));
      let (log_integrals, log_errors): (Vec<f64>, Vec<f64>) = [(0.0, 0.0), (25.0, 200.0)]
        .into_iter()
        .map(|(mu, neg_log_weight)| {
          let middle = f64::midpoint(t, mu);
          let log_scale = -neg_log_weight - (t - mu).powi(2) / 4.0;
          let log_integral = log_scale + ln_gauss_integral(lo - middle, hi - middle);
          let log_error = log_scale + ((hi - lo) * dx * dx / 6.0).ln();
          (log_integral, log_error)
        })
        .unzip();
      let log_sum = |terms: &[f64]| {
        let max = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        max + terms.iter().map(|term| (term - max).exp()).sum::<f64>().ln()
      };
      let (log_integral, log_error) = (log_sum(&log_integrals), log_sum(&log_errors));
      let relative_error = (log_error - log_integral).exp();
      ClosedForm {
        neg_log: -log_integral,
        discretization_bound: -(-relative_error).ln_1p(),
      }
    }
  }
}
