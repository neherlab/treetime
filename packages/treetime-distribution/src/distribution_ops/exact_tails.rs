use crate::distribution_ops::convolve::{ResolvedNegLog, SampleTrace};
use eyre::Report;
use itertools::{Itertools, iproduct, izip};
use ndarray::{Array1, ArrayView1, Zip, s};
use ordered_float::OrderedFloat;
use std::collections::VecDeque;
use std::ops::RangeInclusive;
use treetime_ops::{FftConvolver, FftPlans, convolve_direct_neg_log_at, rounding_gamma, trapezoid_weight};
use treetime_utils::array::ndarray::min_or;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_internal_error, make_internal_report};

const NEGLIGIBLE_PAIR_DEPTH: f64 = 50.0;

const WINDOW_ROUND_OFF_NATS: f64 = 1e-6;

#[derive(Clone, Copy)]
pub(crate) struct ExactTailOperands<'a> {
  pub(crate) a: &'a Array1<f64>,
  pub(crate) b: &'a Array1<f64>,
  pub(crate) dx: f64,
  pub(crate) fft_plans: &'a FftPlans,
}

pub(crate) fn neg_log_on_whole_support(
  operands: ExactTailOperands<'_>,
  resolved: &ResolvedNegLog,
  target: IndexWindow,
  trace: &mut impl SampleTrace,
) -> Result<Option<ResolvedNegLog>, Report> {
  let mut tails = ExactTails::new(operands, resolved, target);
  tails.resolve_outside_support(trace)?;
  let (full_a, full_b) = (
    IndexWindow::of_len(operands.a.len()),
    IndexWindow::of_len(operands.b.len()),
  );
  let mut pending = tails
    .gaps_within(target)
    .into_iter()
    .map(|samples| Gap {
      samples,
      window_a: full_a,
      window_b: full_b,
      bound: f64::INFINITY,
    })
    .collect_vec();
  while let Some(gap) = pending.pop() {
    tails.resolve_gap(gap, &mut pending, trace)?;
  }
  tails.into_resolved()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct IndexWindow {
  pub(crate) first: usize,
  pub(crate) last: usize,
}

impl IndexWindow {
  pub(crate) fn of_len(len: usize) -> Self {
    Self {
      first: 0,
      last: len - 1,
    }
  }

  pub(crate) fn len(self) -> usize {
    self.last - self.first + 1
  }

  pub(crate) fn contains(self, i: usize) -> bool {
    (self.first..=self.last).contains(&i)
  }

  fn iter(self) -> RangeInclusive<usize> {
    self.first..=self.last
  }

  fn middle(self) -> usize {
    self.first.midpoint(self.last)
  }
}

struct ExactTails<'a> {
  operands: ExactTailOperands<'a>,
  first: usize,
  values: Vec<Option<f64>>,
}

impl<'a> ExactTails<'a> {
  fn new(operands: ExactTailOperands<'a>, resolved: &ResolvedNegLog, target: IndexWindow) -> Self {
    let mut values = vec![None; target.len()];
    resolved
      .neg_log
      .iter()
      .enumerate()
      .for_each(|(i, &value)| values[resolved.start + i - target.first] = Some(value));
    Self {
      operands,
      first: target.first,
      values,
    }
  }

  fn value(&self, k: usize) -> Option<f64> {
    self.values[k - self.first]
  }

  fn set(&mut self, k: usize, value: f64) {
    self.values[k - self.first] = Some(value);
  }

  fn target(&self) -> IndexWindow {
    IndexWindow {
      first: self.first,
      last: self.first + self.values.len() - 1,
    }
  }

  fn resolve_outside_support(&mut self, trace: &mut impl SampleTrace) -> Result<(), Report> {
    let support = result_support(self.operands.a, self.operands.b);
    let target = self.target();
    let in_support = |k: usize| {
      let i = support.partition_point(|interval| interval.last < k);
      support.get(i).is_some_and(|interval| interval.first <= k)
    };
    for k in target.iter() {
      if self.value(k).is_none() && !in_support(k) {
        self.set(k, f64::INFINITY);
        trace.quadrature_sample(k, 0);
      }
    }
    let full = Gap {
      samples: target,
      window_a: IndexWindow::of_len(self.operands.a.len()),
      window_b: IndexWindow::of_len(self.operands.b.len()),
      bound: f64::INFINITY,
    };
    for k in support.iter().flat_map(|interval| [interval.first, interval.last]) {
      if target.contains(k) && self.value(k).is_none() {
        self.resolve_by_quadrature(k, &full, trace)?;
      }
    }
    Ok(())
  }

  fn gaps_within(&self, range: IndexWindow) -> Vec<IndexWindow> {
    range
      .iter()
      .chunk_by(|&k| self.value(k).is_none())
      .into_iter()
      .filter(|(unresolved, _)| *unresolved)
      .filter_map(|(_, run)| run.minmax().into_option())
      .map(|(first, last)| IndexWindow { first, last })
      .collect()
  }

  fn resolve_gap(&mut self, gap: Gap, pending: &mut Vec<Gap>, trace: &mut impl SampleTrace) -> Result<(), Report> {
    let mut bounds = self.single_pair_bounds(&gap);
    let bound = gap.bound.min(bounds.iter().copied().fold(f64::NEG_INFINITY, f64::max));
    let Some((window_a, window_b)) = self.contributing_windows(&gap, bound) else {
      return make_internal_error!(
        "No operand pair contributes to the convolution samples {} to {} inside the support",
        gap.samples.first,
        gap.samples.last
      );
    };
    let gap = Gap {
      window_a,
      window_b,
      bound,
      ..gap
    };
    let sections = gap.sections();
    if self.quadrature_is_cheaper(&gap, &sections) {
      return gap
        .samples
        .iter()
        .try_for_each(|k| self.resolve_by_quadrature(k, &gap, trace));
    }
    let median_index = bounds.len().div_ceil(2) - 1;
    let median = *bounds.select_nth_unstable_by_key(median_index, |&b| OrderedFloat(b)).1;
    let resolved = if bound - median <= NEGLIGIBLE_PAIR_DEPTH {
      self.resolve_by_windowed_tilt(&gap, &sections, trace)?
    } else {
      0
    };
    if resolved == 0 {
      return self.split_at_middle(gap, pending, trace);
    }
    for samples in self.gaps_within(gap.samples) {
      let rest = Gap { samples, ..gap };
      if 2 * samples.len() > gap.samples.len() {
        self.split_at_middle(rest, pending, trace)?;
      } else {
        pending.push(rest);
      }
    }
    Ok(())
  }

  fn split_at_middle(&mut self, gap: Gap, pending: &mut Vec<Gap>, trace: &mut impl SampleTrace) -> Result<(), Report> {
    let middle = gap.samples.middle();
    self.resolve_by_quadrature(middle, &gap, trace)?;
    let halves = [
      (gap.samples.first < middle).then(|| IndexWindow {
        first: gap.samples.first,
        last: middle - 1,
      }),
      (middle < gap.samples.last).then(|| IndexWindow {
        first: middle + 1,
        last: gap.samples.last,
      }),
    ];
    pending.extend(halves.into_iter().flatten().map(|samples| Gap { samples, ..gap }));
    Ok(())
  }

  fn single_pair_bounds(&self, gap: &Gap) -> Vec<f64> {
    let samples = gap.samples;
    let probes = [samples.first, samples.middle(), samples.last]
      .into_iter()
      .dedup()
      .filter_map(|p| Some((p, self.best_pair(p, gap)?)))
      .collect_vec();
    samples
      .iter()
      .map(|k| {
        let full = self.full_pairs(k);
        probes
          .iter()
          .flat_map(|&(p, i)| [i, (i + k).saturating_sub(p)])
          .map(|i| self.pair_cost(k, i.clamp(full.first, full.last)))
          .fold(f64::INFINITY, f64::min)
      })
      .collect()
  }

  fn best_pair(&self, p: usize, gap: &Gap) -> Option<usize> {
    let (a, b) = (self.operands.a, self.operands.b);
    self
      .pairs_in_windows(p, gap)?
      .iter()
      .min_by_key(|&i| OrderedFloat(a[i] + b[p - i]))
  }

  fn pair_cost(&self, k: usize, i: usize) -> f64 {
    let (a, b) = (self.operands.a, self.operands.b);
    a[i] + b[k - i] - trapezoid_weight(a.len(), b.len(), k, i).ln()
  }

  fn full_pairs(&self, k: usize) -> IndexWindow {
    IndexWindow {
      first: k.saturating_sub(self.operands.b.len() - 1),
      last: k.min(self.operands.a.len() - 1),
    }
  }

  fn pairs_in_windows(&self, k: usize, gap: &Gap) -> Option<IndexWindow> {
    let full = self.full_pairs(k);
    let first = full
      .first
      .max(gap.window_a.first)
      .max(k.saturating_sub(gap.window_b.last));
    let last = full.last.min(gap.window_a.last).min(k.checked_sub(gap.window_b.first)?);
    (first <= last).then_some(IndexWindow { first, last })
  }

  fn contributing_windows(&self, gap: &Gap, bound: f64) -> Option<(IndexWindow, IndexWindow)> {
    let limit = bound + NEGLIGIBLE_PAIR_DEPTH;
    let (a, b) = (self.operands.a, self.operands.b);
    let window_a = contributing_window(a, b, gap.samples, gap.window_a, gap.window_b, limit)?;
    let window_b = contributing_window(b, a, gap.samples, gap.window_b, window_a, limit)?;
    Some((window_a, window_b))
  }

  fn quadrature_is_cheaper(&self, gap: &Gap, sections: &[Section]) -> bool {
    if sections
      .iter()
      .any(|section| section.a.len() < 2 || section.b.len() < 2)
    {
      return true;
    }
    let terms: usize = gap
      .samples
      .iter()
      .filter_map(|k| self.pairs_in_windows(k, gap))
      .map(IndexWindow::len)
      .sum();
    let fft_cost: f64 = sections
      .iter()
      .map(|section| {
        let points = (section.a.len() + section.b.len()).to_f64_exact();
        points * points.log2()
      })
      .sum();
    terms.to_f64_exact() <= fft_cost
  }

  fn resolve_by_quadrature(&mut self, k: usize, gap: &Gap, trace: &mut impl SampleTrace) -> Result<(), Report> {
    let Some(pairs) = self.pairs_in_windows(k, gap) else {
      return make_internal_error!("Convolution sample {k} inside the support has no operand pair in its windows");
    };
    let value = convolve_direct_neg_log_at(
      self.operands.dx,
      self.operands.a,
      self.operands.b,
      k,
      pairs.first..=pairs.last,
    )?;
    self.set(k, value);
    trace.quadrature_sample(k, pairs.len());
    Ok(())
  }

  fn resolve_by_windowed_tilt(
    &mut self,
    gap: &Gap,
    sections: &[Section],
    trace: &mut impl SampleTrace,
  ) -> Result<usize, Report> {
    let ExactTailOperands { a, b, dx, fft_plans } = self.operands;
    let samples = gap.samples;
    let (Some(left), Some(right)) = (self.value(samples.first - 1), self.value(samples.last + 1)) else {
      return make_internal_error!(
        "Convolution gap {} to {} is not bracketed by resolved samples",
        samples.first,
        samples.last
      );
    };
    if !(left.is_finite() && right.is_finite()) {
      return Ok(0);
    }
    let theta = (right - left) / ((samples.len() + 1).to_f64_exact() * dx);
    let tilted = sections
      .iter()
      .map(|section| {
        (
          TiltedWindow::new(a, section.a, theta, dx),
          TiltedWindow::new(b, section.b, theta, dx),
        )
      })
      .collect_vec();
    let reference = tilted
      .iter()
      .map(|(tilted_a, tilted_b)| tilted_a.min + tilted_b.min)
      .fold(f64::INFINITY, f64::min);
    let mut sums = Array1::<f64>::zeros(samples.len());
    let mut round_off = Array1::<f64>::zeros(samples.len());
    for (tilted_a, tilted_b) in &tilted {
      let offset = tilted_a.window.first + tilted_b.window.first;
      let outputs = IndexWindow {
        first: offset.max(samples.first),
        last: (tilted_a.window.last + tilted_b.window.last).min(samples.last),
      };
      if outputs.first > outputs.last || !(tilted_a.min + tilted_b.min).is_finite() {
        continue;
      }
      let gap_outputs = s![outputs.first - samples.first..=outputs.last - samples.first];
      let scale = (reference - (tilted_a.min + tilted_b.min)).exp();
      if scale < f64::MIN_POSITIVE {
        let largest_sum = tilted_a.plain.len().min(tilted_b.plain.len()).to_f64_exact();
        round_off
          .slice_mut(gap_outputs)
          .mapv_inplace(|bound| bound + f64::MIN_POSITIVE * largest_sum);
        continue;
      }
      let mut convolver = FftConvolver::new(tilted_a.plain.len(), tilted_b.plain.len(), fft_plans)?;
      let section_sums = convolver.overlap_sums(&tilted_a.plain, &tilted_b.plain)?;
      let section_round_off = scale * convolver.overlap_sums_round_off(&tilted_a.plain, &tilted_b.plain);
      Zip::from(sums.slice_mut(gap_outputs))
        .and(round_off.slice_mut(gap_outputs))
        .and(section_sums.slice(s![outputs.first - offset..=outputs.last - offset]))
        .for_each(|sum, bound, &section_sum| {
          *sum += scale * section_sum;
          *bound += section_round_off;
        });
    }
    let summation = rounding_gamma(sections.len() + 3);
    let evaluation = 1.0 + rounding_gamma(sections.len() + 8);
    let accepted_share = -(-WINDOW_ROUND_OFF_NATS).exp_m1();
    let tilted_at = |neg_log: &Array1<f64>, i: usize| neg_log[i] - theta * (i.to_f64_exact() * dx);
    let mut resolved = 0;
    for (k, &gross, &bound) in izip!(samples.iter(), &sums, &round_off) {
      let full = self.full_pairs(k);
      let excess: f64 = [full.first, full.last]
        .into_iter()
        .dedup()
        .filter(|&i| gap.window_a.contains(i) && gap.window_b.contains(k - i))
        .map(|i| {
          let product = (reference - (tilted_at(a, i) + tilted_at(b, k - i))).exp();
          (1.0 - trapezoid_weight(a.len(), b.len(), k, i)) * product
        })
        .sum();
      let sum = gross - excess;
      let error = (1.0 + 3.0 * summation) * bound + summation * (gross + excess);
      if sum > 0.0 && error * evaluation <= sum * accepted_share {
        self.set(k, -(sum * dx).ln() + reference + theta * (k.to_f64_exact() * dx));
        trace.windowed_fft_sample(k);
        resolved += 1;
      }
    }
    Ok(resolved)
  }

  fn into_resolved(self) -> Result<Option<ResolvedNegLog>, Report> {
    let values: Vec<f64> = self
      .values
      .iter()
      .map(|value| value.ok_or_else(|| make_internal_report!("Convolution sample left unresolved")))
      .try_collect()?;
    let (Some(first), Some(last)) = (
      values.iter().position(|v| v.is_finite()),
      values.iter().rposition(|v| v.is_finite()),
    ) else {
      return Ok(None);
    };
    Ok(Some(ResolvedNegLog {
      start: self.first + first,
      neg_log: values[first..=last].iter().copied().collect(),
    }))
  }
}

#[derive(Clone, Copy)]
struct Gap {
  samples: IndexWindow,
  window_a: IndexWindow,
  window_b: IndexWindow,
  bound: f64,
}

impl Gap {
  #[expect(
    clippy::integer_division,
    reason = "integer division splits the window into near-equal parts"
  )]
  fn sections(&self) -> Vec<Section> {
    let split_a = self.window_a.len() >= self.window_b.len();
    let (long, short) = if split_a {
      (self.window_a, self.window_b)
    } else {
      (self.window_b, self.window_a)
    };
    let parts = (long.len() / self.samples.len().max(2)).max(1);
    (0..parts)
      .filter_map(|p| {
        let part = IndexWindow {
          first: long.first + p * long.len() / parts,
          last: long.first + (p + 1) * long.len() / parts - 1,
        };
        let partners = IndexWindow {
          first: short.first.max(self.samples.first.saturating_sub(part.last)),
          last: short.last.min(self.samples.last.checked_sub(part.first)?),
        };
        (partners.first <= partners.last).then_some(if split_a {
          Section { a: part, b: partners }
        } else {
          Section { a: partners, b: part }
        })
      })
      .collect()
  }
}

#[derive(Clone, Copy)]
struct Section {
  a: IndexWindow,
  b: IndexWindow,
}

fn result_support(a: &Array1<f64>, b: &Array1<f64>) -> Vec<IndexWindow> {
  let (runs_a, runs_b) = (finite_runs(a), finite_runs(b));
  iproduct!(&runs_a, &runs_b)
    .map(|(run_a, run_b)| IndexWindow {
      first: run_a.first + run_b.first,
      last: run_a.last + run_b.last,
    })
    .sorted_by_key(|interval| interval.first)
    .coalesce(|x, y| {
      if y.first <= x.last + 1 {
        Ok(IndexWindow {
          first: x.first,
          last: x.last.max(y.last),
        })
      } else {
        Err((x, y))
      }
    })
    .collect()
}

fn finite_runs(neg_log: &Array1<f64>) -> Vec<IndexWindow> {
  neg_log
    .iter()
    .enumerate()
    .chunk_by(|(_, value)| value.is_finite())
    .into_iter()
    .filter(|(finite, _)| *finite)
    .filter_map(|(_, run)| run.map(|(i, _)| i).minmax().into_option())
    .map(|(first, last)| IndexWindow { first, last })
    .collect()
}

fn contributing_window(
  f: &Array1<f64>,
  g: &Array1<f64>,
  samples: IndexWindow,
  window_f: IndexWindow,
  window_g: IndexWindow,
  limit: f64,
) -> Option<IndexWindow> {
  let first = window_f.first.max(samples.first.saturating_sub(window_g.last));
  let last = window_f.last.min(samples.last.checked_sub(window_g.first)?);
  if first > last {
    return None;
  }
  let mut minima = SlidingMinimum::new(g);
  (first..=last)
    .rev()
    .filter(|&i| {
      let partners = IndexWindow {
        first: window_g.first.max(samples.first.saturating_sub(i)),
        last: window_g.last.min(samples.last - i),
      };
      f[i].is_finite() && f[i] + minima.over(partners) <= limit
    })
    .minmax()
    .into_option()
    .map(|(first, last)| IndexWindow { first, last })
}

struct SlidingMinimum<'a> {
  values: ArrayView1<'a, f64>,
  candidates: VecDeque<usize>,
  next: usize,
}

impl<'a> SlidingMinimum<'a> {
  fn new(values: &'a Array1<f64>) -> Self {
    Self {
      values: values.view(),
      candidates: VecDeque::new(),
      next: 0,
    }
  }

  fn over(&mut self, window: IndexWindow) -> f64 {
    self.next = self.next.max(window.first);
    while self.next <= window.last {
      let value = self.values[self.next];
      while self.candidates.back().is_some_and(|&j| self.values[j] >= value) {
        self.candidates.pop_back();
      }
      self.candidates.push_back(self.next);
      self.next += 1;
    }
    while self.candidates.front().is_some_and(|&j| j < window.first) {
      self.candidates.pop_front();
    }
    self.candidates.front().map_or(f64::INFINITY, |&j| self.values[j])
  }
}

struct TiltedWindow {
  window: IndexWindow,
  min: f64,
  plain: Array1<f64>,
}

impl TiltedWindow {
  fn new(neg_log: &Array1<f64>, window: IndexWindow, theta: f64, dx: f64) -> Self {
    let tilted = Zip::indexed(neg_log.slice(s![window.first..=window.last]))
      .map_collect(|i, &value| value - theta * ((window.first + i).to_f64_exact() * dx));
    let min = min_or(&tilted, f64::INFINITY);
    let plain = tilted.mapv(|value| (min - value).exp());
    Self { window, min, plain }
  }
}
