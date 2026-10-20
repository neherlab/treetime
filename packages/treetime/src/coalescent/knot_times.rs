use itertools::Either;
use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};

#[derive(Clone, Debug)]
pub(crate) struct KnotTimes {
  runs: Vec<KnotRun>,
  len: usize,
}

impl KnotTimes {
  pub(crate) fn new(times: &[f64], grid: Option<KnotGrid>) -> Self {
    let mut runs: Vec<KnotRun> = Vec::new();
    for (position, &time) in times.iter().enumerate() {
      let index = grid.and_then(|grid| grid.index_of(time).map(|index| (grid, index)));
      match (runs.last_mut(), index) {
        (Some(KnotRun::Grid { start, count, .. }), Some((_, index))) if *start + *count == index => *count += 1,
        (_, Some((grid, index))) => runs.push(KnotRun::Grid {
          first: position,
          grid,
          start: index,
          count: 1,
        }),
        (Some(KnotRun::Listed { times, .. }), None) => times.push(time),
        (_, None) => runs.push(KnotRun::Listed {
          first: position,
          times: vec![time],
        }),
      }
    }
    for run in &mut runs {
      if let KnotRun::Listed { times, .. } = run {
        times.shrink_to_fit();
      }
    }
    runs.shrink_to_fit();
    Self { runs, len: times.len() }
  }

  pub(crate) fn listed(times: Vec<f64>) -> Self {
    let len = times.len();
    Self {
      runs: vec![KnotRun::Listed { first: 0, times }],
      len,
    }
  }

  pub(crate) fn len(&self) -> usize {
    self.len
  }

  pub(crate) fn first(&self) -> f64 {
    self.at(0)
  }

  pub(crate) fn last(&self) -> f64 {
    self.at(self.len - 1)
  }

  pub(crate) fn at(&self, position: usize) -> f64 {
    self.runs[self.run_of(position)].at(position)
  }

  pub(crate) fn segment(&self, position: usize) -> (f64, f64) {
    self.segment_in(self.run_of(position), position)
  }

  pub(crate) fn iter(&self) -> impl Iterator<Item = f64> + '_ {
    self.runs.iter().flat_map(KnotRun::iter)
  }

  pub(crate) fn partition_point(&self, pred: impl Fn(f64) -> bool) -> usize {
    self.seek(0, 0, pred).position
  }

  pub(crate) fn advance_while(&self, cursor: KnotCursor, pred: impl Fn(f64) -> bool) -> KnotCursor {
    let next = self.seek(cursor.run, cursor.position + 1, pred);
    let position = next.position - 1;
    let run = if self.runs.get(next.run).is_some_and(|run| run.first() <= position) {
      next.run
    } else {
      next.run - 1
    };
    KnotCursor { run, position }
  }

  pub(crate) fn segment_at(&self, cursor: KnotCursor) -> (f64, f64) {
    self.segment_in(cursor.run, cursor.position)
  }

  fn seek(&self, run: usize, start: usize, pred: impl Fn(f64) -> bool) -> KnotCursor {
    let run = gallop(run, self.runs.len(), |run| {
      pred(self.runs[run].at(self.runs[run].end() - 1))
    });
    let position = self.runs.get(run).map_or(self.len, |knots| {
      knots.partition_point_from(start.max(knots.first()), &pred)
    });
    KnotCursor { run, position }
  }

  fn segment_in(&self, run: usize, position: usize) -> (f64, f64) {
    let next = if position + 1 < self.runs[run].end() {
      run
    } else {
      run + 1
    };
    (self.runs[run].at(position), self.runs[next].at(position + 1))
  }

  fn run_of(&self, position: usize) -> usize {
    self.runs.partition_point(|run| run.first() <= position) - 1
  }
}

impl PartialEq for KnotTimes {
  fn eq(&self, other: &Self) -> bool {
    self.len == other.len && self.iter().eq(other.iter())
  }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct KnotCursor {
  run: usize,
  position: usize,
}

impl KnotCursor {
  pub(crate) fn position(self) -> usize {
    self.position
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct KnotGrid {
  pub origin: f64,
  pub step: f64,
}

impl KnotGrid {
  fn time(self, index: usize) -> f64 {
    self.origin + self.step * index.to_f64_exact()
  }

  fn index_of(self, time: f64) -> Option<usize> {
    let offset = ((time - self.origin) / self.step).round();
    if !(offset.is_finite() && offset >= 0.0) {
      return None;
    }
    let index = f64_to_usize_saturating(offset);
    (self.time(index).to_bits() == time.to_bits()).then_some(index)
  }
}

#[derive(Clone, Debug)]
enum KnotRun {
  Listed {
    first: usize,
    times: Vec<f64>,
  },
  Grid {
    first: usize,
    grid: KnotGrid,
    start: usize,
    count: usize,
  },
}

impl KnotRun {
  fn first(&self) -> usize {
    match self {
      Self::Listed { first, .. } | Self::Grid { first, .. } => *first,
    }
  }

  fn end(&self) -> usize {
    match self {
      Self::Listed { first, times } => first + times.len(),
      Self::Grid { first, count, .. } => first + count,
    }
  }

  fn partition_point_from(&self, start: usize, pred: impl Fn(f64) -> bool) -> usize {
    gallop(start, self.end(), |position| pred(self.at(position)))
  }

  fn at(&self, position: usize) -> f64 {
    match self {
      Self::Listed { first, times } => times[position - first],
      Self::Grid { first, grid, start, .. } => grid.time(start + position - first),
    }
  }

  fn iter(&self) -> impl Iterator<Item = f64> + '_ {
    match self {
      Self::Listed { times, .. } => Either::Left(times.iter().copied()),
      Self::Grid { grid, start, count, .. } => Either::Right((*start..start + count).map(|index| grid.time(index))),
    }
  }
}

fn gallop(low: usize, high: usize, holds: impl Fn(usize) -> bool) -> usize {
  let (mut low, mut high, mut step) = (low, high, 1);
  while low < high {
    let probe = (low + step - 1).min(high - 1);
    if !holds(probe) {
      high = probe;
      break;
    }
    low = probe + 1;
    step *= 2;
  }
  while low < high {
    let mid = usize::midpoint(low, high);
    if holds(mid) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }
  low
}
