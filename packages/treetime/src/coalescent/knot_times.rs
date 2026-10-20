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
    let run = self.run_of(position);
    let next = if position + 1 < self.runs[run].end() {
      run
    } else {
      run + 1
    };
    (self.runs[run].at(position), self.runs[next].at(position + 1))
  }

  pub(crate) fn iter(&self) -> impl Iterator<Item = f64> + '_ {
    self.runs.iter().flat_map(KnotRun::iter)
  }

  pub(crate) fn partition_point(&self, pred: impl Fn(f64) -> bool) -> usize {
    self.partition_point_from(0, pred)
  }

  pub(crate) fn partition_point_from(&self, start: usize, pred: impl Fn(f64) -> bool) -> usize {
    self
      .runs
      .iter()
      .filter(|run| run.end() > start)
      .find(|run| !pred(run.at(run.end() - 1)))
      .map_or(self.len, |run| run.partition_point_from(start.max(run.first()), &pred))
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
    match self {
      Self::Listed { first, times } => start + times[start - first..].partition_point(|&time| pred(time)),
      Self::Grid {
        first,
        grid,
        start: index,
        count,
      } => {
        let (mut low, mut high) = (start, first + count);
        while low < high {
          let mid = usize::midpoint(low, high);
          if pred(grid.time(index + mid - first)) {
            low = mid + 1;
          } else {
            high = mid;
          }
        }
        low
      },
    }
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
