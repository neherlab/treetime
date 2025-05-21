use gcollections::ops::Bounded;
use interval::IntervalSet;
use interval::interval_set::ToIntervalSet;
use itertools::Itertools;

pub type Range = (usize, usize);
pub type RangeCollection = Vec<Range>;

pub fn range_contains(ranges: &[(usize, usize)], pos: usize) -> bool {
  range_contains_iter(ranges.iter(), pos)
}

pub fn range_contains_iter<'a>(mut ranges: impl Iterator<Item = &'a (usize, usize)> + 'a, pos: usize) -> bool {
  ranges.any(|(start, end)| *start <= pos && pos < *end)
}

pub fn to_interval_set(range_set: impl AsRef<[(usize, usize)]>) -> IntervalSet<usize> {
  range_set
    .as_ref()
    .iter()
    .map(|(start, end)| (*start, end.saturating_sub(1)))
    .collect_vec()
    .to_interval_set()
}

pub fn to_interval_sets<'a>(
  range_sets: impl Iterator<Item = &'a Vec<(usize, usize)>> + 'a,
) -> impl Iterator<Item = IntervalSet<usize>> + 'a {
  range_sets.map(to_interval_set)
}

pub fn from_interval_set(interval_set: IntervalSet<usize>) -> impl Iterator<Item = (usize, usize)> {
  interval_set
    .into_iter()
    .map(|interval| (interval.lower(), interval.upper().saturating_add(1)))
    .filter(|(start, end)| end > start)
}
