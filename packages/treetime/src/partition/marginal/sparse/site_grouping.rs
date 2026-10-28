use crate::alphabet::alphabet::Alphabet;
use crate::partition::marginal::sparse::candidate_sites::CandidateSites;
use crate::partition::storage::site_groups::SiteGroups;
use crate::partition::storage::sparse::INLINE_STATES;
use indexmap::IndexMap;
use itertools::Itertools;
use smallvec::SmallVec;
use treetime_primitives::{AsciiChar, AsciiTable, Seq};
use treetime_utils::interval::range_difference::range_difference;

pub(crate) fn group_sites(
  layers: &[Vec<LabeledRange>],
  excluded: &[(usize, usize)],
  sites: &CandidateSites,
  sequence: &Seq,
  alphabet: &Alphabet,
) -> Vec<GroupSite> {
  if layers.iter().all(Vec::is_empty) {
    return vec![];
  }
  let segments = segment_labels(layers, excluded);
  segments
    .into_iter()
    .filter_map(|(labels, ranges)| count_states(labels, &ranges, sites, sequence, alphabet))
    .collect()
}

pub(crate) fn message_layer(groups: &SiteGroups) -> Vec<LabeledRange> {
  groups
    .ranges()
    .map(|(start, end, group)| LabeledRange {
      start,
      end,
      label: Label::Group(group),
    })
    .collect()
}

pub(crate) fn missing_layer(missing: &[(usize, usize)], groups: &SiteGroups) -> Vec<LabeledRange> {
  let group_ranges = groups.ranges().map(|(start, end, _)| (start, end)).collect_vec();
  range_difference(missing, &group_ranges)
    .into_iter()
    .map(|(start, end)| LabeledRange {
      start,
      end,
      label: Label::Missing,
    })
    .merge_by(message_layer(groups), |a, b| a.start <= b.start)
    .collect()
}

pub(crate) fn positions_outside<'a>(
  ranges: &'a [(usize, usize)],
  sites: &'a CandidateSites,
) -> impl Iterator<Item = usize> + 'a {
  let candidates = sites.sites();
  let mut next = 0;
  ranges.iter().flat_map(|&(start, end)| start..end).filter(move |&pos| {
    while candidates.get(next).is_some_and(|&(site, _)| site < pos) {
      next += 1;
    }
    candidates.get(next).is_none_or(|&(site, _)| site != pos)
  })
}

pub(crate) fn ranges_contain(ranges: &[(usize, usize)], pos: usize) -> bool {
  let index = ranges.partition_point(|&(_, end)| end <= pos);
  ranges.get(index).is_some_and(|&(start, _)| start <= pos)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Label {
  Fixed,
  Missing,
  Group(usize),
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct LabeledRange {
  start: usize,
  end: usize,
  label: Label,
}

#[derive(Clone, Debug)]
pub(crate) struct GroupSite {
  pub(crate) labels: Labels,
  pub(crate) ranges: Vec<(usize, usize)>,
  pub(crate) counts: SmallVec<[(AsciiChar, usize); INLINE_STATES]>,
}

pub(crate) type Labels = SmallVec<[Label; 4]>;

fn segment_labels(layers: &[Vec<LabeledRange>], excluded: &[(usize, usize)]) -> IndexMap<Labels, Vec<(usize, usize)>> {
  let boundaries = layers
    .iter()
    .flatten()
    .flat_map(|range| [range.start, range.end])
    .chain(excluded.iter().flat_map(|&range| <[usize; 2]>::from(range)))
    .sorted_unstable()
    .dedup()
    .collect_vec();
  let mut cursors: SmallVec<[usize; 4]> = SmallVec::from_elem(0, layers.len());
  let mut segments: IndexMap<Labels, Vec<(usize, usize)>> = IndexMap::new();
  for (&start, &end) in boundaries.iter().tuple_windows() {
    if ranges_contain(excluded, start) {
      continue;
    }
    let labels: Labels = layers
      .iter()
      .zip(cursors.iter_mut())
      .map(|(layer, cursor)| {
        *cursor += layer[*cursor..].partition_point(|range| range.end <= start);
        layer
          .get(*cursor)
          .filter(|range| range.start <= start)
          .map_or(Label::Fixed, |range| range.label)
      })
      .collect();
    if labels.iter().all(|label| *label == Label::Fixed) {
      continue;
    }
    let ranges = segments.entry(labels).or_default();
    match ranges.last_mut() {
      Some(last) if last.1 == start => last.1 = end,
      _ => ranges.push((start, end)),
    }
  }
  segments
}

fn count_states(
  labels: Labels,
  ranges: &[(usize, usize)],
  sites: &CandidateSites,
  sequence: &Seq,
  alphabet: &Alphabet,
) -> Option<GroupSite> {
  let mut counts = AsciiTable::<usize>::default();
  let mut kept: Vec<(usize, usize)> = Vec::with_capacity(ranges.len());
  for pos in positions_outside(ranges, sites) {
    let state = sequence[pos];
    if !alphabet.is_canonical(state) {
      continue;
    }
    counts[state] += 1;
    match kept.last_mut() {
      Some(last) if last.1 == pos => last.1 = pos + 1,
      _ => kept.push((pos, pos + 1)),
    }
  }
  let counts: SmallVec<[(AsciiChar, usize); INLINE_STATES]> = alphabet
    .canonical()
    .map(|state| (state, counts[state]))
    .filter(|&(_, count)| count > 0)
    .collect();
  (!counts.is_empty()).then_some(GroupSite {
    labels,
    ranges: kept,
    counts,
  })
}
