use crate::alphabet::alphabet::Alphabet;
use crate::partition::marginal::sparse::candidate_sites::CandidateSites;
use crate::partition::storage::site_groups::SiteGroups;
use crate::partition::storage::sparse::INLINE_STATES;
use indexmap::IndexMap;
use itertools::Itertools;
use smallvec::SmallVec;
use treetime_primitives::{AsciiChar, AsciiTable, Seq};
use treetime_utils::interval::range::range_contains_sorted;
use treetime_utils::interval::range_difference::range_difference;

pub(crate) fn group_sites(
  layers: &[Vec<SourceRange>],
  excluded: &[(usize, usize)],
  sites: &CandidateSites,
  sequence: &Seq,
  alphabet: &Alphabet,
) -> Vec<GroupedSites> {
  if layers.iter().all(Vec::is_empty) {
    return vec![];
  }
  let segments = segment_sources(layers, excluded);
  segments
    .into_iter()
    .filter_map(|(sources, ranges)| count_states(sources, &ranges, sites, sequence, alphabet))
    .collect()
}

pub(crate) fn group_layer(groups: &SiteGroups) -> Vec<SourceRange> {
  groups
    .ranges()
    .map(|(start, end, group)| SourceRange {
      start,
      end,
      source: ProfileSource::Group(group),
    })
    .collect()
}

pub(crate) fn missing_layer(missing: &[(usize, usize)], groups: &SiteGroups) -> Vec<SourceRange> {
  let group_ranges = groups.ranges().map(|(start, end, _)| (start, end)).collect_vec();
  range_difference(missing, &group_ranges)
    .into_iter()
    .map(|(start, end)| SourceRange {
      start,
      end,
      source: ProfileSource::Missing,
    })
    .merge_by(group_layer(groups), |a, b| a.start <= b.start)
    .collect()
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum ProfileSource {
  Fixed,
  Missing,
  Group(usize),
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct SourceRange {
  start: usize,
  end: usize,
  source: ProfileSource,
}

#[derive(Clone, Debug)]
pub(crate) struct GroupedSites {
  pub(crate) sources: ProfileSources,
  pub(crate) ranges: Vec<(usize, usize)>,
  pub(crate) counts: SmallVec<[(AsciiChar, usize); INLINE_STATES]>,
}

pub(crate) type ProfileSources = SmallVec<[ProfileSource; 4]>;

fn segment_sources(
  layers: &[Vec<SourceRange>],
  excluded: &[(usize, usize)],
) -> IndexMap<ProfileSources, Vec<(usize, usize)>> {
  let boundaries = layers
    .iter()
    .flatten()
    .flat_map(|range| [range.start, range.end])
    .chain(excluded.iter().flat_map(|&range| <[usize; 2]>::from(range)))
    .sorted_unstable()
    .dedup()
    .collect_vec();
  let mut cursors: SmallVec<[usize; 4]> = SmallVec::from_elem(0, layers.len());
  let mut segments: IndexMap<ProfileSources, Vec<(usize, usize)>> = IndexMap::new();
  for (&start, &end) in boundaries.iter().tuple_windows() {
    if range_contains_sorted(excluded, start) {
      continue;
    }
    let sources: ProfileSources = layers
      .iter()
      .zip(cursors.iter_mut())
      .map(|(layer, cursor)| {
        *cursor += layer[*cursor..].partition_point(|range| range.end <= start);
        layer
          .get(*cursor)
          .filter(|range| range.start <= start)
          .map_or(ProfileSource::Fixed, |range| range.source)
      })
      .collect();
    if sources.iter().all(|source| *source == ProfileSource::Fixed) {
      continue;
    }
    let ranges = segments.entry(sources).or_default();
    match ranges.last_mut() {
      Some(last) if last.1 == start => last.1 = end,
      _ => ranges.push((start, end)),
    }
  }
  segments
}

fn count_states(
  sources: ProfileSources,
  ranges: &[(usize, usize)],
  sites: &CandidateSites,
  sequence: &Seq,
  alphabet: &Alphabet,
) -> Option<GroupedSites> {
  let mut counts = AsciiTable::<usize>::default();
  let mut kept: Vec<(usize, usize)> = Vec::with_capacity(ranges.len());
  for pos in sites.non_candidate_positions(ranges) {
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
  (!counts.is_empty()).then_some(GroupedSites {
    sources,
    ranges: kept,
    counts,
  })
}
