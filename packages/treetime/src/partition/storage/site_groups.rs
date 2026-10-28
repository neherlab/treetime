use crate::partition::storage::sparse::{FixedProfile, FixedProfiles, INLINE_STATES, SiteProfile};
use itertools::Itertools;
use smallvec::SmallVec;
use treetime_primitives::{AsciiChar, Seq};

#[derive(Clone, Debug, Default)]
pub struct SiteGroups {
  ranges: Vec<GroupRange>,
  profiles: Vec<FixedProfiles>,
}

impl SiteGroups {
  pub(crate) fn from_groups(groups: impl IntoIterator<Item = (Vec<(usize, usize)>, FixedProfiles)>) -> Self {
    let mut ranges = vec![];
    let mut profiles = vec![];
    for (group, (group_ranges, group_profiles)) in groups.into_iter().enumerate() {
      ranges.extend(
        group_ranges
          .into_iter()
          .map(|(start, end)| GroupRange { start, end, group }),
      );
      profiles.push(group_profiles);
    }
    ranges.sort_unstable_by_key(|range| range.start);
    Self { ranges, profiles }
  }

  pub(crate) fn is_empty(&self) -> bool {
    self.profiles.is_empty()
  }

  pub(crate) fn len(&self) -> usize {
    self.profiles.len()
  }

  pub(crate) fn profiles(&self, group: usize) -> &FixedProfiles {
    &self.profiles[group]
  }

  pub(crate) fn group_at(&self, pos: usize) -> Option<usize> {
    let index = self.ranges.partition_point(|range| range.end <= pos);
    self
      .ranges
      .get(index)
      .filter(|range| range.start <= pos)
      .map(|range| range.group)
  }

  pub(crate) fn entry(&self, pos: usize, state: AsciiChar) -> Option<&FixedProfile> {
    self.group_at(pos).and_then(|group| self.profiles[group].get(state))
  }

  pub(crate) fn first_position(&self, group: usize) -> Option<usize> {
    self
      .ranges
      .iter()
      .find(|range| range.group == group)
      .map(|range| range.start)
  }

  pub(crate) fn ranges(&self) -> impl Iterator<Item = (usize, usize, usize)> + '_ {
    self.ranges.iter().map(|range| (range.start, range.end, range.group))
  }

  pub(crate) fn map_profiles(&self, transform: impl Fn(&SiteProfile) -> SiteProfile) -> Self {
    Self {
      ranges: self.ranges.clone(),
      profiles: self
        .profiles
        .iter()
        .map(|profiles| profiles.map_dis(&transform))
        .collect(),
    }
  }

  pub(crate) fn entries<'a>(&'a self, sequence: &'a Seq) -> impl Iterator<Item = (usize, &'a FixedProfile)> + 'a {
    self.ranges.iter().flat_map(move |range| {
      (range.start..range.end)
        .filter_map(move |pos| self.profiles[range.group].get(sequence[pos]).map(|entry| (pos, entry)))
    })
  }

  pub(crate) fn entries_where<'a>(
    &'a self,
    sequence: &Seq,
    select: impl Fn(&FixedProfile) -> bool,
  ) -> Vec<(usize, &'a FixedProfile)> {
    let selected: Vec<SmallVec<[AsciiChar; INLINE_STATES]>> = self
      .profiles
      .iter()
      .map(|profiles| {
        profiles
          .iter()
          .filter(|entry| select(entry))
          .map(|entry| entry.state)
          .collect()
      })
      .collect_vec();
    self
      .ranges
      .iter()
      .filter(|range| !selected[range.group].is_empty())
      .flat_map(|range| {
        let states = &selected[range.group];
        let profiles = &self.profiles[range.group];
        (range.start..range.end).filter_map(move |pos| {
          let state = sequence[pos];
          profiles
            .get(state)
            .filter(|_| states.contains(&state))
            .map(|entry| (pos, entry))
        })
      })
      .collect()
  }
}

#[derive(Clone, Copy, Debug)]
struct GroupRange {
  start: usize,
  end: usize,
  group: usize,
}
