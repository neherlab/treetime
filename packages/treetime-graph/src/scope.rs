use crate::edge::GraphEdgeKey;
use crate::key_map::EdgeMap;
use crate::node::GraphNodeKey;
use crate::schedule::GraphSchedule;
use eyre::Report;
use itertools::Itertools;
use smallvec::SmallVec;
use treetime_utils::interval::range::range_contains;
use treetime_utils::interval::range_difference::range_difference;
use treetime_utils::interval::range_intersection::range_intersection;
use treetime_utils::interval::range_union::range_union;
use treetime_utils::{make_error, make_internal_error};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum SiteScope {
  #[default]
  Empty,
  Full,
  Ranges(SmallVec<[(usize, usize); 2]>),
}

impl SiteScope {
  pub fn from_ranges(ranges: impl IntoIterator<Item = (usize, usize)>) -> Self {
    let ranges = ranges.into_iter().filter(|(start, end)| start < end).collect_vec();
    Self::normalized(range_union(&[ranges]))
  }

  fn normalized(ranges: Vec<(usize, usize)>) -> Self {
    if ranges.is_empty() {
      Self::Empty
    } else {
      Self::Ranges(ranges.into())
    }
  }

  pub fn is_empty(&self) -> bool {
    matches!(self, Self::Empty)
  }

  pub fn is_full(&self) -> bool {
    matches!(self, Self::Full)
  }

  pub fn contains(&self, site: usize) -> bool {
    match self {
      Self::Empty => false,
      Self::Full => true,
      Self::Ranges(ranges) => range_contains(ranges, site),
    }
  }

  pub fn overlaps(&self, other: &Self) -> bool {
    match (self, other) {
      (Self::Empty, _) | (_, Self::Empty) => false,
      (Self::Full, _) | (_, Self::Full) => true,
      (Self::Ranges(a), Self::Ranges(b)) => !range_intersection(&[a.to_vec(), b.to_vec()]).is_empty(),
    }
  }

  #[must_use]
  pub fn union(&self, other: &Self) -> Self {
    match (self, other) {
      (Self::Full, _) | (_, Self::Full) => Self::Full,
      (Self::Empty, scope) | (scope, Self::Empty) => scope.clone(),
      (Self::Ranges(a), Self::Ranges(b)) => Self::normalized(range_union(&[a.to_vec(), b.to_vec()])),
    }
  }

  pub fn is_subset_of(&self, other: &Self) -> bool {
    match (self, other) {
      (Self::Empty, _) | (_, Self::Full) => true,
      (Self::Full, _) | (Self::Ranges(_), Self::Empty) => false,
      (Self::Ranges(a), Self::Ranges(b)) => range_difference(a, b).is_empty(),
    }
  }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EdgeScopes {
  segment_count: usize,
  by_edge: EdgeMap<SmallVec<[SiteScope; 4]>>,
}

impl EdgeScopes {
  pub fn new(segment_count: usize) -> Self {
    Self {
      segment_count,
      by_edge: EdgeMap::new(),
    }
  }

  pub fn segment_count(&self) -> usize {
    self.segment_count
  }

  pub fn get(&self, edge: GraphEdgeKey, segment: usize) -> &SiteScope {
    const EMPTY: &SiteScope = &SiteScope::Empty;
    self
      .by_edge
      .get(edge)
      .and_then(|scopes| scopes.get(segment))
      .unwrap_or(EMPTY)
  }

  pub fn carries(&self, edge: GraphEdgeKey, segment: usize) -> bool {
    !self.get(edge, segment).is_empty()
  }

  pub fn set(&mut self, edge: GraphEdgeKey, segment: usize, scope: SiteScope) -> Result<(), Report> {
    if segment >= self.segment_count {
      return make_internal_error!(
        "Segment {segment} is out of range: the edge scopes have {} segments",
        self.segment_count
      );
    }
    let count = self.segment_count;
    let scopes = self
      .by_edge
      .entry(edge)
      .or_insert_with(|| SmallVec::from_elem(SiteScope::Empty, count));
    scopes[segment] = scope;
    Ok(())
  }

  pub fn validate(&self, schedule: &GraphSchedule<'_>, name: impl Fn(GraphNodeKey) -> String) -> Result<(), Report> {
    for &node in schedule.preorder() {
      let parents = schedule.parents(node);
      let children = schedule.children(node);
      for segment in 0..self.segment_count {
        for (i, &(first, first_edge)) in parents.iter().enumerate() {
          for &(second, second_edge) in &parents[i + 1..] {
            if self.get(first_edge, segment).overlaps(self.get(second_edge, segment)) {
              return make_error!(
                "Node '{}' inherits sites of segment {} from two parents, '{}' and '{}', but every site has one parent",
                name(node),
                segment + 1,
                name(first),
                name(second)
              );
            }
          }
        }
        if children.is_empty() {
          continue;
        }
        let held = children
          .iter()
          .fold(SiteScope::Empty, |held, &(_, edge)| held.union(self.get(edge, segment)));
        if let Some(&(parent, _)) = parents
          .iter()
          .find(|&&(_, edge)| !self.get(edge, segment).is_subset_of(&held))
        {
          return make_error!(
            "The edge '{}' -> '{}' carries sites of segment {} that no child edge of '{}' carries",
            name(parent),
            name(node),
            segment + 1,
            name(node)
          );
        }
      }
    }
    Ok(())
  }
}
