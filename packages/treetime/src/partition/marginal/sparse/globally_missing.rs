use crate::alphabet::alphabet::Alphabet;
use crate::partition::marginal::posterior_states::most_likely_state;
use crate::partition::marginal::sample::{Resolve, resolve_profile};
use crate::partition::storage::sparse::{SiteProfile, VarPos};
use itertools::{EitherOrBoth, Itertools};
use ndarray::{Array1, aview1};
use treetime_primitives::{AsciiChar, Seq};

#[derive(Clone, Debug, Default)]
pub(crate) struct GloballyMissing {
  ranges: Vec<(usize, usize)>,
  profile: SiteProfile,
  state: AsciiChar,
}

impl GloballyMissing {
  pub(crate) fn new(ranges: Vec<(usize, usize)>, pi: &Array1<f64>, alphabet: &Alphabet) -> Self {
    Self {
      ranges,
      profile: pi.iter().copied().collect(),
      state: alphabet.char(most_likely_state(pi).unwrap_or(0)),
    }
  }

  pub(crate) fn contains(&self, pos: usize) -> bool {
    let index = self.ranges.partition_point(|&(_, end)| end <= pos);
    self.ranges.get(index).is_some_and(|&(start, _)| start <= pos)
  }

  pub(crate) fn get(&self, pos: usize) -> Option<VarPos<'_>> {
    self.contains(pos).then(|| VarPos {
      dis: &self.profile,
      state: self.state,
    })
  }

  pub(crate) fn ranges(&self) -> &Vec<(usize, usize)> {
    &self.ranges
  }

  pub(crate) fn write_most_likely(&self, seq: &mut Seq, gap: AsciiChar) {
    for &(start, end) in &self.ranges {
      for state in seq[start..end].iter_mut().filter(|state| **state != gap) {
        *state = self.state;
      }
    }
  }

  pub(crate) fn most_likely_state_at(&self, pos: usize) -> Option<AsciiChar> {
    self.contains(pos).then_some(self.state)
  }

  pub(crate) fn profile(&self) -> &[f64] {
    &self.profile
  }

  pub(crate) fn resolve(&self, alphabet: &Alphabet, resolve: &mut Resolve<'_>) -> AsciiChar {
    match resolve {
      Resolve::Argmax => self.state,
      Resolve::Sample(_) => alphabet.char(resolve_profile(aview1(&self.profile), resolve)),
    }
  }

  pub(crate) fn merge_into<'a>(
    &'a self,
    explicit: impl Iterator<Item = (&'a usize, VarPos<'a>)> + 'a,
  ) -> impl Iterator<Item = (usize, SitePosterior<'a>)> + 'a {
    explicit
      .merge_join_by(self.positions(), |&(&pos, _), &missing| pos.cmp(&missing))
      .map(|entry| match entry {
        EitherOrBoth::Left((&pos, posterior)) => (pos, SitePosterior::Explicit(posterior)),
        EitherOrBoth::Right(pos) | EitherOrBoth::Both(_, pos) => (pos, SitePosterior::Equilibrium),
      })
  }

  fn positions(&self) -> impl Iterator<Item = usize> + '_ {
    self.ranges.iter().flat_map(|&(start, end)| start..end)
  }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum SitePosterior<'a> {
  Explicit(VarPos<'a>),
  Equilibrium,
}
