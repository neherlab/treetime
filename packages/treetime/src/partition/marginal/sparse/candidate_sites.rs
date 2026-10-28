use smallvec::{SmallVec, smallvec};
use treetime_primitives::AsciiChar;

#[derive(Clone, Debug)]
pub(crate) struct CandidateSiteEvents {
  events: Vec<SiteEvent>,
  n_messages: usize,
}

impl CandidateSiteEvents {
  pub(crate) fn with_capacity(n_messages: usize, capacity: usize) -> Self {
    Self {
      events: Vec::with_capacity(capacity),
      n_messages,
    }
  }

  pub(crate) fn insert(&mut self, pos: usize, slot: SiteSlot, value: AsciiChar) {
    self.push(pos, slot, value, true);
  }

  pub(crate) fn or_insert(&mut self, pos: usize, slot: SiteSlot, value: AsciiChar) {
    self.push(pos, slot, value, false);
  }

  pub(crate) fn finish(mut self) -> CandidateSites {
    self.events.sort_by_key(|event| event.pos);
    let n_positions = self.events.chunk_by(|a, b| a.pos == b.pos).count();
    let mut sites = CandidateSites {
      sites: Vec::with_capacity(n_positions),
      references: Vec::with_capacity(n_positions * self.n_messages),
      n_messages: self.n_messages,
    };
    let mut values: SmallVec<[Option<AsciiChar>; 4]> = smallvec![None; self.n_messages + 1];
    for group in self.events.chunk_by(|a, b| a.pos == b.pos) {
      values.fill(None);
      for event in group {
        let value = &mut values[event.slot];
        if event.overwrite || value.is_none() {
          *value = Some(event.value);
        }
      }
      if let Some(state) = values[0] {
        sites.sites.push((group[0].pos, state));
        sites.references.extend_from_slice(&values[1..]);
      }
    }
    sites
  }

  fn push(&mut self, pos: usize, slot: SiteSlot, value: AsciiChar, overwrite: bool) {
    let slot = slot.index();
    self.events.push(SiteEvent {
      pos,
      slot,
      overwrite,
      value,
    });
  }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum SiteSlot {
  State,
  Reference(usize),
}

impl SiteSlot {
  pub(crate) fn index(self) -> usize {
    match self {
      Self::State => 0,
      Self::Reference(index) => index + 1,
    }
  }
}

#[derive(Clone, Debug)]
pub(crate) struct CandidateSites {
  sites: Vec<(usize, AsciiChar)>,
  references: Vec<Option<AsciiChar>>,
  n_messages: usize,
}

impl CandidateSites {
  pub(crate) fn len(&self) -> usize {
    self.sites.len()
  }

  pub(crate) fn contains(&self, pos: usize) -> bool {
    self.sites.binary_search_by_key(&pos, |&(site, _)| site).is_ok()
  }

  pub(crate) fn iter(&self) -> impl Iterator<Item = (usize, AsciiChar, &[Option<AsciiChar>])> {
    let n = self.n_messages;
    self
      .sites
      .iter()
      .enumerate()
      .map(move |(index, &(pos, state))| (pos, state, &self.references[index * n..(index + 1) * n]))
  }

  pub(crate) fn relabel(&mut self, label: impl Fn(usize, AsciiChar) -> AsciiChar) {
    for (pos, state) in &mut self.sites {
      *state = label(*pos, *state);
    }
  }

  pub(crate) fn fill_missing_references(&mut self, fill: impl Fn(usize, usize, AsciiChar) -> AsciiChar) {
    let n = self.n_messages;
    for (index, &(pos, state)) in self.sites.iter().enumerate() {
      for (message, reference) in self.references[index * n..(index + 1) * n].iter_mut().enumerate() {
        if reference.is_none() {
          *reference = Some(fill(message, pos, state));
        }
      }
    }
  }
}

#[derive(Clone, Debug)]
struct SiteEvent {
  pos: usize,
  slot: usize,
  overwrite: bool,
  value: AsciiChar,
}
