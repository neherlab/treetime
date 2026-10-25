use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::char_slots::CharSlots;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use std::sync::Arc;
use treetime_primitives::{AsciiChar, AsciiTable};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Composition {
  counts: Box<[usize]>,
  slots: Arc<CharSlots>,
  gap: AsciiChar,
}

impl Composition {
  pub(crate) fn new(alphabet: &Alphabet) -> Self {
    let slots = Arc::clone(alphabet.char_slots());
    Self {
      counts: vec![0; slots.len()].into_boxed_slice(),
      slots,
      gap: alphabet.gap(),
    }
  }

  pub(crate) fn with_seq(sequence: impl AsRef<[AsciiChar]>, alphabet: &Alphabet) -> Self {
    let mut this = Self::new(alphabet);
    this.add_seq(sequence);
    this
  }

  pub(crate) fn get(&self, c: AsciiChar) -> Option<usize> {
    self.slots.slot(c).map(|slot| self.counts[slot])
  }

  pub(crate) fn add_seq(&mut self, sequence: impl AsRef<[AsciiChar]>) {
    for (c, &count) in char_histogram(sequence.as_ref())
      .iter()
      .filter(|(_, count)| **count > 0)
    {
      if let Some(slot) = self.slot(c) {
        self.counts[slot] += count;
      }
    }
  }

  pub(crate) fn add_sub(&mut self, sub: &Sub) {
    self.adjust_count(sub.reff(), -1);
    self.adjust_count(sub.qry(), 1);
  }

  pub(crate) fn add_indel(&mut self, indel: &InDel) {
    let adjust_by = if indel.is_deletion() { -1 } else { 1 };
    for nuc in &indel.seq {
      self.adjust_count(*nuc, adjust_by);
      self.adjust_count(self.gap, -adjust_by);
    }
  }

  pub(crate) fn adjust_count(&mut self, nuc: AsciiChar, change: isize) {
    if let Some(slot) = self.slot(nuc) {
      let count = &mut self.counts[slot];
      *count = count.saturating_add_signed(change);
    }
  }

  fn slot(&self, c: AsciiChar) -> Option<usize> {
    let slot = self.slots.slot(c);
    debug_assert!(
      slot.is_some(),
      "Composition counts only characters of its alphabet, got '{c}'"
    );
    slot
  }
}

fn char_histogram(sequence: &[AsciiChar]) -> AsciiTable<usize> {
  let mut lanes: [AsciiTable<usize>; 4] = Default::default();
  let (chunks, remainder) = sequence.as_chunks::<4>();
  for chunk in chunks {
    for (lane, &c) in lanes.iter_mut().zip(chunk) {
      lane[c] += 1;
    }
  }
  for &c in remainder {
    lanes[0][c] += 1;
  }
  AsciiTable::from_fn(|c| lanes.iter().map(|lane| lane[c]).sum())
}
