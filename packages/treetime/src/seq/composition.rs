use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use std::collections::BTreeMap;
use treetime_primitives::AsciiChar;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Composition {
  counts: BTreeMap<AsciiChar, usize>,
  gap: AsciiChar,
}

impl Composition {
  pub(crate) fn new<I>(alphabet_chars: I, gap: AsciiChar) -> Self
  where
    I: IntoIterator<Item = AsciiChar>,
  {
    let counts = alphabet_chars.into_iter().map(|c| (c, 0)).collect();
    Self { counts, gap }
  }

  pub(crate) fn get(&self, c: AsciiChar) -> Option<usize> {
    self.counts.get(&c).copied()
  }

  pub(crate) fn counts(&self) -> &BTreeMap<AsciiChar, usize> {
    &self.counts
  }

  pub(crate) fn with_seq(
    sequence: impl AsRef<[AsciiChar]>,
    alphabet_chars: impl IntoIterator<Item = AsciiChar>,
    gap: AsciiChar,
  ) -> Self {
    let mut this = Self::new(alphabet_chars, gap);
    this.add_seq(sequence);
    this
  }

  pub(crate) fn add_seq(&mut self, sequence: impl AsRef<[AsciiChar]>) {
    self.add_histogram(&char_histogram(sequence.as_ref()));
  }

  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  fn add_histogram(&mut self, additions: &[usize; 128]) {
    for (index, &count) in additions.iter().enumerate().filter(|(_, count)| **count > 0) {
      let c = AsciiChar::from_byte_unchecked(index as u8);
      *self.counts.entry(c).or_default() += count;
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
    let count = self.counts.entry(nuc).or_default();
    *count = count.saturating_add_signed(change);
  }
}

fn char_histogram(sequence: &[AsciiChar]) -> [usize; 128] {
  let mut lanes = [[0_usize; 128]; 4];
  let (chunks, remainder) = sequence.as_chunks::<4>();
  for chunk in chunks {
    for (lane, &c) in lanes.iter_mut().zip(chunk) {
      lane[usize::from(c)] += 1;
    }
  }
  for &c in remainder {
    lanes[0][usize::from(c)] += 1;
  }
  std::array::from_fn(|index| lanes.iter().map(|lane| lane[index]).sum())
}
