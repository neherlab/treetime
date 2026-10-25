use treetime_primitives::{AsciiChar, AsciiTable};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct CharSlots {
  slots: AsciiTable<Option<usize>>,
  len: usize,
}

impl CharSlots {
  pub(crate) fn new(chars: impl IntoIterator<Item = AsciiChar>) -> Self {
    let mut slots = AsciiTable::<Option<usize>>::default();
    let mut len = 0;
    for c in chars {
      if slots[c].is_none() {
        slots[c] = Some(len);
        len += 1;
      }
    }
    Self { slots, len }
  }

  pub(crate) fn slot(&self, c: AsciiChar) -> Option<usize> {
    self.slots[c]
  }

  pub(crate) const fn len(&self) -> usize {
    self.len
  }
}
