use crate::seq::indel::{InDel, InDelKind};
use treetime_primitives::Seq;

pub(crate) fn insertion(range: (usize, usize), seq: impl Into<Seq>) -> InDel {
  let seq = seq.into();
  assert_eq!(
    range.1 - range.0,
    seq.len(),
    "insertion range must match its sequence length"
  );
  InDel {
    range,
    seq,
    kind: InDelKind::Insertion,
  }
}
