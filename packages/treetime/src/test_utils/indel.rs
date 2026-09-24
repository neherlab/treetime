use crate::seq::indel::{InDel, InDelKind};
use treetime_primitives::Seq;

pub(crate) fn insertion(range: (usize, usize), seq: impl Into<Seq>) -> InDel {
  InDel::new(range, seq, InDelKind::Insertion).expect("insertion range must match its sequence length")
}
