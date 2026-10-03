use crate::seq::indel::{InDel, InDelKind};
use eyre::Report;
use treetime_primitives::Seq;
use treetime_utils::make_error;

pub(crate) fn insertion(range: (usize, usize), seq: impl Into<Seq>) -> InDel {
  indel(range, seq, InDelKind::Insertion).expect("insertion range must match its sequence length")
}

pub(crate) fn deletion(range: (usize, usize), seq: impl Into<Seq>) -> Result<InDel, Report> {
  indel(range, seq, InDelKind::Deletion)
}

fn indel(range: (usize, usize), seq: impl Into<Seq>, kind: InDelKind) -> Result<InDel, Report> {
  let seq = seq.into();
  if range.1.checked_sub(range.0) != Some(seq.len()) || seq.is_empty() {
    return make_error!(
      "Indel fixture range {}..{} must be non-empty and match its sequence length {}",
      range.0,
      range.1,
      seq.len()
    );
  }
  Ok(InDel { range, seq, kind })
}
