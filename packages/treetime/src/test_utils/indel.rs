use crate::seq::indel::{InDel, InDelKind};
use eyre::Report;
use treetime_primitives::Seq;
use treetime_utils::make_error;

pub(crate) fn insertion(range: (usize, usize), seq: impl Into<Seq>) -> InDel {
  indel(range, seq, InDelKind::Insertion).expect("insertion range must match its sequence length")
}

pub(crate) fn deletion(range: (usize, usize), seq: impl Into<Seq>) -> InDel {
  indel(range, seq, InDelKind::Deletion).expect("deletion range must match its sequence length")
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

#[cfg(test)]
mod tests {
  use super::indel;
  use crate::seq::indel::InDelKind;
  use rstest::rstest;
  use treetime_primitives::Seq;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::empty(   (2, 2), "",   "Indel fixture range 2..2 must be non-empty and match its sequence length 0")]
  #[case::short(   (2, 4), "A",  "Indel fixture range 2..4 must be non-empty and match its sequence length 1")]
  #[case::reversed((4, 2), "AC", "Indel fixture range 4..2 must be non-empty and match its sequence length 2")]
  #[trace]
  fn test_indel_fixture_rejects_range_that_does_not_match_sequence(
    #[case] range: (usize, usize),
    #[case] seq: &str,
    #[case] message: &str,
  ) {
    let seq = Seq::try_from_str(seq).unwrap();
    assert_error!(indel(range, seq, InDelKind::Deletion), message);
  }
}
