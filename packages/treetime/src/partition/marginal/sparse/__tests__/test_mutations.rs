#[cfg(test)]
mod tests {
  use self::helpers::seq_with_changes;
  use crate::partition::marginal::sparse::mutations::differing_positions;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::sync::Arc;

  #[rustfmt::skip]
  #[rstest]
  #[case::empty(              0, vec![])]
  #[case::identical(        130, vec![])]
  #[case::single_block(      20, vec![0, 7, 19])]
  #[case::block_edges(      130, vec![0, 63, 64, 127, 128, 129])]
  #[case::exact_blocks(     128, vec![62, 63, 64, 65, 127])]
  #[case::every_position(    70, (0..70).collect_vec())]
  #[trace]
  fn test_mutations_differing_positions(#[case] len: usize, #[case] changes: Vec<usize>) {
    let parent = Arc::new(seq_with_changes(len, &[]));
    let child = Arc::new(seq_with_changes(len, &changes));

    let actual = differing_positions(&parent, &child).collect_vec();

    assert_eq!(changes, actual);
  }

  mod helpers {
    use treetime_primitives::Seq;

    pub(crate) fn seq_with_changes(len: usize, changes: &[usize]) -> Seq {
      let text: String = (0..len)
        .map(|pos| if changes.contains(&pos) { 'C' } else { 'A' })
        .collect();
      Seq::try_from_str(&text).expect("an ACGT text is a valid sequence")
    }
  }
}
