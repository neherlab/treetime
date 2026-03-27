#[cfg(test)]
mod tests {
  use eyre::Result;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_primitives::AsciiChar;

  use crate::seq::mutation::{Sub, compose_substitutions};

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
    Sub::new(c(reff), pos, c(qry)).unwrap()
  }

  #[test]
  fn test_compose_substitutions_empty_both() -> Result<()> {
    let result = compose_substitutions(&[], &[])?;
    assert_eq!(Vec::<Sub>::new(), result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_empty_parent() -> Result<()> {
    let child = vec![sub(b'A', 5, b'G'), sub(b'C', 10, b'T')];
    let result = compose_substitutions(&[], &child)?;
    assert_eq!(child, result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_empty_child() -> Result<()> {
    let parent = vec![sub(b'A', 5, b'G'), sub(b'C', 10, b'T')];
    let result = compose_substitutions(&parent, &[])?;
    assert_eq!(parent, result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_no_overlap() -> Result<()> {
    let parent = vec![sub(b'A', 5, b'G')];
    let child = vec![sub(b'C', 10, b'T')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(vec![sub(b'A', 5, b'G'), sub(b'C', 10, b'T')], result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_chain() -> Result<()> {
    // A5G + G5T = A5T (composition at same position)
    let parent = vec![sub(b'A', 5, b'G')];
    let child = vec![sub(b'G', 5, b'T')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(vec![sub(b'A', 5, b'T')], result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_cancellation() -> Result<()> {
    // A5G + G5A = no mutation (cancellation)
    let parent = vec![sub(b'A', 5, b'G')];
    let child = vec![sub(b'G', 5, b'A')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(Vec::<Sub>::new(), result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_mixed() -> Result<()> {
    // Position 5: A->G + G->T = A->T (chain)
    // Position 10: C->T (parent only)
    // Position 15: G->A (child only)
    // Position 20: T->C + C->T = no mutation (cancellation)
    let parent = vec![sub(b'A', 5, b'G'), sub(b'C', 10, b'T'), sub(b'T', 20, b'C')];
    let child = vec![sub(b'G', 5, b'T'), sub(b'G', 15, b'A'), sub(b'C', 20, b'T')];
    let result = compose_substitutions(&parent, &child)?;
    // Expected: pos 5 chained, pos 10 kept, pos 15 added, pos 20 cancelled
    assert_eq!(
      vec![sub(b'A', 5, b'T'), sub(b'C', 10, b'T'), sub(b'G', 15, b'A')],
      result
    );
    Ok(())
  }

  #[rstest]
  #[case(b'A', b'G', b'T', Some((b'A', b'T')))] // chain: A->G->T = A->T
  #[case(b'A', b'G', b'A', None)] // cancel: A->G->A = none
  #[case(b'C', b'T', b'G', Some((b'C', b'G')))] // chain: C->T->G = C->G
  #[case(b'G', b'A', b'G', None)] // cancel: G->A->G = none
  fn test_compose_substitutions_single_position(
    #[case] parent_reff: u8,
    #[case] intermediate: u8,
    #[case] child_qry: u8,
    #[case] expected: Option<(u8, u8)>,
  ) -> Result<()> {
    let parent = vec![sub(parent_reff, 0, intermediate)];
    let child = vec![sub(intermediate, 0, child_qry)];
    let result = compose_substitutions(&parent, &child)?;
    let expected = expected.map_or_else(Vec::new, |(r, q)| vec![sub(r, 0, q)]);
    assert_eq!(expected, result);
    Ok(())
  }
}
