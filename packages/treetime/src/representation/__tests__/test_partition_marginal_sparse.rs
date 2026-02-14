#[cfg(test)]
mod tests {
  use eyre::Result;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  use crate::representation::partition::marginal_sparse::compose_substitutions;
  use crate::seq::mutation::Sub;

  fn sub(reff: char, pos: usize, qry: char) -> Sub {
    Sub::new(reff, pos, qry).unwrap()
  }

  #[test]
  fn test_compose_substitutions_empty_both() -> Result<()> {
    let result = compose_substitutions(&[], &[])?;
    assert_eq!(Vec::<Sub>::new(), result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_empty_parent() -> Result<()> {
    let child = vec![sub('A', 5, 'G'), sub('C', 10, 'T')];
    let result = compose_substitutions(&[], &child)?;
    assert_eq!(child, result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_empty_child() -> Result<()> {
    let parent = vec![sub('A', 5, 'G'), sub('C', 10, 'T')];
    let result = compose_substitutions(&parent, &[])?;
    assert_eq!(parent, result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_no_overlap() -> Result<()> {
    let parent = vec![sub('A', 5, 'G')];
    let child = vec![sub('C', 10, 'T')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(vec![sub('A', 5, 'G'), sub('C', 10, 'T')], result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_chain() -> Result<()> {
    // A5G + G5T = A5T (composition at same position)
    let parent = vec![sub('A', 5, 'G')];
    let child = vec![sub('G', 5, 'T')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(vec![sub('A', 5, 'T')], result);
    Ok(())
  }

  #[test]
  fn test_compose_substitutions_cancellation() -> Result<()> {
    // A5G + G5A = no mutation (cancellation)
    let parent = vec![sub('A', 5, 'G')];
    let child = vec![sub('G', 5, 'A')];
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
    let parent = vec![sub('A', 5, 'G'), sub('C', 10, 'T'), sub('T', 20, 'C')];
    let child = vec![sub('G', 5, 'T'), sub('G', 15, 'A'), sub('C', 20, 'T')];
    let result = compose_substitutions(&parent, &child)?;
    // Expected: pos 5 chained, pos 10 kept, pos 15 added, pos 20 cancelled
    assert_eq!(vec![sub('A', 5, 'T'), sub('C', 10, 'T'), sub('G', 15, 'A')], result);
    Ok(())
  }

  #[rstest]
  #[case('A', 'G', 'T', Some(('A', 'T')))] // chain: A->G->T = A->T
  #[case('A', 'G', 'A', None)] // cancel: A->G->A = none
  #[case('C', 'T', 'G', Some(('C', 'G')))] // chain: C->T->G = C->G
  #[case('G', 'A', 'G', None)] // cancel: G->A->G = none
  fn test_compose_substitutions_single_position(
    #[case] parent_reff: char,
    #[case] intermediate: char,
    #[case] child_qry: char,
    #[case] expected: Option<(char, char)>,
  ) -> Result<()> {
    let parent = vec![sub(parent_reff, 0, intermediate)];
    let child = vec![sub(intermediate, 0, child_qry)];
    let result = compose_substitutions(&parent, &child)?;
    let expected = expected.map_or_else(Vec::new, |(r, q)| vec![sub(r, 0, q)]);
    assert_eq!(expected, result);
    Ok(())
  }
}
