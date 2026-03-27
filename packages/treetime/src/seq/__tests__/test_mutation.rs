#[cfg(test)]
mod tests {
  use crate::seq::mutation::{Sub, compose_substitutions};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use treetime_primitives::AsciiChar;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
    Sub::new(c(reff), pos, c(qry)).unwrap()
  }

  #[test]
  fn test_mutation_compose_substitutions_both_empty() -> Result<(), Report> {
    let result = compose_substitutions(&[], &[])?;
    assert_eq!(result, vec![]);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_parent_empty() -> Result<(), Report> {
    let child = vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')];
    let result = compose_substitutions(&[], &child)?;
    assert_eq!(result, child);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_child_empty() -> Result<(), Report> {
    let parent = vec![sub(b'A', 0, b'T'), sub(b'G', 5, b'C')];
    let result = compose_substitutions(&parent, &[])?;
    assert_eq!(result, parent);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_non_overlapping() -> Result<(), Report> {
    let parent = vec![sub(b'A', 0, b'T'), sub(b'A', 1, b'C')];
    let child = vec![sub(b'G', 2, b'T'), sub(b'C', 3, b'A')];
    let result = compose_substitutions(&parent, &child)?;
    let expected = vec![
      sub(b'A', 0, b'T'),
      sub(b'A', 1, b'C'),
      sub(b'G', 2, b'T'),
      sub(b'C', 3, b'A'),
    ];
    assert_eq!(result, expected);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_chain() -> Result<(), Report> {
    // Parent: A→G at pos 0, child: G→T at pos 0 => net: A→T
    let parent = vec![sub(b'A', 0, b'G')];
    let child = vec![sub(b'G', 0, b'T')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(result, vec![sub(b'A', 0, b'T')]);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_cancellation() -> Result<(), Report> {
    // Parent: A→G at pos 0, child: G→A at pos 0 => cancel (back to original)
    let parent = vec![sub(b'A', 0, b'G')];
    let child = vec![sub(b'G', 0, b'A')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(result, vec![]);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_mixed() -> Result<(), Report> {
    // pos 0: parent only (passthrough)
    // pos 2: chain (A→G + G→T = A→T)
    // pos 4: child only (passthrough)
    // pos 6: cancellation (C→T + T→C = none)
    let parent = vec![
      sub(b'A', 0, b'T'),
      sub(b'A', 2, b'G'),
      sub(b'C', 6, b'T'),
    ];
    let child = vec![
      sub(b'G', 2, b'T'),
      sub(b'G', 4, b'C'),
      sub(b'T', 6, b'C'),
    ];
    let result = compose_substitutions(&parent, &child)?;
    let expected = vec![
      sub(b'A', 0, b'T'),
      sub(b'A', 2, b'T'),
      sub(b'G', 4, b'C'),
    ];
    assert_eq!(result, expected);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_output_sorted_by_position() -> Result<(), Report> {
    // Interleaved positions to verify merge order
    let parent = vec![sub(b'A', 1, b'T'), sub(b'G', 3, b'C'), sub(b'T', 5, b'A')];
    let child = vec![sub(b'C', 0, b'G'), sub(b'A', 2, b'T'), sub(b'C', 4, b'G')];
    let result = compose_substitutions(&parent, &child)?;
    for w in result.windows(2) {
      assert!(w[0].pos() < w[1].pos(), "output not sorted: pos {} >= {}", w[0].pos(), w[1].pos());
    }
    assert_eq!(result.len(), 6);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_all_cancel() -> Result<(), Report> {
    // Every position cancels
    let parent = vec![sub(b'A', 0, b'G'), sub(b'C', 1, b'T'), sub(b'G', 2, b'A')];
    let child = vec![sub(b'G', 0, b'A'), sub(b'T', 1, b'C'), sub(b'A', 2, b'G')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(result, vec![]);
    Ok(())
  }
}
