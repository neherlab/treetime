#[cfg(test)]
mod tests {
  use crate::seq::mutation::{AlignedMutation, MutationEvent, Sub, compose_substitutions, mutation_event_strings};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use treetime_primitives::AsciiChar;
  use treetime_primitives::seq;

  #[test]
  fn test_mutation_compose_substitutions_both_empty() -> Result<(), Report> {
    let result = compose_substitutions(&[], &[])?;
    assert_eq!(result, vec![]);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_parent_empty() -> Result<(), Report> {
    let child = vec![helpers::sub(b'A', 0, b'T'), helpers::sub(b'G', 5, b'C')];
    let result = compose_substitutions(&[], &child)?;
    assert_eq!(result, child);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_child_empty() -> Result<(), Report> {
    let parent = vec![helpers::sub(b'A', 0, b'T'), helpers::sub(b'G', 5, b'C')];
    let result = compose_substitutions(&parent, &[])?;
    assert_eq!(result, parent);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_non_overlapping() -> Result<(), Report> {
    let parent = vec![helpers::sub(b'A', 0, b'T'), helpers::sub(b'A', 1, b'C')];
    let child = vec![helpers::sub(b'G', 2, b'T'), helpers::sub(b'C', 3, b'A')];
    let result = compose_substitutions(&parent, &child)?;
    let expected = vec![
      helpers::sub(b'A', 0, b'T'),
      helpers::sub(b'A', 1, b'C'),
      helpers::sub(b'G', 2, b'T'),
      helpers::sub(b'C', 3, b'A'),
    ];
    assert_eq!(result, expected);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_chain() -> Result<(), Report> {
    // Parent: A→G at pos 0, child: G→T at pos 0 => net: A→T
    let parent = vec![helpers::sub(b'A', 0, b'G')];
    let child = vec![helpers::sub(b'G', 0, b'T')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(result, vec![helpers::sub(b'A', 0, b'T')]);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_cancellation() -> Result<(), Report> {
    // Parent: A→G at pos 0, child: G→A at pos 0 => cancel (back to original)
    let parent = vec![helpers::sub(b'A', 0, b'G')];
    let child = vec![helpers::sub(b'G', 0, b'A')];
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
      helpers::sub(b'A', 0, b'T'),
      helpers::sub(b'A', 2, b'G'),
      helpers::sub(b'C', 6, b'T'),
    ];
    let child = vec![
      helpers::sub(b'G', 2, b'T'),
      helpers::sub(b'G', 4, b'C'),
      helpers::sub(b'T', 6, b'C'),
    ];
    let result = compose_substitutions(&parent, &child)?;
    let expected = vec![
      helpers::sub(b'A', 0, b'T'),
      helpers::sub(b'A', 2, b'T'),
      helpers::sub(b'G', 4, b'C'),
    ];
    assert_eq!(result, expected);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_output_sorted_by_position() -> Result<(), Report> {
    // Interleaved positions to verify merge order
    let parent = vec![
      helpers::sub(b'A', 1, b'T'),
      helpers::sub(b'G', 3, b'C'),
      helpers::sub(b'T', 5, b'A'),
    ];
    let child = vec![
      helpers::sub(b'C', 0, b'G'),
      helpers::sub(b'A', 2, b'T'),
      helpers::sub(b'C', 4, b'G'),
    ];
    let result = compose_substitutions(&parent, &child)?;
    let expected = vec![
      helpers::sub(b'C', 0, b'G'),
      helpers::sub(b'A', 1, b'T'),
      helpers::sub(b'A', 2, b'T'),
      helpers::sub(b'G', 3, b'C'),
      helpers::sub(b'C', 4, b'G'),
      helpers::sub(b'T', 5, b'A'),
    ];
    assert_eq!(result, expected);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_all_cancel() -> Result<(), Report> {
    // Every position cancels
    let parent = vec![
      helpers::sub(b'A', 0, b'G'),
      helpers::sub(b'C', 1, b'T'),
      helpers::sub(b'G', 2, b'A'),
    ];
    let child = vec![
      helpers::sub(b'G', 0, b'A'),
      helpers::sub(b'T', 1, b'C'),
      helpers::sub(b'A', 2, b'G'),
    ];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(result, vec![]);
    Ok(())
  }

  #[test]
  fn test_mutation_event_strings_expand_aligned_deletion() -> Result<(), Report> {
    // Oracle: TreeTime v0 spells gap transitions as one-based per-position substitutions.
    let event = MutationEvent::Deletion(AlignedMutation::new((1, 3), seq![helpers::c(b'C'), helpers::c(b'G')])?);
    let actual = mutation_event_strings(&event)?;
    let expected = vec!["C2-".to_owned(), "G3-".to_owned()];
    assert_eq!(expected, actual);
    Ok(())
  }

  proptest! {
    #[test]
    fn test_prop_mutation_event_strings_preserve_range(start in 0_usize..10_000, length in 1_usize..32) {
      let sequence = seq![helpers::c(b'A'); length];
      let event = MutationEvent::Insertion(AlignedMutation::new((start, start + length), sequence).unwrap());
      let actual = mutation_event_strings(&event).unwrap();
      let expected = (start + 1..=start + length)
        .map(|position| format!("-{position}A"))
        .collect::<Vec<_>>();
      prop_assert_eq!(expected, actual);
    }
  }

  mod helpers {
    use super::*;

    pub fn c(b: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(b)
    }

    pub fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
      Sub::new(c(reff), pos, c(qry)).unwrap()
    }
  }
}
