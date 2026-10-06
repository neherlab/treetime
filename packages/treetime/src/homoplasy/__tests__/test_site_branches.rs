#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::homoplasy::classify::MutationClass;
  use crate::homoplasy::site_branches::{SiteBranches, sites_by_branch_count};
  use crate::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub};
  use crate::test_utils::deletion;
  use eyre::Report;
  use helpers::events;
  use pretty_assertions::assert_eq;
  use treetime_primitives::Seq;

  #[test]
  fn test_sites_by_branch_count_ranks_by_branches_then_position() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let branches = vec![
      events(&["C20T", "A10G"])?,
      events(&["C20A"])?,
      events(&["A10G", "G30A"])?,
      events(&["T5C"])?,
    ];
    let actual = sites_by_branch_count(&branches, &alphabet, MutationClass::Substitution);
    let expected = vec![
      SiteBranches {
        position: 9,
        branches: 2,
      },
      SiteBranches {
        position: 19,
        branches: 2,
      },
      SiteBranches {
        position: 4,
        branches: 1,
      },
      SiteBranches {
        position: 29,
        branches: 1,
      },
    ];
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_sites_by_branch_count_keeps_only_the_requested_class() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let deleted = Mutation::indel(
      MutationTrack::Nucleotide,
      &deletion((9, 11), Seq::try_from_slice(b"AC")?),
    )?;
    let branches = vec![
      [events(&["G10R", "A6N"])?, vec![deleted.event]].concat(),
      events(&["G10R", "T7C"])?,
    ];
    let actual = (
      sites_by_branch_count(&branches, &alphabet, MutationClass::Substitution),
      sites_by_branch_count(&branches, &alphabet, MutationClass::Ambiguous),
    );
    let expected = (
      vec![SiteBranches {
        position: 6,
        branches: 1,
      }],
      vec![
        SiteBranches {
          position: 9,
          branches: 2,
        },
        SiteBranches {
          position: 5,
          branches: 1,
        },
      ],
    );
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_sites_by_branch_count_counts_a_branch_once_per_site() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let branches = vec![events(&["A10G", "A10G"])?, events(&["A10G"])?];
    let actual = sites_by_branch_count(&branches, &alphabet, MutationClass::Substitution);
    assert_eq!(
      vec![SiteBranches {
        position: 9,
        branches: 2
      }],
      actual
    );
    Ok(())
  }

  #[test]
  fn test_sites_by_branch_count_without_mutations_is_empty() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let branches: Vec<Vec<MutationEvent>> = vec![vec![], vec![]];
    let actual = sites_by_branch_count(&branches, &alphabet, MutationClass::Substitution);
    assert_eq!(Vec::<SiteBranches>::new(), actual);
    Ok(())
  }

  mod helpers {
    use super::*;
    use std::str::FromStr;

    pub(super) fn events(texts: &[&str]) -> Result<Vec<MutationEvent>, Report> {
      texts
        .iter()
        .map(|text| Ok(MutationEvent::Substitution(Sub::from_str(text)?)))
        .collect()
    }
  }
}
