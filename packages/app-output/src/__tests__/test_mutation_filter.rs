#[cfg(test)]
mod tests {
  use crate::mutation_filter::UnknownMutationFilter;
  use eyre::Report;
  use helpers::{reported, substitutions};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_primitives::AsciiChar;

  #[rustfmt::skip]
  #[rstest]
  #[case::default_drops_unknown(  false, vec!["C3K", "A4T"])]
  #[case::report_ambiguous(       true,  vec!["C1N", "N2A", "C3K", "A4T"])]
  #[trace]
  fn test_mutation_filter_nucleotide_unknown(
    #[case] report_ambiguous: bool,
    #[case] expected: Vec<&str>,
  ) -> Result<(), Report> {
    let filter = UnknownMutationFilter::new(AsciiChar::try_new(b'N')?, report_ambiguous);
    let actual = reported(filter, substitutions(&["C1N", "N2A", "C3K", "A4T"])?)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_mutation_filter_amino_acid_unknown() -> Result<(), Report> {
    let filter = UnknownMutationFilter::new(AsciiChar::try_new(b'X')?, false);
    let actual = reported(filter, substitutions(&["A1X", "X2N", "K3B"])?)?;
    assert_eq!(vec!["K3B"], actual);
    Ok(())
  }

  mod helpers {
    use crate::mutation_filter::UnknownMutationFilter;
    use eyre::Report;
    use std::str::FromStr;
    use treetime::seq::mutation::{MutationEvent, Sub, mutation_event_strings};

    pub(super) fn substitutions(strings: &[&str]) -> Result<Vec<MutationEvent>, Report> {
      strings
        .iter()
        .map(|string| Ok(MutationEvent::Substitution(Sub::from_str(string)?)))
        .collect()
    }

    pub(super) fn reported(filter: UnknownMutationFilter, events: Vec<MutationEvent>) -> Result<Vec<String>, Report> {
      Ok(
        filter
          .reported_events(events)
          .iter()
          .map(mutation_event_strings)
          .collect::<Result<Vec<_>, _>>()?
          .into_iter()
          .flatten()
          .collect(),
      )
    }
  }
}
