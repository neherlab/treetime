#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::partition::fitch::site_changes::fitch_site_change_counts;
  use crate::test_utils::leaf_seq_inputs;
  use eyre::Report;
  use helpers::record;
  use pretty_assertions::assert_eq;
  use treetime_io::nwk::nwk_read;

  #[test]
  fn test_fitch_site_change_counts_count_the_parsimony_changes_of_each_site() -> Result<(), Report> {
    let parsed = nwk_read(b"((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root;".as_slice())?;
    let sequences = vec![
      record("A", "AAAAA")?,
      record("B", "ACCAC")?,
      record("C", "AAACG")?,
      record("D", "AACCT")?,
    ];
    let inputs = leaf_seq_inputs(&parsed.graph, &parsed.names, sequences);

    let actual = fitch_site_change_counts(&parsed.graph, Alphabet::new(AlphabetName::Nuc)?, inputs)?;

    let expected = vec![0, 1, 2, 1, 3];
    assert_eq!(
      expected, actual,
      "Fitch parsimony scores of the five sites, counted by hand on the tree"
    );
    Ok(())
  }

  mod helpers {
    use eyre::Report;
    use treetime_primitives::{AlignmentRecord, Seq};

    pub(super) fn record(name: &str, sequence: &str) -> Result<AlignmentRecord, Report> {
      Ok(AlignmentRecord {
        name: name.to_owned(),
        seq: Seq::try_from_str(sequence)?,
      })
    }
  }
}
