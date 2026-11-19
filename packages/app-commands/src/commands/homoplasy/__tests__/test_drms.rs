#[cfg(test)]
mod tests {
  use crate::commands::homoplasy::drms::{DrmRow, DrmTable};
  use crate::commands::homoplasy::result::DrmAnnotation;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_io::csv::{TableFormat, csv_read};
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::listed_alt_base(               2, "T", (Some(helpers::rt_nrti(Some("M41L"))), true))]
  #[case::second_alt_base_of_a_position( 2, "C", (Some(helpers::rt_nrti(Some("M41I"))), true))]
  #[case::unlisted_alt_base(             2, "G", (Some(helpers::rt_nrti(None)),         true))]
  #[case::position_outside_the_table(    3, "T", (None,                                 false))]
  #[trace]
  fn test_drms_annotate(
    #[case] position: usize,
    #[case] base: &str,
    #[case] expected: (Option<DrmAnnotation>, bool),
  ) -> Result<(), Report> {
    let table = helpers::table()?;
    assert_eq!(expected, (table.annotate(position, base), table.contains(position)));
    Ok(())
  }

  #[test]
  fn test_drms_reject_position_zero() -> Result<(), Report> {
    let rows: Vec<DrmRow> = csv_read(
      indoc! {"
        GENOMIC_POSITION\tALT_BASE\tDRUG\tGENE\tSUBSTITUTION
        0\tT\tNRTI\tRT\tM41L
      "}
      .as_bytes(),
      TableFormat::Tsv,
    )?;
    assert_error!(
      DrmTable::from_rows(rows),
      "GENOMIC_POSITION counts from 1, but a row has position 0"
    );
    Ok(())
  }

  mod helpers {
    use crate::commands::homoplasy::drms::{DrmRow, DrmTable};
    use crate::commands::homoplasy::result::DrmAnnotation;
    use eyre::Report;
    use indoc::indoc;
    use treetime::o;
    use treetime_io::csv::{TableFormat, csv_read};

    pub(super) fn rt_nrti(substitution: Option<&str>) -> DrmAnnotation {
      DrmAnnotation {
        gene: o!("RT"),
        drug: o!("NRTI"),
        substitution: substitution.map(str::to_owned),
      }
    }

    pub(super) fn table() -> Result<DrmTable, Report> {
      let rows: Vec<DrmRow> = csv_read(
        indoc! {"
          GENOMIC_POSITION\tALT_BASE\tDRUG\tGENE\tSUBSTITUTION
          3\tT\tNRTI\tRT\tM41L
          3\tC\tNRTI\tRT\tM41I
          10\tA\tNNRTI\tRT\tK103N
        "}
        .as_bytes(),
        TableFormat::Tsv,
      )?;
      DrmTable::from_rows(rows)
    }
  }
}
