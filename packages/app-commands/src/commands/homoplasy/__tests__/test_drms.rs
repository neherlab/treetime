#[cfg(test)]
mod tests {
  use crate::commands::homoplasy::drms::{DrmRow, DrmTable};
  use crate::commands::homoplasy::result::DrmAnnotation;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use treetime::o;
  use treetime_io::csv::{TableFormat, csv_read};
  use treetime_utils::assert_error;

  #[test]
  fn test_drms_annotate_listed_alt_base() -> Result<(), Report> {
    let table = helpers::table()?;
    let expected = DrmAnnotation {
      gene: o!("RT"),
      drug: o!("NRTI"),
      substitution: Some(o!("M41L")),
    };
    assert_eq!(Some(expected), table.annotate(2, "T"));
    Ok(())
  }

  #[test]
  fn test_drms_annotate_second_alt_base_of_a_position() -> Result<(), Report> {
    let table = helpers::table()?;
    assert_eq!(Some(o!("M41I")), table.annotate(2, "C").and_then(|drm| drm.substitution));
    Ok(())
  }

  #[test]
  fn test_drms_annotate_unlisted_alt_base_without_substitution() -> Result<(), Report> {
    let table = helpers::table()?;
    let expected = DrmAnnotation {
      gene: o!("RT"),
      drug: o!("NRTI"),
      substitution: None,
    };
    assert_eq!(Some(expected), table.annotate(2, "G"));
    Ok(())
  }

  #[test]
  fn test_drms_annotate_position_outside_the_table() -> Result<(), Report> {
    let table = helpers::table()?;
    assert_eq!((None, false), (table.annotate(3, "T"), table.contains(3)));
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
    use eyre::Report;
    use indoc::indoc;
    use treetime_io::csv::{TableFormat, csv_read};

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
