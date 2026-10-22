#[cfg(test)]
mod tests {
  use crate::commands::clock::args::{TreetimeClockArgs, TreetimeClockArgsRaw};
  use crate::commands::clock::run::run_clock;
  use crate::commands::shared::metadata::DateColumnArgs;
  use crate::commands::shared::output_args::{ClockOutputSelection, OutputCoreArgs};
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::fs;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_utils::datetime::options::DateParserOptions;
  use treetime_utils::datetime::parse_date::parse_date;
  use treetime_utils::datetime::year_fraction::date_to_year_fraction;

  const TREE: &str = "((A:0.1,B:0.2):0.1,(C:0.3,D:0.4):0.1);";

  #[rustfmt::skip]
  #[rstest]
  #[case::built_in_reads_day_first( ["03/06/2016",  "07/04/2017",  "01/02/2018",  "11/12/2019"],  None,             ["2016-06-03", "2017-04-07", "2018-02-01", "2019-12-11"])]
  #[case::format_reads_month_first( ["03/06/2016",  "07/04/2017",  "01/02/2018",  "11/12/2019"],  Some("%m/%d/%Y"), ["2016-03-06", "2017-07-04", "2018-01-02", "2019-11-12"])]
  #[case::format_reads_month_name(  ["13 Jun 2016", "04 Jul 2017", "02 Jan 2018", "12 Nov 2019"], Some("%d %b %Y"), ["2016-06-13", "2017-07-04", "2018-01-02", "2019-11-12"])]
  #[trace]
  fn test_clock_reads_sampling_dates_in_the_date_format(
    #[case] dates: [&str; 4],
    #[case] date_format: Option<&str>,
    #[case] expected_iso: [&str; 4],
  ) -> Result<(), Report> {
    let dir = tempfile::tempdir()?;
    let tree = dir.path().join("tree.nwk");
    let metadata = dir.path().join("metadata.tsv");
    fs::write(&tree, format!("{TREE}\n"))?;
    let rows = ["A", "B", "C", "D"].iter().zip(dates).map(|(name, date)| format!("{name}\t{date}\n")).join("");
    fs::write(&metadata, format!("strain\tdate\n{rows}"))?;
    let date_column = date_format.map_or_else(DateColumnArgs::default, |date_format| DateColumnArgs {
      date_format: date_format.to_owned(),
      ..DateColumnArgs::default()
    });
    let raw = TreetimeClockArgsRaw {
      tree: Some(tree),
      metadata: Some(metadata),
      date_column,
      keep_root: true,
      output: OutputCoreArgs {
        output_all: Some(dir.path().join("output")),
        ..OutputCoreArgs::default()
      },
      output_selection: vec![ClockOutputSelection::Nwk],
      ..TreetimeClockArgsRaw::default()
    };

    let result = run_clock(&TreetimeClockArgs::try_from(raw)?, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let expected = ["A", "B", "C", "D"]
      .iter()
      .zip(expected_iso)
      .map(|(name, iso)| Ok(((*name).to_owned(), Some(date_to_year_fraction(&parse_date(iso, &DateParserOptions::default())?)))))
      .collect::<Result<Vec<_>, Report>>()?;
    let actual = result
      .regression_results
      .into_iter()
      .filter(|row| row.is_leaf)
      .map(|row| (row.name, row.date))
      .sorted_by(|left, right| left.0.cmp(&right.0))
      .collect_vec();
    assert_eq!(expected, actual);
    Ok(())
  }
}
