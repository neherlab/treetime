#[cfg(test)]
mod tests {
  use crate::output_format::table_format;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::path::Path;
  use treetime_io::csv::TableFormat;

  #[rustfmt::skip]
  #[rstest]
  #[case::tsv(           "rtt.tsv",     Some(TableFormat::Tsv))]
  #[case::csv(           "rtt.csv",     Some(TableFormat::Csv))]
  #[case::compressed_tsv("rtt.tsv.gz",  Some(TableFormat::Tsv))]
  #[case::compressed_csv("rtt.csv.zst", Some(TableFormat::Csv))]
  #[case::stdout(        "-",           Some(TableFormat::Tsv))]
  #[case::semicolons(    "rtt.ssv",     None)]
  #[case::text(          "rtt.txt",     None)]
  #[trace]
  fn test_table_output_format_follows_the_extension(#[case] path: &str, #[case] expected: Option<TableFormat>) {
    assert_eq!(expected, table_format(Path::new(path)));
  }
}
