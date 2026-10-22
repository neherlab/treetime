#[cfg(test)]
mod tests {
  use crate::output_format::table_format;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::path::Path;
  use treetime::timetree::confidence::NodeConfidenceInterval;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::csv::{CsvWriter, TableFormat};
  use treetime_utils::o;

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

  #[test]
  fn test_table_output_confidence_tsv_omits_internal_key_column() {
    let intervals = [NodeConfidenceInterval {
      key: GraphNodeKey(0),
      name: o!("named"),
      date: 2020.0,
      lower: 2019.5,
      upper: 2020.5,
    }];

    let mut buf = Vec::new();
    let mut csv = CsvWriter::new(&mut buf, TableFormat::Tsv);
    intervals
      .iter()
      .try_for_each(|interval| csv.write_row(interval))
      .unwrap();
    csv.into_inner().unwrap();
    let output = String::from_utf8(buf).unwrap();

    let header = output.lines().next().unwrap();
    assert_eq!("name\tdate\tlower\tupper", header);
  }
}
