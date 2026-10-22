#[cfg(test)]
mod tests {
  use crate::output_plan::OutputSelection;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::timetree::confidence::NodeConfidenceInterval;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::csv::{CsvWriter, TableFormat};
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::confidence_tsv(OutputSelection::ConfidenceTsv, Some(TableFormat::Tsv))]
  #[case::coalescent_tsv(OutputSelection::Coalescent, Some(TableFormat::Tsv))]
  #[case::confidence_csv(OutputSelection::TraitProbabilities, Some(TableFormat::Tsv))]
  #[case::traits_csv(    OutputSelection::Traits,     Some(TableFormat::Tsv))]
  #[case::clock_csv(     OutputSelection::Rtt,      Some(TableFormat::Tsv))]
  #[case::tracelog(      OutputSelection::Tracelog,      Some(TableFormat::Tsv))]
  #[case::coalescent_json(OutputSelection::CoalescentJson, None)]
  #[case::auspice(       OutputSelection::Auspice,       None)]
  #[case::nuc_fasta(     OutputSelection::Sequences, None)]
  #[trace]
  fn test_table_output_format_of_selection(#[case] selection: OutputSelection, #[case] expected: Option<TableFormat>) {
    assert_eq!(expected, selection.table_format());
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
    let mut csv = CsvWriter::new(&mut buf, OutputSelection::ConfidenceTsv.table_format().unwrap());
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
