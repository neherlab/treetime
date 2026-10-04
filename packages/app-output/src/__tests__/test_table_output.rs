#[cfg(test)]
mod tests {
  use crate::output_plan::OutputSelection;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime::timetree::confidence::extract_confidence_intervals;
  use treetime::timetree::inference::result::NodePosterior;
  use treetime_graph::graph::Graph;
  use treetime_io::csv::{CsvWriter, TableFormat};

  #[rustfmt::skip]
  #[rstest]
  #[case::confidence_tsv(OutputSelection::ConfidenceTsv, Some(TableFormat::Tsv))]
  #[case::coalescent_tsv(OutputSelection::CoalescentTsv, Some(TableFormat::Tsv))]
  #[case::confidence_csv(OutputSelection::ConfidenceCsv, Some(TableFormat::Csv))]
  #[case::traits_csv(    OutputSelection::TraitsCsv,     Some(TableFormat::Csv))]
  #[case::clock_csv(     OutputSelection::ClockCsv,      Some(TableFormat::Csv))]
  #[case::tracelog(      OutputSelection::Tracelog,      Some(TableFormat::Csv))]
  #[case::coalescent_csv(OutputSelection::CoalescentCsv, Some(TableFormat::Csv))]
  #[case::coalescent_json(OutputSelection::CoalescentJson, None)]
  #[case::auspice(       OutputSelection::Auspice,       None)]
  #[case::nuc_fasta(     OutputSelection::ReconstructedNucFasta, None)]
  #[trace]
  fn test_table_output_format_of_selection(#[case] selection: OutputSelection, #[case] expected: Option<TableFormat>) {
    assert_eq!(expected, selection.table_format());
  }

  #[test]
  fn test_table_output_confidence_tsv_omits_internal_key_column() {
    let mut graph = Graph::new();
    let key = graph.add_node();
    graph.build().unwrap();
    let names = btreemap! { key => Some("named".to_owned()) };
    let posterior = btreemap! {
      key => NodePosterior {
        time: Some(2020.0),
        ..NodePosterior::default()
      },
    };
    let intervals = extract_confidence_intervals(&graph, &posterior, &BTreeMap::new(), &names);

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
