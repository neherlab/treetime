use eyre::Report;
use std::io::Write;
use std::path::Path;
use treetime::timetree::confidence::NodeConfidenceInterval;
use treetime_io::csv::CsvStructWriter;
use treetime_utils::io::file::write_file_with;

pub fn write_confidence_intervals_file(intervals: &[NodeConfidenceInterval], filepath: &Path) -> Result<(), Report> {
  write_file_with(filepath, |file| write_confidence_intervals(intervals, file))
}

pub(crate) fn write_confidence_intervals(
  intervals: &[NodeConfidenceInterval],
  writer: impl Write + Send,
) -> Result<(), Report> {
  let mut csv = CsvStructWriter::new(writer, b'\t')?;
  intervals.iter().try_for_each(|ci| csv.write(ci))?;
  csv.into_inner()?;
  Ok(())
}

#[cfg(test)]
mod tests {
  use super::write_confidence_intervals;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime::timetree::confidence::extract_confidence_intervals;
  use treetime::timetree::inference::result::NodePosterior;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;

  #[test]
  fn test_write_confidence_intervals_omits_internal_key_column() {
    let mut graph = Graph::new();
    let mut names = BTreeMap::new();
    let key = add_named(&mut graph, &mut names, Some("named"));
    graph.build().unwrap();
    let posterior = btreemap! {
      key => NodePosterior {
        time: Some(2020.0),
        ..NodePosterior::default()
      },
    };
    let intervals = extract_confidence_intervals(&graph, &posterior, &BTreeMap::new(), &names);

    let mut buf = Vec::new();
    write_confidence_intervals(&intervals, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let header = output.lines().next().unwrap();
    assert_eq!(header, "name\tdate\tlower\tupper");
  }

  fn add_named(
    graph: &mut Graph,
    names: &mut BTreeMap<GraphNodeKey, Option<String>>,
    name: Option<&str>,
  ) -> GraphNodeKey {
    let key = graph.add_node();
    names.insert(key, name.map(|n| n.to_owned()));
    key
  }
}
