use crate::output_plan::OutputSelection;
use crate::table_output::table_create;
use eyre::Report;
use indexmap::IndexMap;
use ndarray::Array1;
use std::collections::BTreeMap;
use std::io::Write;
use std::iter::once;
use std::path::Path;
use treetime::mugration::pipeline::MugrationOutput;
use treetime_graph::node::GraphNodeKey;
use treetime_io::csv::CsvWriter;

#[derive(Debug)]
pub struct MugrationResult {
  pub nodes: BTreeMap<GraphNodeKey, MugrationNodeOut>,
  pub traits: MugrationTraitsOutput,
  pub confidence: MugrationConfidenceOutput,
}

impl MugrationResult {
  pub fn new(
    output: &MugrationOutput,
    input_confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    attribute: &str,
  ) -> Self {
    let nodes: BTreeMap<GraphNodeKey, MugrationNodeOut> = output
      .graph
      .get_nodes()
      .map(|node| {
        let key = node.key();
        (
          key,
          MugrationNodeOut {
            name: names.get(&key).cloned().flatten(),
            branch_support: input_confidences.get(&key).copied().flatten(),
          },
        )
      })
      .collect();
    let assignments = extract_trait_assignments(output, names);
    let traits = MugrationTraitsOutput::new(attribute, assignments);
    let confidence = MugrationConfidenceOutput::new(output, names);

    Self {
      nodes,
      traits,
      confidence,
    }
  }
}

#[derive(Clone, Debug)]
pub struct MugrationConfidenceOutput {
  states: Vec<String>,
  rows: Vec<ConfidenceRow>,
}

impl MugrationConfidenceOutput {
  fn new(output: &MugrationOutput, names: &BTreeMap<GraphNodeKey, Option<String>>) -> Self {
    let states: Vec<String> = output.states.iter().map(|s| s.to_owned()).collect();

    let rows: Vec<ConfidenceRow> = output
      .graph
      .get_nodes()
      .filter_map(|node| {
        let node_key = node.key();
        let node_name = node_name_or_fallback(names, node_key);

        output.confidences[&node_key].clone().map(|profile| ConfidenceRow {
          node: node_name,
          profile,
        })
      })
      .collect();

    Self { states, rows }
  }

  pub fn write_csv_file(&self, filepath: &Path) -> Result<(), Report> {
    let mut csv = table_create(OutputSelection::ConfidenceCsv, filepath)?;
    self.write_rows(&mut csv)?;
    csv.finish()
  }

  fn write_rows<W: Write>(&self, csv: &mut CsvWriter<W>) -> Result<(), Report> {
    csv.write_record(once("node").chain(self.states.iter().map(String::as_str)))?;
    self.rows.iter().try_for_each(|row| {
      let probs = row.profile.iter().map(|p| format!("{p:.6}"));
      csv.write_record(once(row.node.clone()).chain(probs))
    })
  }
}

#[derive(Clone, Debug)]
pub struct ConfidenceRow {
  node: String,
  profile: Array1<f64>,
}

#[derive(Clone, Debug)]
pub struct MugrationTraitsOutput {
  pub(crate) attribute: String,
  pub(crate) assignments: IndexMap<String, String>,
}

impl MugrationTraitsOutput {
  fn new(attribute: &str, assignments: IndexMap<String, String>) -> Self {
    Self {
      attribute: attribute.to_owned(),
      assignments,
    }
  }

  pub fn write_csv_file(&self, filepath: &Path) -> Result<(), Report> {
    let mut csv = table_create(OutputSelection::TraitsCsv, filepath)?;
    self.write_rows(&mut csv)?;
    csv.finish()
  }

  fn write_rows<W: Write>(&self, csv: &mut CsvWriter<W>) -> Result<(), Report> {
    csv.write_record(["node", self.attribute.as_str()])?;
    self
      .assignments
      .iter()
      .map(<[&String; 2]>::from)
      .try_for_each(|record| csv.write_record(record))
  }
}

#[derive(Debug, Clone)]
pub struct MugrationNodeOut {
  pub(crate) name: Option<String>,
  pub(crate) branch_support: Option<f64>,
}

fn extract_trait_assignments(
  output: &MugrationOutput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> IndexMap<String, String> {
  output
    .graph
    .get_nodes()
    .filter_map(|node| {
      let node_key = node.key();
      let node_name = node_name_or_fallback(names, node_key);

      output.reconstructed_traits[&node_key]
        .clone()
        .map(|trait_value| (node_name, trait_value))
    })
    .collect()
}

fn node_name_or_fallback(names: &BTreeMap<GraphNodeKey, Option<String>>, node_key: GraphNodeKey) -> String {
  names[&node_key]
    .clone()
    .unwrap_or_else(|| format!("node_{}", node_key.0))
}

#[cfg(test)]
mod tests {
  use super::{ConfidenceRow, MugrationConfidenceOutput, MugrationTraitsOutput};
  use indexmap::indexmap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_io::csv::{CsvWriter, TableFormat};
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::plain(            ("country",       "A",     "usa"),              "node,country\nA,usa\n")]
  #[case::comma_in_value(   ("country",       "A",     "Congo, DR"),        "node,country\nA,\"Congo, DR\"\n")]
  #[case::quote_in_value(   ("country",       "A",     "Cote \"d'Ivoire\""), "node,country\nA,\"Cote \"\"d'Ivoire\"\"\"\n")]
  #[case::newline_in_value( ("country",       "A",     "a\nb"),             "node,country\nA,\"a\nb\"\n")]
  #[case::comma_in_node(    ("country",       "A/1,2", "usa"),              "node,country\n\"A/1,2\",usa\n")]
  #[case::comma_in_header(  ("region, state", "A",     "usa"),              "node,\"region, state\"\nA,usa\n")]
  #[trace]
  fn test_mugration_result_traits_csv_quotes_fields_per_rfc4180(
    #[case] (attribute, node, value): (&str, &str, &str),
    #[case] expected: &str,
  ) {
    let traits = MugrationTraitsOutput::new(attribute, indexmap! { o!(node) => o!(value) });

    let mut buf = Vec::new();
    let mut csv = CsvWriter::new(&mut buf, TableFormat::Csv);
    traits.write_rows(&mut csv).unwrap();
    csv.into_inner().unwrap();

    assert_eq!(expected, String::from_utf8(buf).unwrap());
  }

  #[test]
  fn test_mugration_result_confidence_csv_quotes_fields_per_rfc4180() {
    let confidence = MugrationConfidenceOutput {
      states: vec![o!("Congo, DR"), o!("usa")],
      rows: vec![
        ConfidenceRow {
          node: o!("A,1"),
          profile: array![0.25, 0.75],
        },
        ConfidenceRow {
          node: o!("B"),
          profile: array![1.0, 0.0],
        },
      ],
    };

    let mut buf = Vec::new();
    let mut csv = CsvWriter::new(&mut buf, TableFormat::Csv);
    confidence.write_rows(&mut csv).unwrap();
    csv.into_inner().unwrap();

    let expected = "\
      node,\"Congo, DR\",usa\n\
      \"A,1\",0.250000,0.750000\n\
      B,1.000000,0.000000\n\
    ";
    assert_eq!(expected, String::from_utf8(buf).unwrap());
  }
}
