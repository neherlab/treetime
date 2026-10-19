#[cfg(test)]
pub(super) mod support {
  use crate::alphabet::alphabet::Alphabet;
  use crate::timetree::inference::result::NodePosterior;
  use deser::Deserialize;
  use eyre::Report;
  use std::collections::BTreeMap;
  use std::fs;
  use std::path::{Path, PathBuf};
  use std::sync::LazyLock;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::dates_csv::{DateConstraint, metadata_read_file};
  use treetime_io::fasta::{FastaRecord, fasta_read_file};
  use treetime_utils::io::json::json_read_str;

  const FIXTURES_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/src/timetree/inference/__tests__/__fixtures__"
  );

  pub(crate) static OUTPUTS: LazyLock<BTreeMap<String, DatasetOutputs>> = LazyLock::new(|| {
    let path = Path::new(FIXTURES_DIR).join("gm_runner_outputs.json");
    let content = fs::read_to_string(&path).expect("Failed to read gm_runner_outputs.json");
    json_read_str(&content).expect("Failed to parse gm_runner_outputs.json")
  });

  pub(crate) static ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  static INPUTS: LazyLock<BTreeMap<String, DatasetInput>> = LazyLock::new(|| {
    let path = Path::new(FIXTURES_DIR).join("gm_runner_inputs.json");
    let content = fs::read_to_string(&path).expect("Failed to read gm_runner_inputs.json");
    json_read_str(&content).expect("Failed to parse gm_runner_inputs.json")
  });

  static PROJECT_ROOT: LazyLock<PathBuf> = LazyLock::new(|| {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .expect("Failed to find project root")
      .to_path_buf()
  });

  #[derive(Debug, Deserialize)]
  pub(crate) struct DatasetOutputs {
    rerooted_tree_nwk: String,
    clock_rate: f64,
    marginal_dense: BTreeMap<String, f64>,
  }

  impl DatasetOutputs {
    pub(crate) fn rerooted_tree_nwk(&self) -> &str {
      &self.rerooted_tree_nwk
    }

    pub(crate) fn clock_rate(&self) -> f64 {
      self.clock_rate
    }

    pub(crate) fn marginal_dense(&self) -> &BTreeMap<String, f64> {
      &self.marginal_dense
    }
  }

  #[derive(Debug, Deserialize)]
  struct DatasetInput {
    aln_path: String,
    metadata_path: String,
    name_column: Option<String>,
  }

  pub(crate) fn load_dates_for_dataset(dataset: &str) -> Result<Vec<(String, Option<DateConstraint>)>, Report> {
    let input = &INPUTS[dataset];
    let metadata_path = PROJECT_ROOT.join(&input.metadata_path);
    metadata_read_file(
      &metadata_path,
      &[',', '\t', ';'],
      &treetime_io::csv::default_name_candidates(),
      input.name_column.as_deref(),
      None,
    )
    .and_then(|table| table.dates())
  }

  pub(crate) fn load_alignment_for_dataset(dataset: &str) -> Result<Vec<FastaRecord>, Report> {
    let input = &INPUTS[dataset];
    let aln_path = PROJECT_ROOT.join(&input.aln_path);
    fasta_read_file(aln_path, &*ALPHABET)
  }

  pub(crate) fn extract_node_times(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    posterior: &BTreeMap<GraphNodeKey, NodePosterior>,
  ) -> BTreeMap<String, f64> {
    graph
      .get_nodes()
      .filter_map(|node_ref| {
        let key = node_ref.key();
        let name = names[&key].clone()?;
        let time = posterior[&key].time?;
        Some((name, time))
      })
      .collect()
  }
}
