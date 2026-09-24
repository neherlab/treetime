#[cfg(test)]
pub(super) mod support {
  use crate::alphabet::alphabet::Alphabet;
  use crate::timetree::timetree_state::TimetreeState;
  use eyre::Report;
  use ndarray::Array1;
  use ordered_float::OrderedFloat;
  use serde::Deserialize;
  use std::collections::BTreeMap;
  use std::fs;
  use std::path::{Path, PathBuf};
  use std::sync::Arc;
  use std::sync::LazyLock;
  use treetime_distribution::{Distribution, DistributionFunction, NegLog};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::dates_csv::{DatesMap, read_dates};
  use treetime_io::fasta::{FastaRecord, read_many_fasta_path};

  const FIXTURES_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/src/timetree/inference/__tests__/__fixtures__"
  );

  const MIN_TIME_MUTATION_FRACTION: f64 = 0.01;

  pub(crate) static OUTPUTS: LazyLock<BTreeMap<String, DatasetOutputs>> = LazyLock::new(|| {
    let path = Path::new(FIXTURES_DIR).join("gm_runner_outputs.json");
    let content = fs::read_to_string(&path).expect("Failed to read gm_runner_outputs.json");
    serde_json::from_str(&content).expect("Failed to parse gm_runner_outputs.json")
  });

  pub(crate) static ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  static INPUTS: LazyLock<BTreeMap<String, DatasetInput>> = LazyLock::new(|| {
    let path = Path::new(FIXTURES_DIR).join("gm_runner_inputs.json");
    let content = fs::read_to_string(&path).expect("Failed to read gm_runner_inputs.json");
    serde_json::from_str(&content).expect("Failed to parse gm_runner_inputs.json")
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
    sequence_length: usize,
    poisson: BTreeMap<String, f64>,
    marginal_dense: BTreeMap<String, f64>,
  }

  impl DatasetOutputs {
    pub(crate) fn rerooted_tree_nwk(&self) -> &str {
      &self.rerooted_tree_nwk
    }

    pub(crate) fn clock_rate(&self) -> f64 {
      self.clock_rate
    }

    pub(crate) fn sequence_length(&self) -> usize {
      self.sequence_length
    }

    pub(crate) fn poisson(&self) -> &BTreeMap<String, f64> {
      &self.poisson
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

  pub(crate) fn load_dates_for_dataset(dataset: &str) -> Result<DatesMap, Report> {
    let input = &INPUTS[dataset];
    let metadata_path = PROJECT_ROOT.join(&input.metadata_path);
    read_dates(
      &metadata_path,
      &[',', '\t', ';'],
      &treetime_io::csv::default_name_candidates(),
      &input.name_column,
      &None,
    )
  }

  pub(crate) fn load_alignment_for_dataset(dataset: &str) -> Result<Vec<FastaRecord>, Report> {
    let input = &INPUTS[dataset];
    let aln_path = PROJECT_ROOT.join(&input.aln_path);
    read_many_fasta_path(&[&aln_path], &*ALPHABET)
  }

  pub(crate) fn extract_node_times(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    state: &TimetreeState,
  ) -> BTreeMap<String, f64> {
    graph
      .get_nodes()
      .filter_map(|node_ref| {
        let key = node_ref.key();
        let name = names[&key].clone()?;
        let time = state.nodes.get(&key).and_then(|node| node.time)?;
        Some((name, time))
      })
      .collect()
  }

  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  pub(crate) fn create_poisson_branch_distributions(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    mu: f64,
    seq_len: usize,
    n_points: usize,
  ) -> Result<BTreeMap<GraphEdgeKey, Arc<Distribution<NegLog>>>, Report> {
    let seq_len_f64 = seq_len as f64;

    let mut distributions = BTreeMap::new();
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();

      if let Some(branch_length) = branch_lengths[&edge_key] {
        let expected_time = branch_length / mu;
        let max_time = 3.0 * expected_time.max(1.0);

        let min_time = MIN_TIME_MUTATION_FRACTION / (mu * seq_len_f64);
        let grid = Array1::linspace(min_time, max_time, n_points);

        let log_p = grid.mapv(|dt| -dt * mu * seq_len_f64 + branch_length * seq_len_f64 * (dt * mu * seq_len_f64).ln());

        let log_p_max = log_p.iter().copied().map(OrderedFloat).max().map_or(0.0, |x| x.0);
        let neg_log = log_p.mapv(|value| log_p_max - value);

        let distribution_fn = DistributionFunction::from_range_values((min_time, max_time), neg_log)?;
        let distribution = Distribution::Function(distribution_fn);
        distributions.insert(edge_key, Arc::new(distribution));
      }
    }

    Ok(distributions)
  }
}
