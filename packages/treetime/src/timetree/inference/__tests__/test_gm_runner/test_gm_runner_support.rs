#[cfg(test)]
pub mod support {
  //! Golden master tests for timetree inference runner.
  //!
  //! Compares Rust v1 timetree inference against Python v0 baseline node times.
  //!
  //! Three algorithms tested:
  //! - Poisson (input branch-length mode): backward/forward pass with Poisson distributions
  //! - Marginal dense: full timetree pipeline with dense sequence partition
  //! - Marginal sparse: same as dense but with sparse representation (no direct v0 analog,
  //!   validated against v0 marginal dense golden values)
  //!
  //! Golden outputs captured via `gm_runner_capture` script from v0 Python TreeTime.
  //!
  //! Tolerances are wide (3e-1 for poisson, 9e-1 for marginal) because v0 and v1 use
  //! different numerical implementations. The root node dominates the max diff in all
  //! cases - non-root nodes typically agree within 1e-2. These serve as regression
  //! guards: any code change that significantly worsens agreement will be caught.

  use crate::alphabet::alphabet::Alphabet;
  use eyre::Report;
  use serde::Deserialize;
  use std::collections::BTreeMap;
  use std::fs;
  use std::path::{Path, PathBuf};
  use std::sync::LazyLock;
  use treetime_io::dates_csv::{DatesMap, read_dates};
  use treetime_io::fasta::{FastaRecord, read_many_fasta};

  // --- Fixture types and loading ---

  const FIXTURES_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/src/timetree/inference/__tests__/__fixtures__"
  );

  pub static OUTPUTS: LazyLock<BTreeMap<String, DatasetOutputs>> = LazyLock::new(|| {
    let path = Path::new(FIXTURES_DIR).join("gm_runner_outputs.json");
    let content = fs::read_to_string(&path).expect("Failed to read gm_runner_outputs.json");
    serde_json::from_str(&content).expect("Failed to parse gm_runner_outputs.json")
  });

  pub static ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

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
  pub struct DatasetOutputs {
    rerooted_tree_nwk: String,
    clock_rate: f64,
    sequence_length: usize,
    poisson: BTreeMap<String, f64>,
    marginal_dense: BTreeMap<String, f64>,
  }

  impl DatasetOutputs {
    pub fn rerooted_tree_nwk(&self) -> &str {
      &self.rerooted_tree_nwk
    }

    pub fn clock_rate(&self) -> f64 {
      self.clock_rate
    }

    pub fn sequence_length(&self) -> usize {
      self.sequence_length
    }

    pub fn poisson(&self) -> &BTreeMap<String, f64> {
      &self.poisson
    }

    pub fn marginal_dense(&self) -> &BTreeMap<String, f64> {
      &self.marginal_dense
    }
  }

  #[derive(Debug, Deserialize)]
  struct DatasetInput {
    #[allow(dead_code)]
    tree_path: String,
    aln_path: String,
    metadata_path: String,
    name_column: Option<String>,
  }

  // --- Helper functions ---

  pub fn load_dates_for_dataset(dataset: &str) -> Result<DatesMap, Report> {
    let input = &INPUTS[dataset];
    let metadata_path = PROJECT_ROOT.join(&input.metadata_path);
    read_dates(&metadata_path, &input.name_column, &None)
  }

  pub fn load_alignment_for_dataset(dataset: &str) -> Result<Vec<FastaRecord>, Report> {
    let input = &INPUTS[dataset];
    let aln_path = PROJECT_ROOT.join(&input.aln_path);
    read_many_fasta(&[&aln_path], &*ALPHABET)
  }
}
