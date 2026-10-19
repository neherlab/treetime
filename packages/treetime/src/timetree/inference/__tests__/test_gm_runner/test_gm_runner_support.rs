#[cfg(test)]
pub(super) mod support {
  use crate::alphabet::alphabet::Alphabet;
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_regression::{ClockFitParams, ClockTree, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_state::ClockInputs;
  use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
  use crate::clock::reroot::RerootParams;
  use crate::coalescent::mean_field::CoalescentPrior;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::indel::IndelPolicy;
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, SparseReconstruction};
  use crate::partition::set::PartitionSet;
  use crate::progress::NoopProgress;
  use crate::test_utils::{dates_by_node, leaf_seq_inputs};
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::inference::bad_branches::bad_leaves;
  use crate::timetree::inference::result::{NodePosterior, PosteriorSummaries, TimeInference, given_times};
  use crate::timetree::inference::runner::{TimeInferenceInputs, run_timetree};
  use crate::timetree::optimization::relaxed_clock::unit_gammas;
  use deser::Deserialize;
  use eyre::Report;
  use itertools::Itertools;
  use std::collections::{BTreeMap, BTreeSet};
  use std::fs;
  use std::path::{Path, PathBuf};
  use std::sync::LazyLock;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::MaxGridPoints;
  use treetime_io::dates_csv::{DateConstraint, metadata_read_file};
  use treetime_io::fasta::{FastaRecord, fasta_read_file};
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;
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

  #[derive(Clone, Copy, Debug)]
  pub(crate) enum SequenceMode {
    Dense,
    Sparse,
  }

  pub(crate) struct MarginalSetup {
    pub(crate) graph: Graph,
    pub(crate) names: BTreeMap<GraphNodeKey, Option<String>>,
    pub(crate) branch_model: BranchModel,
    pub(crate) clock_model: ClockModel,
    pub(crate) constraints: DateConstraints,
    pub(crate) branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  }

  impl MarginalSetup {
    pub(crate) fn run(
      &self,
      indel_policy: IndelPolicy,
      coalescent: Option<&CoalescentPrior>,
      posterior_summaries: PosteriorSummaries,
    ) -> Result<TimeInference, Report> {
      run_timetree(
        &TimeInferenceInputs {
          graph: &self.graph,
          date_constraints: &self.constraints,
          leaf_bad_branches: &bad_leaves(&self.graph, &self.constraints, &BTreeSet::new()),
          gammas: &unit_gammas(&self.graph),
          branch_model: &self.branch_model,
          branch_lengths: &self.branch_lengths,
          names: &self.names,
          clock_model: &self.clock_model,
          clock_rate_fixed: false,
          indel_policy,
          max_grid_points: MaxGridPoints::default(),
          posterior_summaries,
        },
        coalescent,
        &NoopProgress,
      )
    }

    pub(crate) fn node_times(&self, inference: &TimeInference) -> BTreeMap<String, f64> {
      extract_node_times(&self.graph, &self.names, &inference.posterior)
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

  pub(crate) fn marginal_node_times(
    dataset: &str,
    mode: SequenceMode,
    indel_policy: IndelPolicy,
  ) -> Result<BTreeMap<String, f64>, Report> {
    let setup = marginal_setup(dataset, mode)?;
    let inference = setup.run(indel_policy, None, PosteriorSummaries::default())?;
    Ok(setup.node_times(&inference))
  }

  pub(crate) fn marginal_setup(dataset: &str, mode: SequenceMode) -> Result<MarginalSetup, Report> {
    let case = &OUTPUTS[dataset];

    let nwk_parsed = nwk_read(case.rerooted_tree_nwk().as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let dates = load_dates_for_dataset(dataset)?;
    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let aln = load_alignment_for_dataset(dataset)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect_vec();
    let reconstruction = marginal_reconstruction(&graph, &names, aln, mode)?;
    let (partition, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let times = given_times(&graph, &constraints)?;
    let clock_estimate_inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
    let (
      ClockTree {
        graph, branch_lengths, ..
      },
      clock_reroot,
    ) = estimate_clock_model_with_reroot_policy(
      ClockTree {
        graph,
        branch_lengths,
        inputs: clock_estimate_inputs,
      },
      &ClockFitParams {
        outliers: &BTreeSet::new(),
        clock_rate: Some(case.clock_rate()),
        keep_root: true,
        reroot_params: &RerootParams::default(),
        prev_clock_rate: None,
      },
      &names,
      &NoopProgress,
    )?;
    let clock_model = clock_reroot.into_clock_fit()?.model;

    Ok(MarginalSetup {
      graph,
      names,
      branch_model: BranchModel::Marginal(partition),
      clock_model,
      constraints,
      branch_lengths,
    })
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

  fn marginal_reconstruction(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: Vec<AlignmentRecord>,
    mode: SequenceMode,
  ) -> Result<PartitionSet, Report> {
    let gtr = jc69(JC69Params::default())?;
    let leaf_inputs = leaf_seq_inputs(graph, names, aln);
    Ok(match mode {
      SequenceMode::Dense => PartitionSet::Dense(vec![DenseReconstruction::seeded(
        PartitionMarginalDense::new(ALPHABET.clone(), graph, &leaf_inputs)?,
        gtr,
      )]),
      SequenceMode::Sparse => {
        let fitch = create_fitch_partition(graph, ALPHABET.clone(), leaf_inputs)?;
        let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
        PartitionSet::Sparse(vec![SparseReconstruction::seeded(partition, gtr, node_states)])
      },
    })
  }
}
