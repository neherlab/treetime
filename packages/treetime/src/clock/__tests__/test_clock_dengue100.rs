#[cfg(test)]
mod tests {
  use crate::cancel::NoopCancel;
  use crate::clock::assign_dates::assign_dates;
  use crate::clock::clock_filter::clock_filter_inplace;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_regression::{ClockVarianceParams, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_state::{ClockInputs, ClockState};
  use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
  use crate::clock::pipeline::{self, ClockInput, ClockParams};
  use crate::clock::reroot::RerootParams;
  use crate::o;
  use crate::progress::NoopProgress;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::path::Path;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::dates_csv::read_dates;
  use treetime_io::nwk::nwk_read_file;

  const DATA_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/dengue/100");

  fn load_dengue100() -> Result<
    (
      Graph,
      BTreeMap<GraphNodeKey, Option<String>>,
      ClockInputs,
      ClockState,
      BTreeMap<GraphEdgeKey, Option<f64>>,
    ),
    Report,
  > {
    let data_dir = Path::new(DATA_DIR);
    let nwk_parsed = nwk_read_file(data_dir.join("tree.nwk"))?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let dates = read_dates(
      data_dir.join("metadata.tsv"),
      &[',', '\t', ';'],
      &[],
      &Some(o!("genbank_accession")),
      &Some(o!("date")),
    )?;
    let mut inputs = ClockInputs::new(&graph);
    assign_dates(&graph, &dates, &mut inputs, &names)?;
    let state = ClockState::new(&graph);
    Ok((graph, names, inputs, state, branch_lengths))
  }

  fn run_prefilter_pipeline(
    graph: &mut Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    inputs: &mut ClockInputs,
    state: &mut ClockState,
    clock_params: &ClockVarianceParams,
    branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<(ClockModel, i32), Report> {
    let params = BranchPointOptimizationParams::default();

    let prefilter_reroot_params = RerootParams {
      force_positive_rate: false,
      ..RerootParams::default()
    };
    let names_tt_2 = names.clone();
    let (new_state, prefilter_result) = estimate_clock_model_with_reroot_policy(
      graph,
      inputs,
      std::mem::take(state),
      clock_params,
      None,
      false,
      &params,
      &prefilter_reroot_params,
      branch_lengths,
      None,
      &names_tt_2,
    )?;
    *state = new_state;
    let pre_regression = prefilter_result.regression();

    let filter_result = clock_filter_inplace(graph, inputs, state, pre_regression, branch_lengths, 3.0)?;

    let final_reroot_params = RerootParams::default();
    let names_tt_1 = names.clone();
    let (new_state, final_result) = estimate_clock_model_with_reroot_policy(
      graph,
      inputs,
      std::mem::take(state),
      clock_params,
      None,
      false,
      &params,
      &final_reroot_params,
      branch_lengths,
      None,
      &names_tt_1,
    )?;
    *state = new_state;

    Ok((final_result.into_clock_model()?, filter_result.new_outliers))
  }

  fn get_outlier_names(
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    graph: &Graph,
    state: &ClockState,
  ) -> Vec<String> {
    graph
      .get_leaves()
      .filter_map(|leaf| {
        let node = leaf;
        if state.node(node.key()).is_outlier {
          names.get(&node.key()).cloned().flatten()
        } else {
          None
        }
      })
      .sorted()
      .collect()
  }

  #[test]
  fn test_dengue100_clock_pipeline_structural_properties() -> Result<(), Report> {
    let (mut graph, names, mut inputs, mut state, mut branch_lengths) = load_dengue100()?;

    let (clock_model, new_outliers) = run_prefilter_pipeline(
      &mut graph,
      &names,
      &mut inputs,
      &mut state,
      &ClockVarianceParams::default(),
      &mut branch_lengths,
    )?;
    let outlier_names = get_outlier_names(&names, &graph, &state);

    assert!(
      clock_model.clock_rate() > 0.0,
      "Final clock rate should be positive, got {:.6e}",
      clock_model.clock_rate()
    );

    assert!(
      clock_model.clock_rate() > 1e-5 && clock_model.clock_rate() < 1e-2,
      "Clock rate {:.6e} outside plausible dengue range [1e-5, 1e-2]",
      clock_model.clock_rate()
    );

    let r_val = clock_model.r_val().expect("should have r_val");
    assert!(
      r_val > 0.5,
      "R value should indicate meaningful temporal signal, got {r_val:.4}"
    );

    assert!(new_outliers > 0, "Should detect at least one outlier");
    assert!(
      outlier_names.len() < 50,
      "Should not flag more than half the leaves as outliers, got {}",
      outlier_names.len()
    );

    #[rustfmt::skip]
    let v0_outliers = [o!("GQ398257"), o!("GQ398268"), o!("HQ891024"), o!("KF704357"),
      o!("KY586699"), o!("OR389309"), o!("OR389321"), o!("OR389326")];
    let shared: Vec<_> = v0_outliers.iter().filter(|name| outlier_names.contains(name)).collect();
    assert!(
      shared.len() >= 6,
      "Should share at least 6 of 8 v0 outliers, got {}: {shared:?}",
      shared.len()
    );

    Ok(())
  }

  #[test]
  fn test_dengue100_clock_pipeline_golden_master() -> Result<(), Report> {
    let (mut graph, names, mut inputs, mut state, mut branch_lengths) = load_dengue100()?;

    let (clock_model, _) = run_prefilter_pipeline(
      &mut graph,
      &names,
      &mut inputs,
      &mut state,
      &ClockVarianceParams::default(),
      &mut branch_lengths,
    )?;
    let outlier_names = get_outlier_names(&names, &graph, &state);

    assert_abs_diff_eq!(clock_model.clock_rate(), 6.787225349993138e-04, epsilon = 1e-10);
    assert_abs_diff_eq!(clock_model.intercept(), -1.116032990518721, epsilon = 1e-6);

    let r_val = clock_model.r_val().expect("should have r_val");
    assert_abs_diff_eq!(r_val, 0.810694218745354, epsilon = 1e-6);

    let chisq = clock_model.chisq().expect("should have chisq");
    assert_abs_diff_eq!(chisq, 3.297162543922308e-03, epsilon = 1e-9);

    #[rustfmt::skip]
    let expected_outliers = vec![
      o!("EF105383"), o!("EF105387"), o!("GQ398268"), o!("HQ891024"), o!("KF704357"),
      o!("KY586699"), o!("MW946564"), o!("OR389309"), o!("OR389321"), o!("OR389326"),
    ];
    assert_eq!(outlier_names, expected_outliers);

    Ok(())
  }

  #[test]
  fn test_dengue100_clock_pipeline_prefilter_uses_supplied_clock_params() -> Result<(), Report> {
    let custom_params = ClockVarianceParams {
      variance_factor: 1e-3,
      variance_offset: 0.0,
      variance_offset_leaf: 1e-4,
    };

    let (mut expected_graph, expected_names, mut expected_inputs, mut expected_state, mut expected_branch_lengths) =
      load_dengue100()?;
    let (_expected_clock_model, _) = run_prefilter_pipeline(
      &mut expected_graph,
      &expected_names,
      &mut expected_inputs,
      &mut expected_state,
      &custom_params,
      &mut expected_branch_lengths,
    )?;
    let expected_outliers = get_outlier_names(&expected_names, &expected_graph, &expected_state);

    let (mut default_graph, default_names, mut default_inputs, mut default_state, mut default_branch_lengths) =
      load_dengue100()?;
    let (_default_clock_model, _) = run_prefilter_pipeline(
      &mut default_graph,
      &default_names,
      &mut default_inputs,
      &mut default_state,
      &ClockVarianceParams::default(),
      &mut default_branch_lengths,
    )?;
    let default_outliers = get_outlier_names(&default_names, &default_graph, &default_state);
    assert_ne!(default_outliers, expected_outliers);

    let data_dir = Path::new(DATA_DIR);
    let nwk_parsed = nwk_read_file(data_dir.join("tree.nwk"))?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let dates = read_dates(
      data_dir.join("metadata.tsv"),
      &[',', '\t', ';'],
      &[],
      &Some(o!("genbank_accession")),
      &Some(o!("date")),
    )?;
    let params = ClockParams {
      clock_params: custom_params,
      clock_filter: 3.0,
      keep_root: false,
      allow_negative_rate: false,
      branch_params: BranchPointOptimizationParams::default(),
      reroot_spec: RerootSpec::default(),
    };
    let output = pipeline::run(
      &params,
      ClockInput {
        graph,
        dates,
        branch_lengths,
      },
      &names,
      &NoopCancel,
      &NoopProgress,
    )?;
    let actual_outliers = get_outlier_names(&names, &output.graph, &output.state);

    assert_eq!(expected_outliers, actual_outliers);
    Ok(())
  }

  #[test]
  fn test_dengue100_clock_pipeline_keep_root_allows_negative_rate() -> Result<(), Report> {
    let data_dir = Path::new(DATA_DIR);
    let nwk_parsed = nwk_read_file(data_dir.join("tree.nwk"))?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let dates = read_dates(
      data_dir.join("metadata.tsv"),
      &[',', '\t', ';'],
      &[],
      &Some(o!("genbank_accession")),
      &Some(o!("date")),
    )?;
    let params = ClockParams {
      clock_params: ClockVarianceParams::default(),
      clock_filter: 0.0,
      keep_root: true,
      allow_negative_rate: false,
      branch_params: BranchPointOptimizationParams::default(),
      reroot_spec: RerootSpec::default(),
    };
    let output = pipeline::run(
      &params,
      ClockInput {
        graph,
        dates,
        branch_lengths,
      },
      &names,
      &NoopCancel,
      &NoopProgress,
    )?;
    assert!(
      output.clock_model.clock_rate() < 0.0,
      "keep-root on dengue/100 should yield a negative rate, got {:.6e}",
      output.clock_model.clock_rate()
    );
    Ok(())
  }
}
