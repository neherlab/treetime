#[cfg(test)]
mod tests {
  use crate::commands::clock::assign_dates::assign_dates;
  use crate::commands::clock::clock_filter::clock_filter_inplace;
  use crate::commands::clock::clock_graph::GraphClock;
  use crate::commands::clock::clock_model::ClockModel;
  use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
  use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::commands::clock::reroot::RerootParams;
  use crate::o;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use itertools::Itertools;
  use std::path::Path;
  use treetime_graph::node::{Named, Outlier};
  use treetime_io::dates_csv::read_dates;
  use treetime_io::nwk::nwk_read_file;

  const DATA_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/dengue/100");

  /// Load dengue/100 graph with dates assigned.
  fn load_dengue100() -> Result<GraphClock, Report> {
    let data_dir = Path::new(DATA_DIR);
    let graph: GraphClock = nwk_read_file(data_dir.join("tree.nwk"))?;
    let dates = read_dates(
      data_dir.join("metadata.tsv"),
      &Some(o!("genbank_accession")),
      &Some(o!("date")),
    )?;
    assign_dates(&graph, &dates)?;
    Ok(graph)
  }

  /// Run the full prefilter pipeline: pre-filter with force_positive=false,
  /// IQD-based outlier filtering, then final regression with force_positive=true.
  fn run_prefilter_pipeline(graph: &mut GraphClock) -> Result<(ClockModel, i32), Report> {
    let params = BranchPointOptimizationParams::default();

    // Pre-filter: allow negative rates (matching v1's estimate_clock_model_with_prefilter)
    let prefilter_reroot_params = RerootParams {
      force_positive_rate: false,
      ..RerootParams::default()
    };
    let prefilter_result = estimate_clock_model_with_reroot_policy(
      graph,
      &ClockParams::default(),
      None,
      false,
      &params,
      &prefilter_reroot_params,
      None,
    )?;
    let pre_clock_model = prefilter_result.clock_model;

    // Filter outliers
    let filter_result = clock_filter_inplace(graph, &pre_clock_model, 3.0);

    // Final regression: require positive rate
    let final_reroot_params = RerootParams::default();
    let final_result = estimate_clock_model_with_reroot_policy(
      graph,
      &ClockParams::default(),
      None,
      false,
      &params,
      &final_reroot_params,
      None,
    )?;

    Ok((final_result.clock_model, filter_result.new_outliers))
  }

  fn get_outlier_names(graph: &GraphClock) -> Vec<String> {
    graph
      .get_leaves()
      .iter()
      .filter_map(|leaf| {
        let node = leaf.read_arc();
        let payload = node.payload().read_arc();
        if payload.is_outlier() {
          payload.name().map(|n| n.as_ref().to_owned())
        } else {
          None
        }
      })
      .sorted()
      .collect()
  }

  /// Assertion-based regression test: verify the clock pipeline produces
  /// scientifically reasonable results on dengue/100.
  ///
  /// Dengue/100 is the motivating dataset for the force_positive_rate fix:
  /// all 198 root positions have negative estimated rate before outlier filtering.
  /// The pre-filter step must use force_positive_rate=false to proceed.
  #[test]
  fn test_dengue100_clock_pipeline_structural_properties() -> Result<(), Report> {
    let mut graph = load_dengue100()?;

    let (clock_model, new_outliers) = run_prefilter_pipeline(&mut graph)?;
    let outlier_names = get_outlier_names(&graph);

    // Pipeline completes with positive rate
    assert!(
      clock_model.clock_rate() > 0.0,
      "Final clock rate should be positive, got {:.6e}",
      clock_model.clock_rate()
    );

    // Rate is in plausible range for dengue (1e-5 to 1e-2 subs/site/year)
    assert!(
      clock_model.clock_rate() > 1e-5 && clock_model.clock_rate() < 1e-2,
      "Clock rate {:.6e} outside plausible dengue range [1e-5, 1e-2]",
      clock_model.clock_rate()
    );

    // R-squared is meaningful (R > 0.5 means R^2 > 0.25)
    let r_val = clock_model.r_val().expect("should have r_val");
    assert!(
      r_val > 0.5,
      "R value should indicate meaningful temporal signal, got {r_val:.4}"
    );

    // Outliers were detected
    assert!(new_outliers > 0, "Should detect at least one outlier");
    assert!(
      outlier_names.len() < 50,
      "Should not flag more than half the leaves as outliers, got {}",
      outlier_names.len()
    );

    // v0 reference outliers (7 of 8 shared with v1, 1 exclusive to v0)
    // Differences tracked in M-clock-filter-residual-parity
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

  /// Golden master test: pin v1's current output on dengue/100.
  /// When the clock filter parity issue (M-clock-filter-residual-parity) is fixed,
  /// these values should converge toward v0 and this test should be updated.
  #[test]
  fn test_dengue100_clock_pipeline_golden_master() -> Result<(), Report> {
    let mut graph = load_dengue100()?;

    let (clock_model, _) = run_prefilter_pipeline(&mut graph)?;
    let outlier_names = get_outlier_names(&graph);

    // v1 golden master values (captured from current implementation)
    assert_abs_diff_eq!(clock_model.clock_rate(), 6.787225349993138e-04, epsilon = 1e-10);
    assert_abs_diff_eq!(clock_model.intercept(), -1.116032990518721, epsilon = 1e-6);

    let r_val = clock_model.r_val().expect("should have r_val");
    assert_abs_diff_eq!(r_val, 0.810694218745354, epsilon = 1e-6);

    let chisq = clock_model.chisq().expect("should have chisq");
    assert_abs_diff_eq!(chisq, 3.297162543922308e-03, epsilon = 1e-9);

    // v1 outlier set (golden master)
    #[rustfmt::skip]
    let expected_outliers = vec![
      o!("EF105383"), o!("EF105387"), o!("GQ398268"), o!("HQ891024"), o!("KF704357"),
      o!("KY586699"), o!("MW946564"), o!("OR389309"), o!("OR389321"), o!("OR389326"),
    ];
    assert_eq!(outlier_names, expected_outliers);

    Ok(())
  }
}
