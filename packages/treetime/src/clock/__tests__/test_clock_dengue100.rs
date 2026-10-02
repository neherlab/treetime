#[cfg(test)]
mod tests {
  use crate::cancel::NoopCancel;
  use crate::clock::clock_regression::ClockVarianceParams;
  use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootSpec};
  use crate::clock::pipeline::{self, ClockInput, ClockOutput, ClockParams};
  use crate::o;
  use crate::progress::NoopProgress;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::path::Path;
  use treetime_io::dates_csv::read_dates;
  use treetime_io::nwk::nwk_read_file;

  const DATA_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/dengue/100");

  fn run_clock(clock_params: ClockVarianceParams, clock_filter: f64, keep_root: bool) -> Result<ClockOutput, Report> {
    let data_dir = Path::new(DATA_DIR);
    let nwk_parsed = nwk_read_file(data_dir.join("tree.nwk"))?;
    let names = nwk_parsed.names();
    let dates = read_dates(
      data_dir.join("metadata.tsv"),
      &[',', '\t', ';'],
      &[],
      &Some(o!("genbank_accession")),
      &Some(o!("date")),
    )?;
    let params = ClockParams {
      clock_params,
      clock_filter,
      keep_root,
      allow_negative_rate: false,
      branch_params: BranchPointOptimizationParams::default(),
      reroot_spec: RerootSpec::default(),
    };
    let input = ClockInput {
      graph: nwk_parsed.graph,
      dates,
      branch_lengths: nwk_parsed.branch_lengths,
    };
    Ok(pipeline::run(
      &params,
      input,
      &names,
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )?)
  }

  fn get_outlier_names(output: &ClockOutput) -> Vec<String> {
    output
      .outliers
      .iter()
      .map(|key| output.names[key].clone().unwrap())
      .sorted()
      .collect()
  }

  #[test]
  fn test_dengue100_clock_pipeline_structural_properties() -> Result<(), Report> {
    let output = run_clock(ClockVarianceParams::default(), 3.0, false)?;
    let clock_model = &output.clock_model;
    let outlier_names = get_outlier_names(&output);

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

    assert!(!outlier_names.is_empty(), "Should detect at least one outlier");
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
    let output = run_clock(ClockVarianceParams::default(), 3.0, false)?;
    let clock_model = &output.clock_model;
    let outlier_names = get_outlier_names(&output);

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

    let custom_outliers = get_outlier_names(&run_clock(custom_params, 3.0, false)?);
    let default_outliers = get_outlier_names(&run_clock(ClockVarianceParams::default(), 3.0, false)?);

    assert_ne!(default_outliers, custom_outliers);
    Ok(())
  }

  #[test]
  fn test_dengue100_clock_pipeline_keep_root_allows_negative_rate() -> Result<(), Report> {
    let output = run_clock(ClockVarianceParams::default(), 0.0, true)?;
    assert!(
      output.clock_model.clock_rate() < 0.0,
      "keep-root on dengue/100 should yield a negative rate, got {:.6e}",
      output.clock_model.clock_rate()
    );
    Ok(())
  }
}
