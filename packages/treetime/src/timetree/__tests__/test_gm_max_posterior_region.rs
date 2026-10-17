#[cfg(test)]
mod tests {
  use crate::timetree::confidence::{CI_FRACTION, confidence_bounds};
  use crate::timetree::max_posterior_region::{
    MaxPosteriorRegion, PeakLocation, PosteriorRegion, grid_max_posterior_region,
  };
  use deser::Deserialize;
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::path::Path;
  use treetime_utils::adapters::ArrayVec;
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::pretty_assert_abs_diff_eq;

  const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/src/timetree/__tests__/__fixtures__");

  #[rustfmt::skip]
  #[rstest]
  #[case::ebola_interior_narrow_grid( "ebola_20", "NODE_0000018")]
  #[case::ebola_interior_mid_grid(    "ebola_20", "NODE_0000017")]
  #[case::ebola_interior_wide_grid(   "ebola_20", "NODE_0000011")]
  #[case::rsv_interior_internal_node( "rsv_a_20", "NODE_0000018")]
  #[case::rsv_earliest_two_points(    "rsv_a_20", "KJ643565")]
  #[case::rsv_earliest_ok649593(      "rsv_a_20", "OK649593")]
  #[case::rsv_earliest_ok649679(      "rsv_a_20", "OK649679")]
  #[case::rsv_earliest_kp258733(      "rsv_a_20", "KP258733")]
  #[case::rsv_earliest_kp258734(      "rsv_a_20", "KP258734")]
  #[case::rsv_latest_ok649668(        "rsv_a_20", "OK649668")]
  #[case::rsv_latest_kj627732(        "rsv_a_20", "KJ627732")]
  #[case::rsv_latest_mg642040(        "rsv_a_20", "MG642040")]
  #[case::rsv_latest_ku316121(        "rsv_a_20", "KU316121")]
  #[case::rsv_latest_kj627697(        "rsv_a_20", "KJ627697")]
  #[trace]
  fn test_gm_max_posterior_region_matches_v0(#[case] dataset: &str, #[case] node: &str) -> Result<(), Report> {
    let outputs: BTreeMap<String, BTreeMap<String, GoldenNode>> =
      json_read_file(Path::new(FIXTURES_DIR).join("gm_max_posterior_region_outputs.json"))?;
    let expected = &outputs[dataset][node];
    let to_numdate = |t: f64| expected.numdate_at_zero + t / expected.clock_rate;
    let limits = (to_numdate(expected.t[0]), to_numdate(expected.t[expected.t.len() - 1]));

    let region = grid_max_posterior_region(&expected.t, &expected.neglog, CI_FRACTION)?;
    let region = PosteriorRegion {
      lower: to_numdate(region.lower),
      upper: to_numdate(region.upper),
      ..region
    };
    let without_rate = confidence_bounds(expected.numdate, limits, Some(MaxPosteriorRegion::Grid(region)), None);
    let with_rate = confidence_bounds(
      expected.numdate,
      limits,
      Some(MaxPosteriorRegion::Grid(region)),
      Some(&expected.rate_dates),
    );

    assert_eq!(expected.peak.location(), region.peak);
    pretty_assert_abs_diff_eq!(expected.region_without_rate.0, without_rate.0, epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(expected.region_without_rate.1, without_rate.1, epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(expected.region.0, with_rate.0, epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(expected.region.1, with_rate.1, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gm_max_posterior_region_shallow_posterior_holds_fraction_where_v0_returns_support() -> Result<(), Report> {
    let outputs: BTreeMap<String, BTreeMap<String, GoldenNode>> =
      json_read_file(Path::new(FIXTURES_DIR).join("gm_max_posterior_region_outputs.json"))?;
    let expected = &outputs["rsv_a_20"]["KP258696"];
    let to_numdate = |t: f64| expected.numdate_at_zero + t / expected.clock_rate;
    let support = (to_numdate(expected.t[0]), to_numdate(expected.t[expected.t.len() - 1]));
    let v0_cdf = expected
      .cdf
      .as_ref()
      .expect("interior-peak fixture node carries the v0 CDF");

    let region = grid_max_posterior_region(&expected.t, &expected.neglog, CI_FRACTION)?;

    let v0_mass =
      (interpolate(&expected.t, v0_cdf, region.upper) - interpolate(&expected.t, v0_cdf, region.lower)).abs();
    assert_eq!(PeakLocation::Interior, region.peak);
    pretty_assert_abs_diff_eq!(support.0, expected.region_without_rate.0, epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(support.1, expected.region_without_rate.1, epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(CI_FRACTION, v0_mass, epsilon = 1e-11);
    assert!(
      to_numdate(region.lower) > support.0 + 0.1,
      "region {region:?} must exclude the low-density start of the support"
    );
    Ok(())
  }

  #[derive(Deserialize)]
  struct GoldenNode {
    #[deser(as = ArrayVec)]
    t: Array1<f64>,
    #[deser(as = ArrayVec)]
    neglog: Array1<f64>,
    numdate_at_zero: f64,
    clock_rate: f64,
    numdate: f64,
    rate_dates: [f64; 3],
    peak: GoldenPeak,
    region: (f64, f64),
    region_without_rate: (f64, f64),
    #[deser(as = Option<ArrayVec>)]
    cdf: Option<Array1<f64>>,
  }

  #[derive(Clone, Copy, Deserialize)]
  #[deser(rename_all = "kebab-case")]
  enum GoldenPeak {
    Earliest,
    Interior,
    Latest,
  }

  fn interpolate(t: &Array1<f64>, y: &Array1<f64>, x: f64) -> f64 {
    let i = t.iter().position(|&ti| ti > x).unwrap_or(t.len() - 1).max(1);
    y[i - 1] + (x - t[i - 1]) / (t[i] - t[i - 1]) * (y[i] - y[i - 1])
  }

  impl GoldenPeak {
    const fn location(self) -> PeakLocation {
      match self {
        Self::Earliest => PeakLocation::Earliest,
        Self::Interior => PeakLocation::Interior,
        Self::Latest => PeakLocation::Latest,
      }
    }
  }
}
