#[cfg(test)]
mod tests {
  use crate::o;
  use eyre::Report;
  use helpers::{get_outlier_names, prefilter_outlier_names, run_clock};
  use pretty_assertions::assert_eq;
  use treetime_utils::pretty_assert_abs_diff_eq;

  const V0_CLOCK_RATE_AT_CLOCK_FILTER_3: f64 = 6.361743529953548e-04;

  #[test]
  fn test_dengue100_clock_pipeline_matches_v0() -> Result<(), Report> {
    let output = run_clock(3.0, false)?;
    let outlier_names = get_outlier_names(&output);

    pretty_assert_abs_diff_eq!(
      V0_CLOCK_RATE_AT_CLOCK_FILTER_3,
      output.clock_model.clock_rate(),
      epsilon = 1e-15
    );

    #[rustfmt::skip]
    let expected_outliers = vec![
      o!("GQ398257"), o!("GQ398268"), o!("HQ891024"), o!("KF704357"), o!("KY586699"),
      o!("OQ101603"), o!("OR389309"), o!("OR389321"), o!("OR389326"),
    ];
    assert_eq!(expected_outliers, outlier_names);

    Ok(())
  }

  #[test]
  fn test_dengue100_clock_pipeline_outliers_are_those_of_the_prefilter_fit() -> Result<(), Report> {
    let expected_outliers = prefilter_outlier_names()?;
    let outliers = get_outlier_names(&run_clock(3.0, false)?);

    assert_eq!(expected_outliers, outliers);
    Ok(())
  }

  #[test]
  fn test_dengue100_clock_pipeline_keep_root_allows_negative_rate() -> Result<(), Report> {
    let output = run_clock(0.0, true)?;
    assert!(
      output.clock_model.clock_rate() < 0.0,
      "keep-root on dengue/100 should yield a negative rate, got {:.6e}",
      output.clock_model.clock_rate()
    );
    Ok(())
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::clock::assign_dates::assign_dates;
    use crate::clock::clock_filter::clock_filter;
    use crate::clock::clock_regression::ClockTree;
    use crate::clock::clock_state::ClockInputs;
    use crate::clock::pipeline::{self, ClockInput, ClockOutput, ClockParams};
    use crate::clock::prune_tips::PruneTips;
    use crate::progress::NoopProgress;
    use crate::reroot::params::RerootParams;
    use crate::reroot::params::RerootSpec;
    use crate::test_utils::{dates_by_node, reroot_and_fit_clock};
    use eyre::Report;
    use itertools::Itertools;
    use std::collections::BTreeSet;
    use std::path::Path;
    use treetime_graph::key_map::NodeMap;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::node_names::NodeNames;
    use treetime_io::dates_csv::metadata_read_file;
    use treetime_io::nwk::{NwkParse, nwk_read_file};
    use treetime_primitives::date::DateConstraint;
    use treetime_utils::datetime::options::DateParserOptions;

    const DATA_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/dengue/100");

    pub(super) fn run_clock(clock_filter: f64, keep_root: bool) -> Result<ClockOutput, Report> {
      let (nwk_parsed, dates) = load()?;
      let names = nwk_parsed.names.clone();
      let params = ClockParams {
        clock_filter,
        allow_negative_rate: false,
        reroot_spec: (!keep_root).then(RerootSpec::default),
        prune_tips: PruneTips::default(),
      };
      let input = ClockInput {
        graph: nwk_parsed.graph,
        names,
        dates,
        branch_lengths: nwk_parsed.branch_lengths,
      };
      Ok(pipeline::run(
        &params,
        input,
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?)
    }

    pub(super) fn prefilter_outlier_names() -> Result<Vec<String>, Report> {
      let (nwk_parsed, dates) = load()?;
      let mut names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let mut inputs = ClockInputs::new(&graph);
      assign_dates(&graph, &dates, &mut inputs)?;
      let mut tree = ClockTree {
        graph,
        branch_lengths: nwk_parsed.branch_lengths,
        inputs,
      };
      let (_, result) = reroot_and_fit_clock(
        &mut tree,
        &BTreeSet::new(),
        None,
        &RerootParams::new(RerootSpec::default(), true),
        &mut names,
      )?;
      let filtered = clock_filter(
        &tree.graph,
        &tree.inputs,
        result.regression().unwrap(),
        &tree.branch_lengths,
        3.0,
        &NoopProgress,
      )?;
      Ok(sorted_names(&names, &filtered.outliers))
    }

    pub(super) fn get_outlier_names(output: &ClockOutput) -> Vec<String> {
      sorted_names(&output.names, &output.outliers)
    }

    fn load() -> Result<(NwkParse, NodeMap<DateConstraint>), Report> {
      let data_dir = Path::new(DATA_DIR);
      let nwk_parsed = nwk_read_file(data_dir.join("tree.nwk"))?;
      let dates = metadata_read_file(
        data_dir.join("metadata.tsv"),
        &[',', '\t', ';'],
        &[],
        Some("genbank_accession"),
        Some("date"),
      )
      .and_then(|table| table.dates(&DateParserOptions::default()))?;
      let dates = dates_by_node(dates, &nwk_parsed.graph, &nwk_parsed.names);
      Ok((nwk_parsed, dates))
    }

    fn sorted_names(names: &NodeNames, keys: &BTreeSet<GraphNodeKey>) -> Vec<String> {
      keys.iter().map(|key| names[*key].to_owned()).sorted().collect()
    }
  }
}
