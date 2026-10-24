#[cfg(test)]
mod tests {
  use crate::reroot::params::{RerootMethod, RerootSpec};
  use eyre::Report;
  use helpers::{CLOCK_FILTER_IQD, root_child_lengths, run_clock};
  use treetime_utils::pretty_assert_map_ulps_eq;

  #[test]
  fn test_pipeline_reroot_oldest_gives_the_same_root_with_and_without_the_outlier_filter() -> Result<(), Report> {
    let spec = RerootSpec::Method(RerootMethod::Oldest);

    let without_filter = run_clock(&spec, 0.0)?;
    let with_filter = run_clock(&spec, CLOCK_FILTER_IQD)?;

    pretty_assert_map_ulps_eq!(
      root_child_lengths(&without_filter)?,
      root_child_lengths(&with_filter)?,
      max_ulps = 4
    );
    Ok(())
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::clock::pipeline::{self, ClockInput, ClockOutput, ClockParams};
    use crate::progress::NoopProgress;
    use crate::reroot::params::RerootSpec;
    use crate::test_utils::dates_by_node;
    use eyre::Report;
    use std::collections::BTreeMap;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::date::DateConstraint;

    pub(super) const CLOCK_FILTER_IQD: f64 = 3.0;

    const TREE: &str = "(((A:0.010,B:0.020)AB:0.010,(C:0.030,D:0.020)CD:0.015)X:0.020,\
                        ((E:0.040,F:0.050)EF:0.010,(G:0.035,H:0.060)GH:0.020)Y:0.010)R;";

    const DATES: [(&str, f64); 8] = [
      ("A", 2000.0),
      ("B", 2002.0),
      ("C", 2004.0),
      ("D", 2003.0),
      ("E", 2006.0),
      ("F", 2008.0),
      ("G", 2007.0),
      ("H", 2011.0),
    ];

    pub(super) fn run_clock(spec: &RerootSpec, clock_filter: f64) -> Result<ClockOutput, Report> {
      let parsed = nwk_read(TREE.as_bytes())?;
      let dates = DATES
        .iter()
        .map(|&(name, date)| (name.to_owned(), Some(DateConstraint::exact(date))));
      let input = ClockInput {
        dates: dates_by_node(dates, &parsed.graph, &parsed.names),
        graph: parsed.graph,
        names: parsed.names,
        branch_lengths: parsed.branch_lengths,
      };
      let params = ClockParams {
        clock_filter,
        allow_negative_rate: false,
        reroot_spec: Some(spec.clone()),
      };
      Ok(pipeline::run(
        &params,
        input,
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?)
    }

    pub(super) fn root_child_lengths(output: &ClockOutput) -> Result<BTreeMap<String, f64>, Report> {
      let schedule = GraphSchedule::new(&output.graph)?;
      let root_key = output.graph.get_exactly_one_root()?.key();
      Ok(
        schedule
          .children(root_key)
          .iter()
          .map(|&(child_key, edge_key)| {
            let length = output.branch_lengths[&edge_key].expect("every edge has a length");
            (output.names[child_key].to_owned(), length)
          })
          .collect(),
      )
    }
  }
}
