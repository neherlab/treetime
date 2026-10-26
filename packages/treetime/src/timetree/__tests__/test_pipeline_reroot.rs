#[cfg(test)]
mod tests {
  use crate::error::OperationError;
  use crate::reroot::params::{RerootMethod, RerootSpec};
  use eyre::Report;
  use helpers::{
    MISSING_LENGTH_TREE, NEGATIVE_RATE_TREE, SYMMETRIC_TREE, X_TIPS, find_tips, outlier_names, root_child_lengths,
    run_timetree, run_timetree_operation,
  };
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_utils::error::report_to_string;
  use treetime_utils::{o, pretty_assert_map_ulps_eq, pretty_assert_ulps_eq};

  #[test]
  fn test_pipeline_reroot_date_free_keeps_the_clock_refits_of_the_pre_loop() -> Result<(), Report> {
    let least_squares = run_timetree(SYMMETRIC_TREE, |_| RerootSpec::Method(RerootMethod::LeastSquares))?;
    let outgroup = run_timetree(SYMMETRIC_TREE, |tree| RerootSpec::Tips(find_tips(tree, &X_TIPS)))?;

    let input_root = btreemap! { o!("X") => 0.010, o!("Y") => 0.010 };
    pretty_assert_map_ulps_eq!(input_root, root_child_lengths(&least_squares)?, max_ulps = 4);
    pretty_assert_map_ulps_eq!(input_root, root_child_lengths(&outgroup)?, max_ulps = 4);
    assert!(
      !least_squares.outliers.is_empty(),
      "the outlier filter must flag a tip in this case"
    );
    assert_eq!(outlier_names(&least_squares), outlier_names(&outgroup));
    pretty_assert_ulps_eq!(
      least_squares.clock_model.clock_rate(),
      outgroup.clock_model.clock_rate(),
      max_ulps = 4
    );
    pretty_assert_ulps_eq!(
      least_squares.clock_model.intercept(),
      outgroup.clock_model.intercept(),
      max_ulps = 4
    );
    Ok(())
  }

  #[test]
  fn test_pipeline_reroot_date_free_root_without_a_positive_regression_rate_starts_from_the_gaussian_rate()
  -> Result<(), Report> {
    let output = run_timetree(NEGATIVE_RATE_TREE, |tree| {
      RerootSpec::Tips(find_tips(tree, &["A", "B"]))
    })?;

    let input_root = btreemap! { o!("X") => 0.010, o!("Y") => 0.010 };
    pretty_assert_map_ulps_eq!(input_root, root_child_lengths(&output)?, max_ulps = 4);
    assert!(
      output.clock_model.clock_rate() > 0.0,
      "the clock refits of the pre-loop on the kept root fall back to the positive Gaussian rate, got {}",
      output.clock_model.clock_rate()
    );
    Ok(())
  }

  #[test]
  fn test_pipeline_reroot_min_dev_rejects_a_missing_branch_length_as_an_input_error() -> Result<(), Report> {
    let result = run_timetree_operation(MISSING_LENGTH_TREE, |_| RerootSpec::Method(RerootMethod::MinDev))?;

    let Err(OperationError::InvalidInput(report)) = result else {
      panic!("a missing branch length must be an input error");
    };
    assert_eq!(
      "Failed to infer clock model: --reroot=min-dev places the root from branch lengths and requires every edge to \
       have a finite, non-negative branch length. Invalid edges:\n  AB -> B: missing\nGive the tree branch lengths, \
       for example with `treetime optimize`",
      report_to_string(&report)
    );
    Ok(())
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::error::OperationError;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::progress::NoopProgress;
    use crate::reroot::params::RerootSpec;
    use crate::test_utils::{dates_by_node, find_node_key_by_name, marginal_timetree_params};
    use crate::timetree::params::{BranchLengthMode, TimetreeParams};
    use crate::timetree::pipeline::{self, SequenceLength, TimetreeInput, TimetreeOutput, TimetreeSequences};
    use eyre::Report;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::{NwkParse, nwk_read};
    use treetime_primitives::date::DateConstraint;

    pub(super) const SYMMETRIC_TREE: &str = "((A:0.014,B:0.010,C:0.018,D:0.015,E:0.025,F:0.023,O:0.0575)X:0.010,\
                                             (P:0.014,Q:0.010,S:0.018,T:0.015,U:0.025,V:0.023,W:0.0575)Y:0.010)R;";

    pub(super) const X_TIPS: [&str; 7] = ["A", "B", "C", "D", "E", "F", "O"];

    pub(super) const NEGATIVE_RATE_TREE: &str = "((A:0.040,B:0.030)X:0.010,(C:0.020,D:0.010)Y:0.010)R;";

    pub(super) const MISSING_LENGTH_TREE: &str = "((A:0.010,B)AB:0.010,(C:0.010,D:0.020)CD:0.010)R;";

    const DATES: [(&str, f64); 14] = [
      ("A", 2000.0),
      ("B", 2003.0),
      ("C", 2006.0),
      ("D", 2009.0),
      ("E", 2012.0),
      ("F", 2015.0),
      ("O", 2007.5),
      ("P", 2000.0),
      ("Q", 2003.0),
      ("S", 2006.0),
      ("T", 2009.0),
      ("U", 2012.0),
      ("V", 2015.0),
      ("W", 2007.5),
    ];

    const SEQUENCE_LENGTH: usize = 1000;

    const CLOCK_FILTER_IQD: f64 = 3.0;

    pub(super) fn run_timetree(
      newick: &str,
      spec: impl FnOnce(&NwkParse) -> RerootSpec,
    ) -> Result<TimetreeOutput, Report> {
      run_timetree_operation(newick, spec)?.map_err(OperationError::into_report)
    }

    pub(super) fn run_timetree_operation(
      newick: &str,
      spec: impl FnOnce(&NwkParse) -> RerootSpec,
    ) -> Result<Result<TimetreeOutput, OperationError>, Report> {
      let tree = nwk_read(newick.as_bytes())?;
      let params = TimetreeParams {
        branch_length_mode: BranchLengthMode::Input,
        reroot_spec: Some(spec(&tree)),
        clock_filter: CLOCK_FILTER_IQD,
        prune_filtered_tips: None,
        ..marginal_timetree_params()
      };
      let dates = DATES
        .iter()
        .map(|&(name, date)| (name.to_owned(), Some(DateConstraint::exact(date))));
      let input = TimetreeInput {
        dates: Some(dates_by_node(dates, &tree.graph, &tree.names)),
        sequences: TimetreeSequences::Length(SequenceLength {
          length: Some(SEQUENCE_LENGTH),
          model: GtrModel::JC69(JC69Params::default()),
          alphabet: None,
        }),
        graph: tree.graph,
        names: tree.names,
        branch_lengths: tree.branch_lengths,
      };
      Ok(pipeline::run(
        &params,
        input,
        None,
        None,
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      ))
    }

    pub(super) fn find_tips(tree: &NwkParse, tips: &[&str]) -> Vec<GraphNodeKey> {
      tips
        .iter()
        .map(|tip| find_node_key_by_name(&tree.graph, &tree.names, tip).expect("the tip exists"))
        .collect()
    }

    pub(super) fn outlier_names(output: &TimetreeOutput) -> BTreeSet<String> {
      output
        .outliers
        .iter()
        .map(|&key| output.names[key].to_owned())
        .collect()
    }

    pub(super) fn root_child_lengths(output: &TimetreeOutput) -> Result<BTreeMap<String, f64>, Report> {
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
