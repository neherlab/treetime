#[cfg(test)]
mod tests {
  use crate::clock::clock_filter::clock_filter_inplace;
  use crate::clock::clock_graph::GraphClock;
  use crate::clock::clock_model::ClockModel;
  use crate::o;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::node::{Named, Outlier};
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;

  /// Build an 8-leaf balanced tree where 6 leaves follow a positive clock and 2 are extreme
  /// outliers (high divergence, mid-range dates). Outlier branches G=2.0 and H=3.0 produce
  /// root-to-tip distances ~2.02 and ~3.02, far above the clock-expected ~0.15.
  fn setup_outlier_graph() -> Result<GraphClock, Report> {
    let tree = "(((A:0.1,B:0.2):0.01,(C:0.15,D:0.25):0.01):0.01,((E:0.12,F:0.18):0.01,(G:2.0,H:3.0):0.01):0.01)root;";
    let graph: GraphClock = nwk_read_str(tree)?;

    // Good clock: rate=0.01/year, base=2000 → date = div/0.01 + 2000
    // Root-to-tip: A=0.12, B=0.22, C=0.17, D=0.27, E=0.14, F=0.20, G=2.02, H=3.02
    #[rustfmt::skip]
    let dates = btreemap! {
      o!("A") => 2012.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2027.0,
      o!("E") => 2014.0,
      o!("F") => 2020.0,
      o!("G") => 2015.0, // outlier: div=2.02 but date implies div~0.15
      o!("H") => 2015.0, // outlier: div=3.02 but date implies div~0.15
    };

    for n in graph.get_leaves() {
      let name = n.read_arc().payload().read_arc().name().unwrap().as_ref().to_owned();
      n.write_arc().payload().write_arc().time = Some(dates[&name]);
    }

    Ok(graph)
  }

  fn get_outlier_names(graph: &GraphClock) -> Vec<String> {
    let mut names: Vec<String> = graph
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
      .collect();
    names.sort();
    names
  }

  #[test]
  fn test_clock_filter_positive_rate_identifies_outliers() -> Result<(), Report> {
    let graph = setup_outlier_graph()?;
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let result = clock_filter_inplace(&graph, &clock_model, 3.0)?;

    assert!(result.iqd > 0.0, "IQD should be positive");
    let outliers = get_outlier_names(&graph);
    assert_eq!(outliers, vec![o!("G"), o!("H")]);

    Ok(())
  }

  #[test]
  fn test_clock_filter_negative_rate_identifies_same_outliers() -> Result<(), Report> {
    let graph = setup_outlier_graph()?;
    // Negative rate model: slope inverted, intercept adjusted.
    // IQD-based filtering uses |deviation| > IQD*threshold, so the absolute-value
    // comparison makes outlier detection slope-sign-invariant for extreme outliers.
    let clock_model = ClockModel::for_testing(-0.005, 10.5);

    let result = clock_filter_inplace(&graph, &clock_model, 3.0)?;

    assert!(result.iqd > 0.0, "IQD should be positive");
    let outliers = get_outlier_names(&graph);
    assert_eq!(outliers, vec![o!("G"), o!("H")]);

    Ok(())
  }

  #[test]
  fn test_clock_filter_rejects_no_dated_leaves() -> Result<(), Report> {
    let graph = helpers::setup_low_cardinality_graph(0)?;
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let result = clock_filter_inplace(&graph, &clock_model, 3.0);

    assert_error!(result, "Clock filtering requires at least one dated leaf");
    Ok(())
  }

  // The v0 NumPy oracle accepts non-empty residual arrays of every cardinality.
  // packages/legacy/treetime/treetime/clock_filter_methods.py#L14
  #[rustfmt::skip]
  #[rstest]
  #[case::one_dated_leaf(  1)]
  #[case::two_dated_leaves(2)]
  #[case::three_dated_leaves(3)]
  #[trace]
  fn test_clock_filter_accepts_low_cardinality_input(#[case] dated_leaf_count: usize) -> Result<(), Report> {
    let graph = helpers::setup_low_cardinality_graph(dated_leaf_count)?;
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let result = clock_filter_inplace(&graph, &clock_model, 3.0)?;

    assert!(result.iqd.is_finite());
    Ok(())
  }

  mod helpers {
    use crate::clock::clock_graph::GraphClock;
    use eyre::Report;
    use treetime_io::nwk::nwk_read_str;

    pub fn setup_low_cardinality_graph(dated_leaf_count: usize) -> Result<GraphClock, Report> {
      let graph: GraphClock = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;

      graph
        .get_leaves()
        .iter()
        .take(dated_leaf_count)
        .enumerate()
        .for_each(|(index, leaf)| {
          leaf.write_arc().payload().write_arc().time = Some(2000.0 + index as f64);
        });

      Ok(graph)
    }
  }
}
