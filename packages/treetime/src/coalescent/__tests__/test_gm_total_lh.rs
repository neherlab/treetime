#[cfg(test)]
mod tests {
  use super::super::helpers::{coalescent_node_times, setup_graph};
  use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use rstest::rstest;
  use treetime_distribution::Distribution;
  use treetime_graph::graph::Graph;
  use treetime_io::dates_csv::DateConstraint;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::o;

  fn setup_polytomy_graph() -> Result<(Graph, DateConstraints), Report> {
    let dates = btreemap! {
      o!("root") => Some(DateConstraint::exact(2000.0)),
      o!("internal") => Some(DateConstraint::exact(2005.0)),
      o!("leaf1") => Some(DateConstraint::exact(2010.0)),
      o!("leaf2") => Some(DateConstraint::exact(2010.0)),
      o!("leaf3") => Some(DateConstraint::exact(2010.0)),
      o!("leaf4") => Some(DateConstraint::exact(2012.0)),
    };
    let nwk_parsed = nwk_read_str("((leaf1:0.005,leaf2:0.005,leaf3:0.005)internal:0.01,leaf4:0.02)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let constraints = load_date_constraints(&dates, &graph, &names)?;
    Ok((graph, constraints))
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tc_0_1(  0.1, -199.29621752330755)]
  #[case::tc_1(    1.0,  -19.401387711312765)]
  #[case::tc_10(  10.0,   -5.556557897319915)]
  #[case::tc_100(100.0,   -8.316728083308085)]
  #[trace]
  fn test_gm_total_lh_binary(#[case] tc: f64, #[case] expected: f64) -> Result<(), Report> {
    let (graph, names, constraints) = setup_graph()?;
    let actual = compute_coalescent_total_lh(&graph, &Distribution::constant(tc), &coalescent_node_times(&graph, &constraints))?.value();
    assert_abs_diff_eq!(expected, actual, epsilon = 1e-8);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tc_0_1(  0.1, -344.5087257790268)]
  #[case::tc_1(    1.0,  -31.916481061509884)]
  #[case::tc_10(  10.0,   -6.874236340525823)]
  #[case::tc_100(100.0,  -10.586991619508176)]
  #[trace]
  fn test_gm_total_lh_polytomy(#[case] tc: f64, #[case] expected: f64) -> Result<(), Report> {
    let (graph, constraints) = setup_polytomy_graph()?;
    let actual = compute_coalescent_total_lh(&graph, &Distribution::constant(tc), &coalescent_node_times(&graph, &constraints))?.value();
    assert_abs_diff_eq!(expected, actual, epsilon = 1e-8);
    Ok(())
  }
}
