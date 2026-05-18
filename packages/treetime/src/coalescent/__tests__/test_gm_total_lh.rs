/// Golden master tests for coalescent total log-likelihood.
///
/// Reference values captured from v0 `merger_model.total_LH()` with branch
/// lengths overridden to calendar time durations (parent_tbp - child_tbp)
/// to match v1's input semantics.
///
/// Capture script: `__fixtures__/gm_total_lh_capture`
/// Fixture data: `__fixtures__/gm_total_lh.json`
///
/// Binary tree cases use v0's default multiplicity=2 (all formulas agree
/// for binary trees). Polytomy cases use the correct parent-based Kingman
/// formula: m = parent's child count at the merger node.
#[cfg(test)]
mod tests {
  use super::super::helpers::setup_graph;
  use crate::clock::date_constraints::load_date_constraints;
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::representation::partition::timetree::GraphTimetree;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use rstest::rstest;
  use treetime_distribution::Distribution;
  use treetime_io::dates_csv::DateOrRange;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::o;

  fn setup_polytomy_graph() -> Result<GraphTimetree, Report> {
    let dates = btreemap! {
      o!("root") => Some(DateOrRange::YearFraction(2000.0)),
      o!("internal") => Some(DateOrRange::YearFraction(2005.0)),
      o!("leaf1") => Some(DateOrRange::YearFraction(2010.0)),
      o!("leaf2") => Some(DateOrRange::YearFraction(2010.0)),
      o!("leaf3") => Some(DateOrRange::YearFraction(2010.0)),
      o!("leaf4") => Some(DateOrRange::YearFraction(2012.0)),
    };
    let graph: GraphTimetree =
      nwk_read_str("((leaf1:0.005,leaf2:0.005,leaf3:0.005)internal:0.01,leaf4:0.02)root:0.0;")?;
    load_date_constraints(&dates, &graph)?;
    Ok(graph)
  }

  // Binary tree: all formulas agree (m=2 for every edge).
  // Reference: v0 total_LH() with time-based branch lengths.
  #[rustfmt::skip]
  #[rstest]
  #[case::tc_0_1(  0.1, -199.29621752330755)]
  #[case::tc_1(    1.0,  -19.401387711312765)]
  #[case::tc_10(  10.0,   -5.556557897319915)]
  #[case::tc_100(100.0,   -8.316728083308085)]
  #[trace]
  fn test_gm_total_lh_binary(#[case] tc: f64, #[case] expected: f64) -> Result<(), Report> {
    let graph = setup_graph()?;
    let actual = compute_coalescent_total_lh(&graph, &Distribution::constant(tc))?;
    assert_abs_diff_eq!(expected, actual, epsilon = 1e-8);
    Ok(())
  }

  // Polytomy tree: internal node has 3 children. Multiplicity = parent's
  // child count at each merger node (correct Kingman formula).
  // Reference: v0 patched with m = len(node.up.clades).
  #[rustfmt::skip]
  #[rstest]
  #[case::tc_0_1(  0.1, -344.5087257790268)]
  #[case::tc_1(    1.0,  -31.916481061509884)]
  #[case::tc_10(  10.0,   -6.874236340525823)]
  #[case::tc_100(100.0,  -10.586991619508176)]
  #[trace]
  fn test_gm_total_lh_polytomy(#[case] tc: f64, #[case] expected: f64) -> Result<(), Report> {
    let graph = setup_polytomy_graph()?;
    let actual = compute_coalescent_total_lh(&graph, &Distribution::constant(tc))?;
    assert_abs_diff_eq!(expected, actual, epsilon = 1e-8);
    Ok(())
  }
}
