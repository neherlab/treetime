/// Golden master tests for coalescent total log-likelihood.
///
/// Reference values captured from v0 `merger_model.total_LH()` with branch
/// lengths overridden to calendar time durations (parent_tbp - child_tbp)
/// to match v1's input semantics.
///
/// Capture script: `__fixtures__/gm_total_lh_capture`
/// Fixture data: `__fixtures__/gm_total_lh.json`
///
/// Binary tree cases use v0's default multiplicity=2 (identical to v1 for
/// binary trees). Polytomy cases use v0 patched with v1's multiplicity
/// formula (child's child count, 2 for leaves).
#[cfg(test)]
mod tests {
  use crate::commands::clock::date_constraints::load_date_constraints;
  use crate::commands::timetree::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::representation::partition::timetree::GraphTimetree;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use rstest::rstest;
  use treetime_distribution::Distribution;
  use treetime_io::dates_csv::DateOrRange;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::o;

  fn setup_binary_graph() -> Result<GraphTimetree, Report> {
    let dates = btreemap! {
      o!("root") => Some(DateOrRange::YearFraction(2000.0)),
      o!("internal1") => Some(DateOrRange::YearFraction(2005.0)),
      o!("leaf1") => Some(DateOrRange::YearFraction(2010.0)),
      o!("leaf2") => Some(DateOrRange::YearFraction(2010.0)),
      o!("leaf3") => Some(DateOrRange::YearFraction(2012.0)),
    };
    let graph: GraphTimetree = nwk_read_str("((leaf1:0.01,leaf2:0.01)internal1:0.01,leaf3:0.02)root:0.0;")?;
    load_date_constraints(&dates, &graph)?;
    Ok(graph)
  }

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

  // Binary tree: v0 and v1 use identical multiplicity (m=2 for all edges).
  // Reference: v0 total_LH() with time-based branch lengths.
  #[rustfmt::skip]
  #[rstest]
  #[case::tc_0_1(  0.1, -199.29621752330755)]
  #[case::tc_1(    1.0,  -19.401387711312765)]
  #[case::tc_10(  10.0,   -5.556557897319915)]
  #[case::tc_100(100.0,   -8.316728083308085)]
  #[trace]
  fn test_gm_total_lh_binary(#[case] tc: f64, #[case] expected: f64) -> Result<(), Report> {
    let graph = setup_binary_graph()?;
    let actual = compute_coalescent_total_lh(&graph, &Distribution::constant(tc))?;
    assert_abs_diff_eq!(expected, actual, epsilon = 1e-8);
    Ok(())
  }

  // Polytomy tree: internal node has 3 children. v1 uses child's child count
  // (m=3 for internal, m=2 for leaves). Reference: v0 patched with v1's formula.
  #[rustfmt::skip]
  #[rstest]
  #[case::tc_0_1(  0.1, -346.17213387796886)]
  #[case::tc_1(    1.0,  -32.812360796123585)]
  #[case::tc_10(  10.0,   -7.002587710808467)]
  #[case::tc_100(100.0,   -9.947814625459504)]
  #[trace]
  fn test_gm_total_lh_polytomy(#[case] tc: f64, #[case] expected: f64) -> Result<(), Report> {
    let graph = setup_polytomy_graph()?;
    let actual = compute_coalescent_total_lh(&graph, &Distribution::constant(tc))?;
    assert_abs_diff_eq!(expected, actual, epsilon = 1e-8);
    Ok(())
  }
}
