#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::load_date_constraints;
  use crate::coalescent::coalescent::compute_coalescent_contributions;
  use crate::coalescent::time_coordinate::{CalendarTime, Tbp};
  use crate::o;
  use crate::partition::timetree::GraphTimetree;
  use eyre::{Report, WrapErr};
  use indexmap::IndexMap;
  use itertools::Itertools;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde::Deserialize;
  use std::path::Path;
  use treetime_distribution::Distribution;
  use treetime_graph::node::{GraphNodeKey, Named};
  use treetime_grid::grid::Grid;
  use treetime_io::dates_csv::read_dates;
  use treetime_io::nwk::nwk_read_file;
  use treetime_utils::array::serde::{array1_from_vec, indexmap_array1_from_map};
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::make_report;

  #[rustfmt::skip]
  #[rstest]
  #[case::flu_h3n2_20_tc_0_01(      "gm_coalescent_flu_h3n2_20_tc0.01.json")]
  #[case::flu_h3n2_20_tc_0_1(       "gm_coalescent_flu_h3n2_20_tc0.1.json")]
  #[case::flu_h3n2_20_tc_1_0(       "gm_coalescent_flu_h3n2_20_tc1.0.json")]
  #[case::flu_h3n2_20_tc_10_0(      "gm_coalescent_flu_h3n2_20_tc10.0.json")]
  #[case::ebola_20_tc_0_1(          "gm_coalescent_ebola_20_tc0.1.json")]
  #[case::ebola_20_tc_1_0(          "gm_coalescent_ebola_20_tc1.0.json")]
  #[case::ebola_20_tc_10_0(         "gm_coalescent_ebola_20_tc10.0.json")]
  #[case::dengue_20_tc_1_0(         "gm_coalescent_dengue_20_tc1.0.json")]
  #[case::dengue_20_tc_10_0(        "gm_coalescent_dengue_20_tc10.0.json")]
  #[case::dengue_20_tc_100_0(       "gm_coalescent_dengue_20_tc100.0.json")]
  #[case::rsv_a_20_tc_1_0(          "gm_coalescent_rsv_a_20_tc1.0.json")]
  #[case::rsv_a_20_tc_10_0(         "gm_coalescent_rsv_a_20_tc10.0.json")]
  #[case::rsv_a_20_tc_100_0(        "gm_coalescent_rsv_a_20_tc100.0.json")]
  #[case::mpox_clade_ii_20_tc_0_1(  "gm_coalescent_mpox_clade_ii_20_tc0.1.json")]
  #[case::mpox_clade_ii_20_tc_1_0(  "gm_coalescent_mpox_clade_ii_20_tc1.0.json")]
  #[case::mpox_clade_ii_20_tc_10_0( "gm_coalescent_mpox_clade_ii_20_tc10.0.json")]
  fn test_gm_coalescent(#[case] snapshot_file: &str) -> Result<(), Report> {
    let snapshot: Snapshot = json_read_file(Path::new(FIXTURES_DIR).join(snapshot_file))
      .wrap_err_with(|| format!("When reading snapshot {snapshot_file}"))?;

    let graph = load_graph(&snapshot)?;
    let tc = Distribution::constant(snapshot.inputs.tc);

    let actuals = compute_coalescent_contributions(&graph, &tc)?;
    let actuals = convert_map_keys_to_names(&graph, &actuals)?;

    let t_grid = Grid::from_range_n_points(
      snapshot.tbp_grid.start,
      snapshot.tbp_grid.end,
      snapshot.tbp_grid.n_points,
    )?;
    // Contributions are in calendar time; convert TBP grid to calendar for evaluation
    let present_time = CalendarTime::new(snapshot.inputs.present_time);
    let t_grid_calendar = t_grid.to_array().mapv(|t| Tbp::new(t).to_calendar(present_time).value());

    let actuals: IndexMap<String, Array1<f64>> = actuals
      .iter()
      .map(|(name, dist)| {
        let y_resampled = dist
          .eval_many(&t_grid_calendar)
          .wrap_err_with(|| format!("When evaluating contribution distribution for node '{name}'"))?;
        Ok((name.clone(), y_resampled))
      })
      .collect::<Result<_, Report>>()?;

    let expected = &snapshot.node_contributions;

    let expected_keys = expected.keys().sorted().collect_vec();
    let actual_keys = actuals.keys().sorted().collect_vec();
    assert_eq!(expected_keys, actual_keys, "Node key sets must be identical");

    let mut worst_err = 0.0_f64;
    let mut worst_node = String::new();

    for (node, exp) in expected {
      let act = &actuals[node];
      assert_eq!(
        exp.len(),
        act.len(),
        "Array length mismatch for node '{node}'"
      );

      for (i, (e, a)) in exp.iter().zip(act.iter()).enumerate() {
        let diff = (e - a).abs();
        assert!(
          diff.is_finite(),
          "Non-finite difference for node '{node}' at index {i}: expected={e}, actual={a}"
        );
        if diff > worst_err {
          worst_err = diff;
          worst_node = (*node).clone();
        }
      }
    }

    assert!(
      worst_err <= 1e-6,
      "Max absolute error {worst_err:.2e} at node '{worst_node}' exceeds threshold 1e-6"
    );

    Ok(())
  }

  const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/src/coalescent/__tests__/__fixtures__");

  #[derive(Debug, Deserialize)]
  struct Snapshot {
    /// Deserialized for identification in fixture files, not compared
    #[allow(dead_code)]
    description: String,
    /// Deserialized for identification in fixture files, not compared
    #[allow(dead_code)]
    virus_path: String,
    inputs: SnapshotInputs,
    tbp_grid: SnapshotTbpGrid,
    /// Intermediate coalescent value, deserialized for debugging; test validates final node contributions
    #[serde(deserialize_with = "array1_from_vec")]
    #[allow(dead_code)]
    lineage_counts: Array1<f64>,
    /// Intermediate coalescent value, deserialized for debugging; test validates final node contributions
    #[serde(deserialize_with = "array1_from_vec")]
    #[allow(dead_code)]
    integral_merger_rate: Array1<f64>,
    /// Intermediate coalescent value, deserialized for debugging; test validates final node contributions
    #[serde(deserialize_with = "array1_from_vec")]
    #[allow(dead_code)]
    total_merger_rate: Array1<f64>,
    #[serde(deserialize_with = "indexmap_array1_from_map")]
    node_contributions: IndexMap<String, Array1<f64>>,
  }

  #[derive(Debug, Deserialize)]
  struct SnapshotInputs {
    tree_path: String,
    metadata_path: String,
    tc: f64,
    present_time: f64,
  }

  #[derive(Debug, Deserialize)]
  struct SnapshotTbpGrid {
    start: f64,
    end: f64,
    /// Grid construction detail from v0, not needed for v1 grid reconstruction
    #[allow(dead_code)]
    padding: f64,
    n_points: usize,
  }

  fn load_graph(snapshot: &Snapshot) -> Result<GraphTimetree, Report> {
    let fixtures_dir = Path::new(FIXTURES_DIR);
    let tree_path = fixtures_dir.join(&snapshot.inputs.tree_path);
    let metadata_path = fixtures_dir.join(&snapshot.inputs.metadata_path);
    let graph = nwk_read_file(&tree_path)?;
    let dates = read_dates(&metadata_path, &[], &Some(o!("name")), &Some(o!("date")))?;
    load_date_constraints(&dates, &graph)?;
    Ok(graph)
  }

  fn convert_map_keys_to_names<V: Clone>(
    graph: &GraphTimetree,
    map: &IndexMap<GraphNodeKey, V>,
  ) -> Result<IndexMap<String, V>, Report> {
    let mut out = IndexMap::with_capacity(map.len());
    for (key, value) in map {
      let node = graph
        .get_node(*key)
        .ok_or_else(|| make_report!("Node key not found in graph: {key:?}"))?;
      let name = node
        .read_arc()
        .payload()
        .read_arc()
        .name()
        .ok_or_else(|| make_report!("Node has no name: {key:?}"))?
        .as_ref()
        .to_owned();
      out.insert(name, value.clone());
    }
    Ok(out)
  }
}
