#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::load_date_constraints;
  use crate::coalescent::coalescent::compute_coalescent_model;
  use crate::o;
  use crate::partition::timetree::GraphTimetree;
  use eyre::{Report, WrapErr};
  use indexmap::IndexMap;
  use ndarray::Array1;
  use rstest::rstest;
  use serde::Deserialize;
  use std::path::Path;
  use treetime_distribution::Distribution;
  use treetime_graph::node::Named;
  use treetime_grid::grid::Grid;
  use treetime_io::dates_csv::read_dates;
  use treetime_io::nwk::nwk_read_file;
  use treetime_utils::array::serde::{array1_from_vec, indexmap_array1_from_map};
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

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
  fn test_gm_coalescent_model_matches_v0_node_contributions(#[case] snapshot_file: &str) -> Result<(), Report> {
    let snapshot: Snapshot = json_read_file(Path::new(FIXTURES_DIR).join(snapshot_file))
      .wrap_err_with(|| format!("When reading snapshot {snapshot_file}"))?;
    let graph = load_graph(&snapshot)?;
    let model = compute_coalescent_model(&graph, &Distribution::constant(snapshot.inputs.tc))?;
    let tbp_grid = Grid::from_range_n_points(
      snapshot.tbp_grid.start,
      snapshot.tbp_grid.end,
      snapshot.tbp_grid.n_points,
    )?;
    let calendar_grid = tbp_grid
      .to_array()
      .mapv(|time| snapshot.inputs.present_time - time);

    let actual: IndexMap<String, Array1<f64>> = graph
      .get_nodes()
      .iter()
      .filter_map(|node| {
        let node = node.read_arc();
        let name = node.payload().read_arc().name()?.as_ref().to_owned();
        let n_children = node.outbound().len();
        Some((name, n_children))
      })
      .map(|(name, n_children)| {
        let values = calendar_grid
          .iter()
          .map(|&time| {
            if n_children == 0 {
              Ok(model.leaf_contribution(time))
            } else {
              model.internal_contribution(time, n_children)
            }
          })
          .collect::<Result<Vec<_>, Report>>()?;
        Ok((name, Array1::from_vec(values)))
      })
      .collect::<Result<_, Report>>()?;

    // Oracle: packages/legacy/treetime/treetime/merger_models.py, captured by
    // gm_coalescent_capture without rounding deterministic floating-point values.
    pretty_assert_map_abs_diff_eq!(&snapshot.node_contributions, &actual, epsilon = 1e-6);
    Ok(())
  }

  const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/src/coalescent/__tests__/__fixtures__");

  #[derive(Debug, Deserialize)]
  struct Snapshot {
    #[allow(dead_code)]
    description: String,
    #[allow(dead_code)]
    virus_path: String,
    inputs: SnapshotInputs,
    tbp_grid: SnapshotTbpGrid,
    #[serde(deserialize_with = "array1_from_vec")]
    #[allow(dead_code)]
    lineage_counts: Array1<f64>,
    #[serde(deserialize_with = "array1_from_vec")]
    #[allow(dead_code)]
    integral_merger_rate: Array1<f64>,
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
    #[allow(dead_code)]
    padding: f64,
    n_points: usize,
  }

  fn load_graph(snapshot: &Snapshot) -> Result<GraphTimetree, Report> {
    let fixtures_dir = Path::new(FIXTURES_DIR);
    let graph = nwk_read_file(fixtures_dir.join(&snapshot.inputs.tree_path))?;
    let dates = read_dates(
      fixtures_dir.join(&snapshot.inputs.metadata_path),
      &[],
      &Some(o!("name")),
      &Some(o!("date")),
    )?;
    load_date_constraints(&dates, &graph)?;
    Ok(graph)
  }
}
