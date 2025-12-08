use crate::commands::timetree::coalescent::contributions::compute_node_contributions;
use crate::commands::timetree::coalescent::events::collect_tree_events;
use crate::commands::timetree::coalescent::integration::compute_integral_merger_rate;
use crate::commands::timetree::coalescent::lineage_dynamics::compute_lineage_count_distribution;
use crate::distribution::distribution::Distribution;
use crate::graph::node::GraphNodeKey;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::Report;
use indexmap::IndexMap;
use std::sync::Arc;

/// Computes Kingman coalescent prior contributions for all nodes in the phylogenetic tree.
///
/// Returns distributions that encode coalescent likelihood contributions for each node,
/// to be multiplied with node time distributions during backward pass optimization.
///
/// # Meaning
///
/// The Kingman coalescent model provides a probabilistic framework for assessing whether
/// the timing and structure of a phylogenetic tree are consistent with a given population
/// history. Not all tree topologies and divergence times are equally likely - trees with
/// many lineages coalescing simultaneously are less probable than gradual coalescence,
/// especially in large populations.
///
/// This function computes how likely each node's divergence time is under the coalescent
/// model, given the population size history Tc(t) and the number of concurrent lineages k(t).
/// These likelihood contributions act as priors that guide the time tree inference toward
/// more biologically plausible configurations by penalizing unlikely coalescence patterns.
///
/// # Notation and Terms
///
/// - `t` - time (negative values for past, zero at present)
/// - `k(t)` - number of concurrent lineages at time t
/// - `Tc(t)` - coalescence time scale (effective population size) at time t, controls merger rate
/// - `κ(t)` - branch merger rate: rate for one branch to merge with any other = (k(t)-1)/(2*Tc(t))
/// - `λ(t)` - total merger rate: rate for any merger to occur = k(t)*(k(t)-1)/(2*Tc(t))
/// - `I(t)` - cumulative merger rate: I(t) = ∫₀ᵗ κ(t') dt'
///
/// # Kingman Coalescent Probability Density
///
/// For internal nodes with m children (m≥2):
///
/// - Internal nodes: P(t) ∝ λ(t)^(m-1) · exp(-I(t))
///   - λ(t)^(m-1): probability density of m-way merger at time t
///   - exp(-I(t)): probability of no merger before time t
///
/// Leaf nodes do not receive coalescent contributions as they represent observed
/// samples at known times, not coalescence events.
///
/// # Returns
///
/// Map from node keys to distributions representing coalescent prior contributions.
/// Each distribution should be multiplied with the node's time distribution.
pub fn compute_coalescent_contributions(
  graph: &GraphAncestral,
  tc: &Distribution,
) -> Result<IndexMap<GraphNodeKey, Arc<Distribution>>, Report> {
  let (present_time, events_calendar) = collect_tree_events(graph)?;
  let events_tbp: Vec<_> = events_calendar
    .iter()
    .map(|(t, delta)| (present_time - *t, *delta))
    .collect();

  let lineage_counts = compute_lineage_count_distribution(&events_tbp)?;
  let integral_merger_rate = compute_integral_merger_rate(tc, &lineage_counts)?;
  compute_node_contributions(graph, &integral_merger_rate, tc, &lineage_counts, present_time)
}

#[cfg(test)]
mod tests {
  #![allow(dead_code)]
  use super::*;
  use crate::commands::timetree::data::date_constraints::load_date_constraints;
  use crate::graph::node::Named;
  use crate::io::dates_csv::read_dates;
  use crate::io::nwk::nwk_read_file;
  use crate::o;
  use eyre::{Ok, Report};
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use serde::Deserialize;
  use std::path::Path;
  use treetime_convolution::grid::Grid;
  use treetime_io::json::json_read_file;
  use treetime_utils::pretty_assert_ulps_eq;
  use treetime_utils::serde::{array1_from_vec, indexmap_array1_from_map};

  #[derive(Debug, Deserialize)]
  struct Snapshot {
    description: String,
    inputs: SnapshotInputs,
    tbp_grid: SnapshotTbpGrid,
    #[serde(deserialize_with = "array1_from_vec")]
    lineage_counts: Array1<f64>,
    #[serde(deserialize_with = "array1_from_vec")]
    integral_merger_rate: Array1<f64>,
    #[serde(deserialize_with = "array1_from_vec")]
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
    padding: f64,
    n_points: usize,
  }

  fn load_snapshot(filename: &str) -> Result<(Snapshot, GraphAncestral), Report> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
      .join("../../test_scripts/snapshots")
      .join(filename);
    let snapshot = json_read_file::<Snapshot, _>(&path)?;
    let graph = load_graph(&snapshot)?;
    Ok((snapshot, graph))
  }

  fn load_graph(snapshot: &Snapshot) -> Result<GraphAncestral, Report> {
    let tree_path = Path::new(env!("CARGO_MANIFEST_DIR"))
      .join("../../")
      .join(snapshot.inputs.tree_path.clone());
    let dates_path = Path::new(env!("CARGO_MANIFEST_DIR"))
      .join("../../")
      .join(snapshot.inputs.metadata_path.clone());
    let graph = nwk_read_file(tree_path)?;
    let dates = read_dates(dates_path, &Some(o!("name")), &Some(o!("date")))?;
    load_date_constraints(&dates, &graph)?;
    Ok(graph)
  }

  fn convert_map_keys_to_names<V>(graph: &GraphAncestral, map: &IndexMap<GraphNodeKey, V>) -> IndexMap<String, V>
  where
    V: Clone,
  {
    map
      .iter()
      .map(|(key, value)| {
        let node = graph.get_node(*key).unwrap();
        let name = node.read_arc().payload().read_arc().name().unwrap().as_ref().to_owned();
        (name, value.clone())
      })
      .collect()
  }

  #[test]
  fn test_compute_coalescent_contributions_tc0_01() -> Result<(), Report> {
    let (snapshot, graph) = load_snapshot("coalescent_flu_h3n2_20_tc0.01.json")?;
    let tc = Distribution::constant(snapshot.inputs.tc);

    let actuals = compute_coalescent_contributions(&graph, &tc)?;
    let actuals = convert_map_keys_to_names(&graph, &actuals);

    assert_eq!(
      snapshot.node_contributions.keys().collect::<Vec<_>>(),
      actuals.keys().collect::<Vec<_>>()
    );

    let t_grid = Grid::from_range_n_points(
      snapshot.tbp_grid.start,
      snapshot.tbp_grid.end,
      snapshot.tbp_grid.n_points,
    )?;

    let actuals: IndexMap<String, Array1<f64>> = actuals
      .iter()
      .map(|(name, dist)| {
        let y_resampled = dist.eval_many(&t_grid.to_array()).unwrap();
        (name.clone(), y_resampled)
      })
      .collect();

    for (node_name, expected) in &snapshot.node_contributions {
      let actual = &actuals[node_name];
      pretty_assert_ulps_eq!(expected, actual, epsilon = 1e-5, "Node: {node_name}");
    }

    Ok(())
  }
}
