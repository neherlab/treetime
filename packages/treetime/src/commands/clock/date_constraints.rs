use crate::distribution::distribution::Distribution;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named, TimeConstraint};
use crate::io::dates_csv::DatesMap;
use crate::make_error;
use crate::o;
use eyre::Report;
use itertools::Itertools;
use log::{info, warn};
use std::collections::BTreeSet;
use std::sync::Arc;

pub trait DateConstraintNode: GraphNode + Named + TimeConstraint<Arc<Distribution>> {}

pub fn load_date_constraints<N, E, D>(dates: &DatesMap, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: DateConstraintNode,
  E: GraphEdge,
  D: Sync + Send,
{
  let mut good_leaf_count = 0;
  let mut bad_leaf_count = 0;
  let mut internal_constraint_count = 0;
  let mut used_names = BTreeSet::new();

  graph.iter_depth_first_postorder_forward(|node| {
    let mut payload = node.payload;

    let name = payload.name().map(|n| o!(n.as_ref()));
    let has_constraint = name
      .as_ref()
      .and_then(|n| dates.get(n.as_str()))
      .and_then(|d| d.as_ref())
      .is_some();

    if has_constraint {
      let name = name.unwrap();
      let date_or_range = dates[name.as_str()].as_ref().unwrap();

      let dist = Arc::new(Distribution::from_date_or_range(date_or_range));

      payload.set_time_distribution(Some(dist));
      payload.set_bad_branch(false);
      used_names.insert(name);

      if node.is_leaf {
        good_leaf_count += 1;
      } else {
        internal_constraint_count += 1;
      }
    } else if node.is_leaf {
      payload.set_bad_branch(true);
      bad_leaf_count += 1;
    } else {
      let all_children_bad = node
        .children
        .iter()
        .all(|(child_payload, _)| child_payload.read_arc().bad_branch());
      payload.set_bad_branch(all_children_bad);
    }
  });

  warn_unused_date_constraints(dates, &used_names);

  let total_leaf_count = good_leaf_count + bad_leaf_count;
  let coverage_percent = if total_leaf_count > 0 {
    (good_leaf_count as f64 / total_leaf_count as f64) * 100.0
  } else {
    0.0
  };

  validate_minimum_date_constraints(good_leaf_count, total_leaf_count, coverage_percent)?;

  log_date_constraint_summary(
    good_leaf_count,
    bad_leaf_count,
    internal_constraint_count,
    coverage_percent,
    total_leaf_count,
  );

  Ok(())
}

fn warn_unused_date_constraints(dates: &DatesMap, used_names: &BTreeSet<String>) {
  let unused_names: Vec<_> = dates
    .keys()
    .filter(|name| !used_names.contains(name.as_str()))
    .collect();

  if !unused_names.is_empty() {
    let sample = unused_names
      .iter()
      .take(10)
      .map(|s| s.as_str())
      .collect_vec()
      .join(", ");
    let suffix = if unused_names.len() > 10 { "..." } else { "" };
    warn!(
      "Date constraints found for {} names not present in tree: {}{}",
      unused_names.len(),
      sample,
      suffix
    );
  }
}

fn validate_minimum_date_constraints(
  good_leaf_count: usize,
  total_leaf_count: usize,
  coverage_percent: f64,
) -> Result<(), Report> {
  if good_leaf_count < 3 {
    return make_error!(
      "Insufficient dated leaves: found {} out of {} ({:.1}% coverage, minimum 3 required). Need {} more dated leaves.",
      good_leaf_count,
      total_leaf_count,
      coverage_percent,
      3 - good_leaf_count
    );
  }
  Ok(())
}

fn log_date_constraint_summary(
  good_leaf_count: usize,
  bad_leaf_count: usize,
  internal_constraint_count: usize,
  coverage_percent: f64,
  total_leaf_count: usize,
) {
  let bad_percent = if total_leaf_count > 0 {
    (bad_leaf_count as f64 / total_leaf_count as f64) * 100.0
  } else {
    0.0
  };

  info!("Date constraint summary:");
  info!("  - Total leaves: {total_leaf_count}");
  info!("  - Leaves with dates: {good_leaf_count} ({coverage_percent:.1}%)");
  info!("  - Leaves without dates: {bad_leaf_count} ({bad_percent:.1}%)");
  if internal_constraint_count > 0 {
    info!("  - Internal nodes with dates: {internal_constraint_count}");
  }

  if bad_percent > 50.0 {
    warn!("More than half of leaves lack date constraints. This may affect inference quality.");
  }
}
