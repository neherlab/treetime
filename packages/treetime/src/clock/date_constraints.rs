use crate::make_error;
use eyre::Report;
use itertools::Itertools;
use log::{info, warn};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::date::{DateConstraint, DateValue, DatesMap};

/// The per-node date inputs [`load_date_constraints`] derives from the dates metadata, keyed by node.
///
/// Returned as values so the timetree pipeline can seed [`TimetreeState`] from them directly
/// (see [`TimetreeState::seed_from_values`]).
/// `date_constraints` is the fixed input date per node, `time_distributions` its initial posterior
/// (equal to the constraint before any date pass refines it), and `bad_branches` the exclusion flag.
/// Every node of the tree has an entry in each map.
///
/// [`TimetreeState`]: crate::timetree::timetree_state::TimetreeState
/// [`TimetreeState::seed_from_values`]: crate::timetree::timetree_state::TimetreeState::seed_from_values
#[derive(Debug, Clone, Default)]
pub struct DateConstraints {
  pub date_constraints: BTreeMap<GraphNodeKey, Option<Arc<Distribution<NegLog>>>>,
  pub time_distributions: BTreeMap<GraphNodeKey, Option<Arc<Distribution<NegLog>>>>,
  pub bad_branches: BTreeMap<GraphNodeKey, bool>,
}

pub fn date_constraint_to_distribution(constraint: &DateConstraint) -> Distribution<NegLog> {
  // A certain date carries probability 1, whose negative-log ordinate is `-ln(1) = 0`, the
  // multiplicative identity under `NegLog`. Storing `1.0` here would add a spurious constant offset
  // on every multiplication, so the ordinate is `0.0`.
  match &constraint.value {
    DateValue::Exact(d) => Distribution::point(d.value, 0.0),
    DateValue::Uncertain(r) | DateValue::Range(r) => Distribution::range((r.start, r.end), 0.0),
  }
}

#[allow(clippy::as_conversions, clippy::unwrap_used, reason = "count/index numeric cast is exact for the domain range; unwrap on a value an upstream invariant guarantees is present")]
pub fn load_date_constraints(
  dates: &DatesMap,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<DateConstraints, Report> {
  let mut good_leaf_count = 0;
  let mut bad_leaf_count = 0;
  let mut internal_constraint_count = 0;
  let mut used_names = BTreeSet::new();

  // The value maps returned to the caller; every node gets an entry. The timetree pipeline seeds
  // [`TimetreeState`] straight from these maps (see [`TimetreeState::seed_from_values`]).
  let mut date_constraints: BTreeMap<GraphNodeKey, Option<Arc<Distribution<NegLog>>>> = BTreeMap::new();
  let mut time_distributions: BTreeMap<GraphNodeKey, Option<Arc<Distribution<NegLog>>>> = BTreeMap::new();
  let mut bad_branches: BTreeMap<GraphNodeKey, bool> = BTreeMap::new();

  graph.iter_depth_first_postorder_forward(|node| {
    let key = node.key;

    let name = names[&key].clone();
    let has_constraint = name
      .as_ref()
      .and_then(|n| dates.get(n.as_str()))
      .and_then(|d| d.as_ref())
      .is_some();

    if has_constraint {
      let name = name.unwrap();
      let constraint = dates[name.as_str()].as_ref().unwrap();

      let dist = Arc::new(date_constraint_to_distribution(constraint));

      // The constraint is the input, kept as given for the whole run; the time distribution is the
      // current estimate, which starts out as the input and is refined by every inference pass.
      date_constraints.insert(key, Some(Arc::clone(&dist)));
      time_distributions.insert(key, Some(dist));
      bad_branches.insert(key, false);
      used_names.insert(name);

      if node.is_leaf {
        good_leaf_count += 1;
      } else {
        internal_constraint_count += 1;
      }
    } else if node.is_leaf {
      date_constraints.insert(key, None);
      time_distributions.insert(key, None);
      bad_branches.insert(key, true);
      bad_leaf_count += 1;
    } else {
      // Postorder guarantees every child is already recorded in the map.
      let all_children_bad = node.child_keys.iter().all(|(child_key, _)| bad_branches[child_key]);
      date_constraints.insert(key, None);
      time_distributions.insert(key, None);
      bad_branches.insert(key, all_children_bad);
    }
    Ok(())
  })?;

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

  Ok(DateConstraints {
    date_constraints,
    time_distributions,
    bad_branches,
  })
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

#[allow(clippy::as_conversions, reason = "count/index numeric cast is exact for the domain range")]
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
