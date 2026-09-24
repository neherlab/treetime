use crate::make_error;
use crate::progress::ProgressSink;
use crate::{progress_info, progress_warn};
use eyre::Report;
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::date::{DateConstraint, DateValue, DatesMap};

#[allow(
  clippy::as_conversions,
  clippy::unwrap_used,
  reason = "count/index numeric cast is exact for the domain range; unwrap on a value an upstream invariant guarantees is present"
)]
pub fn load_date_constraints(
  dates: &DatesMap,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  progress: &dyn ProgressSink,
) -> Result<DateConstraints, Report> {
  let mut good_leaf_count = 0;
  let mut bad_leaf_count = 0;
  let mut internal_constraint_count = 0;
  let mut used_names = BTreeSet::new();

  let mut date_constraints: BTreeMap<GraphNodeKey, Option<Arc<Distribution<NegLog>>>> = BTreeMap::new();

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

      date_constraints.insert(key, Some(dist));
      used_names.insert(name);

      if node.is_leaf {
        good_leaf_count += 1;
      } else {
        internal_constraint_count += 1;
      }
    } else {
      date_constraints.insert(key, None);
      if node.is_leaf {
        bad_leaf_count += 1;
      }
    }
    Ok(())
  })?;

  warn_unused_date_constraints(dates, &used_names, progress);

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
    progress,
  );

  Ok(DateConstraints { date_constraints })
}

#[derive(Debug, Clone, Default)]
pub struct DateConstraints {
  pub(crate) date_constraints: BTreeMap<GraphNodeKey, Option<Arc<Distribution<NegLog>>>>,
}

impl DateConstraints {
  #[must_use]
  pub(crate) fn date_constraint(&self, key: GraphNodeKey) -> Option<Arc<Distribution<NegLog>>> {
    self.date_constraints.get(&key).cloned().flatten()
  }
}

fn date_constraint_to_distribution(constraint: &DateConstraint) -> Distribution<NegLog> {
  match &constraint.value {
    DateValue::Exact(d) => Distribution::point(d.value, 0.0),
    DateValue::Uncertain(r) | DateValue::Range(r) => Distribution::range((r.start, r.end), 0.0),
  }
}

fn warn_unused_date_constraints(dates: &DatesMap, used_names: &BTreeSet<String>, progress: &dyn ProgressSink) {
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
    progress_warn!(
      progress,
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

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn log_date_constraint_summary(
  good_leaf_count: usize,
  bad_leaf_count: usize,
  internal_constraint_count: usize,
  coverage_percent: f64,
  total_leaf_count: usize,
  progress: &dyn ProgressSink,
) {
  let bad_percent = if total_leaf_count > 0 {
    (bad_leaf_count as f64 / total_leaf_count as f64) * 100.0
  } else {
    0.0
  };

  progress_info!(progress, "Date constraint summary:");
  progress_info!(progress, "  - Total leaves: {total_leaf_count}");
  progress_info!(
    progress,
    "  - Leaves with dates: {good_leaf_count} ({coverage_percent:.1}%)"
  );
  progress_info!(
    progress,
    "  - Leaves without dates: {bad_leaf_count} ({bad_percent:.1}%)"
  );
  if internal_constraint_count > 0 {
    progress_info!(progress, "  - Internal nodes with dates: {internal_constraint_count}");
  }

  if bad_percent > 50.0 {
    progress_warn!(
      progress,
      "More than half of leaves lack date constraints. This may affect inference quality."
    );
  }
}
