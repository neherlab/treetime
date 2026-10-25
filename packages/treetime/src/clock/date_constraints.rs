use crate::error::input_error;
use crate::progress::LogSink;
use crate::{progress_info, progress_warn};
use eyre::Report;
use std::sync::Arc;
use treetime_distribution::{Distribution, NegLog};
use treetime_graph::graph::Graph;
use treetime_graph::key_map::NodeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_primitives::date::{DateConstraint, DateValue};
use treetime_utils::cast::ToF64Exact;

pub fn load_date_constraints(
  dates: &NodeMap<DateConstraint>,
  graph: &Graph,
  log: &dyn LogSink,
) -> Result<DateConstraints, Report> {
  let mut good_leaf_count = 0;
  let mut bad_leaf_count = 0;
  let mut internal_constraint_count = 0;

  let mut date_constraints: NodeMap<Option<Arc<Distribution<NegLog>>>> = NodeMap::new();

  let schedule = GraphSchedule::new(graph)?;
  for &key in schedule.postorder() {
    let is_leaf = schedule.children(key).is_empty();

    if let Some(constraint) = dates.get(key) {
      let dist = Arc::new(date_constraint_to_distribution(constraint));

      date_constraints.insert(key, Some(dist));

      if is_leaf {
        good_leaf_count += 1;
      } else {
        internal_constraint_count += 1;
      }
    } else {
      date_constraints.insert(key, None);
      if is_leaf {
        bad_leaf_count += 1;
      }
    }
  }

  let total_leaf_count = good_leaf_count + bad_leaf_count;
  let coverage_percent = if total_leaf_count > 0 {
    (good_leaf_count.to_f64_exact() / total_leaf_count.to_f64_exact()) * 100.0
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
    log,
  );

  Ok(DateConstraints {
    by_node: date_constraints,
  })
}

#[derive(Debug, Clone, Default)]
pub struct DateConstraints {
  pub(crate) by_node: NodeMap<Option<Arc<Distribution<NegLog>>>>,
}

impl DateConstraints {
  #[must_use]
  pub(crate) fn date_constraint(&self, key: GraphNodeKey) -> Option<Arc<Distribution<NegLog>>> {
    self.by_node.get(key).cloned().flatten()
  }
}

fn date_constraint_to_distribution(constraint: &DateConstraint) -> Distribution<NegLog> {
  match &constraint.value {
    DateValue::Exact(d) => Distribution::point(d.value, 0.0),
    DateValue::Uncertain(r) | DateValue::Range(r) => Distribution::range((r.start, r.end), 0.0),
  }
}

fn validate_minimum_date_constraints(
  good_leaf_count: usize,
  total_leaf_count: usize,
  coverage_percent: f64,
) -> Result<(), Report> {
  if good_leaf_count < 3 {
    return Err(input_error(format!(
      "Insufficient dated leaves: found {} out of {} ({:.1}% coverage, minimum 3 required). Need {} more dated leaves.",
      good_leaf_count,
      total_leaf_count,
      coverage_percent,
      3 - good_leaf_count
    )));
  }
  Ok(())
}

fn log_date_constraint_summary(
  good_leaf_count: usize,
  bad_leaf_count: usize,
  internal_constraint_count: usize,
  coverage_percent: f64,
  total_leaf_count: usize,
  log: &dyn LogSink,
) {
  let bad_percent = if total_leaf_count > 0 {
    (bad_leaf_count.to_f64_exact() / total_leaf_count.to_f64_exact()) * 100.0
  } else {
    0.0
  };

  progress_info!(log, "Date constraint summary:");
  progress_info!(log, "  - Total leaves: {total_leaf_count}");
  progress_info!(log, "  - Leaves with dates: {good_leaf_count} ({coverage_percent:.1}%)");
  progress_info!(log, "  - Leaves without dates: {bad_leaf_count} ({bad_percent:.1}%)");
  if internal_constraint_count > 0 {
    progress_info!(log, "  - Internal nodes with dates: {internal_constraint_count}");
  }

  if bad_percent > 50.0 {
    progress_warn!(
      log,
      "More than half of leaves lack date constraints. This may affect inference quality."
    );
  }
}
