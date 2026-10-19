use crate::progress::{RunWarning, RunWarningKind};
use itertools::Itertools;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use treetime_graph::assign_node_names::node_name_or_key;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::date::{DateConstraint, DateRange, DateValue};
use treetime_utils::fmt::string::{counted, name_list};

pub(crate) fn dates_outside_input_dates(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  dates: &BTreeMap<GraphNodeKey, DateConstraint>,
  node_dates: &BTreeMap<GraphNodeKey, Option<f64>>,
) -> Option<RunWarning> {
  let latest_input_date = dates.values().map(|date| OrderedFloat(date_range(date).end)).max()?.0;
  let name = |key: GraphNodeKey| node_name_or_key(key, names[&key].as_deref());

  let tips = graph
    .get_leaves()
    .filter_map(|leaf| {
      let range = date_range(dates.get(&leaf.key())?);
      let date = node_dates[&leaf.key()]?;
      let excess = (range.start - date).max(date - range.end);
      (excess > 0.0).then(|| (name(leaf.key()), excess))
    })
    .collect_vec();
  let internal_nodes = graph
    .get_internal_nodes()
    .filter_map(|node| {
      let excess = node_dates[&node.key()]? - latest_input_date;
      (excess > 0.0).then(|| (name(node.key()), excess))
    })
    .collect_vec();
  if tips.is_empty() && internal_nodes.is_empty() {
    return None;
  }

  let cases = [
    (
      &tips,
      counted(tips.len(), "tip"),
      "outside the input date range of the tip".to_owned(),
    ),
    (
      &internal_nodes,
      counted(internal_nodes.len(), "internal node"),
      format!("after the latest input date {latest_input_date:.2}"),
    ),
  ];
  let parts = cases
    .iter()
    .filter(|(outside, _, _)| !outside.is_empty())
    .map(|(outside, count, bound)| {
      let largest = outside
        .iter()
        .map(|(_, excess)| OrderedFloat(*excess))
        .max()
        .map_or(0.0, |excess| excess.0);
      let names = outside.iter().map(|(name, _)| name.clone()).sorted().collect_vec();
      format!(
        "{count} dated {bound}, by up to {largest:.2} years ({})",
        name_list(&names)
      )
    })
    .join("; ");
  Some(RunWarning {
    kind: RunWarningKind::DatesOutsideInputDates,
    message: format!(
      "Some final dates lie outside the input dates: {parts}. The clock alone dates these nodes, as it does for the tips that the clock filter flagged and for their ancestors"
    ),
    names: tips
      .into_iter()
      .chain(internal_nodes)
      .map(|(name, _)| name)
      .sorted()
      .dedup()
      .collect_vec(),
  })
}

fn date_range(date: &DateConstraint) -> DateRange {
  match date.value {
    DateValue::Exact(exact) => DateRange {
      start: exact.value,
      end: exact.value,
    },
    DateValue::Uncertain(range) | DateValue::Range(range) => range,
  }
}
