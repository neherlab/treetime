use crate::progress::{RunWarning, RunWarningKind};
use itertools::Itertools;
use ordered_float::OrderedFloat;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::NodeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_primitives::date::{DateConstraint, DateRange, DateValue};
use treetime_utils::fmt::string::{counted, name_list};

pub(crate) fn dates_outside_input_dates(
  graph: &Graph,
  names: &NodeNames,
  dates: &NodeMap<DateConstraint>,
  node_dates: &NodeMap<Option<f64>>,
) -> Option<RunWarning> {
  let sampling_period = sampling_period(dates)?;
  let outside = |key: GraphNodeKey, range: DateRange| {
    let date = node_dates[&key]?;
    let excess = (range.start - date).max(date - range.end);
    (excess > 0.0).then(|| OutsideDate {
      name: names[key].to_owned(),
      date,
      excess,
    })
  };

  let (dated_tips, undated_tips): (Vec<_>, Vec<_>) = graph
    .get_leaves()
    .map(|leaf| (leaf.key(), dates.get(leaf.key()).map(date_range)))
    .partition(|(_, range)| range.is_some());
  let dated_tips = dated_tips
    .into_iter()
    .filter_map(|(key, range)| outside(key, range?))
    .collect_vec();
  let undated_tips = undated_tips
    .into_iter()
    .filter_map(|(key, _)| outside(key, sampling_period))
    .collect_vec();
  let after_latest = DateRange {
    start: f64::NEG_INFINITY,
    end: sampling_period.end,
  };
  let internal_nodes = graph
    .get_internal_nodes()
    .filter_map(|node| outside(node.key(), after_latest))
    .collect_vec();

  let cases = [
    OutsideCase {
      noun: "tip",
      bound: "outside the input date range of the tip".to_owned(),
      outside: dated_tips,
    },
    OutsideCase {
      noun: "undated tip",
      bound: format!(
        "outside the sampling period {:.2} to {:.2}",
        sampling_period.start, sampling_period.end
      ),
      outside: undated_tips,
    },
    OutsideCase {
      noun: "internal node",
      bound: format!("after the latest input date {:.2}", sampling_period.end),
      outside: internal_nodes,
    },
  ]
  .into_iter()
  .filter(|case| !case.outside.is_empty())
  .collect_vec();
  if cases.is_empty() {
    return None;
  }

  let parts = cases.iter().map(OutsideCase::describe).join("; ");
  Some(RunWarning {
    kind: RunWarningKind::DatesOutsideInputDates,
    message: format!(
      "Some final dates lie outside the input dates: {parts}. The clock alone dates these nodes, as it does for the tips that the clock filter flagged and for their ancestors"
    ),
    names: cases
      .into_iter()
      .flat_map(|case| case.outside)
      .map(|outside| outside.name)
      .sorted()
      .dedup()
      .collect_vec(),
  })
}

fn sampling_period(dates: &NodeMap<DateConstraint>) -> Option<DateRange> {
  let (start, end) = dates
    .values()
    .map(date_range)
    .map(|range| (OrderedFloat(range.start), OrderedFloat(range.end)))
    .reduce(|(start, end), (next_start, next_end)| (start.min(next_start), end.max(next_end)))?;
  Some(DateRange {
    start: start.0,
    end: end.0,
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

struct OutsideCase {
  noun: &'static str,
  bound: String,
  outside: Vec<OutsideDate>,
}

impl OutsideCase {
  fn describe(&self) -> String {
    let count = counted(self.outside.len(), self.noun);
    let bound = &self.bound;
    let largest = self
      .outside
      .iter()
      .map(|outside| OrderedFloat(outside.excess))
      .max()
      .map_or(0.0, |excess| excess.0);
    let entries = self
      .outside
      .iter()
      .sorted_by(|a, b| a.name.cmp(&b.name))
      .map(|outside| format!("{} at {:.2}", outside.name, outside.date))
      .collect_vec();
    format!(
      "{count} dated {bound}, by up to {largest:.2} years ({})",
      name_list(&entries)
    )
  }
}

struct OutsideDate {
  name: String,
  date: f64,
  excess: f64,
}
