use crate::commands::clock::clock_model::ClockModel;
use crate::representation::partition::timetree::GraphTimetree;
use itertools::Itertools;
use log::warn;
use ordered_float::OrderedFloat;
use serde::Serialize;
use treetime_graph::node::Named;

#[derive(Debug, Clone, Serialize)]
pub struct OutlierRecord {
  pub name: String,
  pub given_date: f64,
  pub apparent_date: f64,
  pub residual: f64,
}

/// Collect outlier records from the graph for nodes marked as outliers.
pub fn collect_outliers(graph: &GraphTimetree, clock_model: &ClockModel, iqd: f64) -> Vec<OutlierRecord> {
  graph
    .get_leaves()
    .iter()
    .filter_map(|leaf| {
      let node = leaf.read_arc();
      let payload_arc = node.payload();
      let payload = payload_arc.read();
      if !payload.is_outlier {
        return None;
      }
      let name = payload.name().map(|n| n.as_ref().to_owned())?;
      let given_date = payload.time_distribution.as_ref()?.likely_time()?;
      let div = payload.div;
      let apparent_date = clock_model.date(div);
      let clock_deviation = clock_model.clock_deviation(given_date, div);
      let residual = if iqd > 0.0 { clock_deviation / iqd } else { 0.0 };
      Some(OutlierRecord {
        name,
        given_date,
        apparent_date,
        residual,
      })
    })
    .sorted_by_key(|r| (OrderedFloat(r.residual.abs()), r.name.clone()))
    .rev()
    .collect_vec()
}

/// Report outlier branches that violate molecular clock.
pub fn report_bad_branches(graph: &GraphTimetree, clock_model: &ClockModel, iqd: f64) {
  let outliers = collect_outliers(graph, clock_model, iqd);
  if outliers.is_empty() {
    return;
  }

  warn!("Clock filter marked {} outliers:", outliers.len());
  warn!("{:>20} {:>12} {:>14} {:>10}", "name", "given_date", "apparent_date", "residual");
  for r in &outliers {
    warn!(
      "{:>20} {:>12.2} {:>14.2} {:>10.2}",
      truncate_name(&r.name, 20),
      r.given_date,
      r.apparent_date,
      r.residual
    );
  }
}

fn truncate_name(name: &str, max_len: usize) -> String {
  if name.chars().count() <= max_len {
    name.to_owned()
  } else {
    let truncated: String = name.chars().take(max_len.saturating_sub(3)).collect();
    format!("{truncated}...")
  }
}
