use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::distribution::distribution::Distribution;
use crate::io::dates_csv::{DateOrRange, read_dates};
use crate::make_error;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::{Report, WrapErr};
use std::sync::Arc;

pub fn load_date_constraints(args: &TreetimeTimetreeArgs, graph: &GraphAncestral) -> Result<(), Report> {
  let Some(dates_path) = &args.dates else {
    return Ok(());
  };

  let dates = read_dates(dates_path, &args.name_column, &args.date_column).wrap_err("When reading dates")?;

  let mut count = 0;

  for node_ref in graph.get_nodes() {
    let node = node_ref.read_arc();
    let mut payload = node.payload().write_arc();

    let Some(name) = &payload.name else {
      continue;
    };

    let Some(date_or_range) = dates.get(name.as_str()).and_then(|d| d.as_ref()) else {
      continue;
    };

    let dist = match date_or_range {
      DateOrRange::YearFraction(t) => Arc::new(Distribution::point(*t, 1.0)),
      DateOrRange::YearFractionRange((start, end)) => Arc::new(Distribution::range((*start, *end), 1.0)),
    };

    payload.time_distribution = Some(dist);
    count += 1;
  }

  if count == 0 {
    return make_error!("No valid date constraints found in '{}'", dates_path.display());
  }

  Ok(())
}
