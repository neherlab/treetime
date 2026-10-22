use crate::annotated_graph::{AnnotatedTreeView, TreeDates};
use crate::table_output::planned_table_format;
use eyre::{Report, WrapErr};
use std::io::Write;
use std::iter::once;
use std::path::Path;
use treetime_graph::node::GraphNodeKey;
use treetime_io::csv::{CsvWriter, TableFormat};
use treetime_utils::datetime::year_fraction::year_fraction_to_datestring;
use treetime_utils::io::file::write_file_with;
use treetime_utils::make_internal_report;

const CONFIDENCE_COMMENT: &str = "#Lower and upper bound delineate the 90% max posterior region";

const MISSING_DATE: &str = "--";

pub fn write_dates_table_file(tree: &AnnotatedTreeView<'_>, path: &Path) -> Result<(), Report> {
  let format = planned_table_format(path)?;
  write_file_with(path, |writer| write_dates_table(tree, writer, format)).wrap_err("When writing the node dates table")
}

pub(crate) fn write_dates_table<W: Write>(
  tree: &AnnotatedTreeView<'_>,
  mut writer: W,
  format: TableFormat,
) -> Result<(), Report> {
  let graph = tree.graph();
  let dates = graph
    .dates
    .as_ref()
    .ok_or_else(|| make_internal_report!("The node dates table requires inferred node dates"))?;
  let with_confidence = dates.confidence.is_some();
  if with_confidence {
    writeln!(writer, "{CONFIDENCE_COMMENT}").wrap_err("When writing the node dates table comment")?;
  }
  let mut csv = CsvWriter::new(&mut writer, format);
  let bounds: &[&str] = if with_confidence {
    &["lower bound", "upper bound"]
  } else {
    &[]
  };
  csv.write_record(["#node", "date", "numeric date"].iter().chain(bounds))?;
  tree.tree().preorder().iter().try_for_each(|&key| {
    let name = graph.names[key].to_owned();
    csv.write_record(once(name).chain(date_fields(dates, key, with_confidence)))
  })?;
  csv.into_inner()?;
  Ok(())
}

fn date_fields(dates: &TreeDates<'_>, key: GraphNodeKey, with_confidence: bool) -> Vec<String> {
  let numdate = dates.num_date[&key].filter(|_| !dates.excluded.contains(&key));
  let Some(numdate) = numdate else {
    let columns = if with_confidence { 4 } else { 2 };
    return vec![MISSING_DATE.to_owned(); columns];
  };
  let date = [year_fraction_to_datestring(numdate), numeric_date(numdate)];
  if !with_confidence {
    return date.to_vec();
  }
  let bounds = dates
    .confidence
    .and_then(|confidence| confidence.get(&key).copied())
    .map_or_else(
      || [MISSING_DATE.to_owned(), MISSING_DATE.to_owned()],
      |[lower, upper]| [numeric_date(lower), numeric_date(upper)],
    );
  date.into_iter().chain(bounds).collect()
}

fn numeric_date(value: f64) -> String {
  format!("{value:.6}")
}
