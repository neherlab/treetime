use crate::json_float::JsonFloat;
use crate::results::run_results::{CommandResults, results_of_record};
use crate::results::year_date::YearDate;
use crate::runs::record::RunRecord;
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_with::skip_serializing_none;
use std::path::Path;

pub fn run_headline(record: &RunRecord, out_dir: &Path) -> Result<RunHeadline, Report> {
  let results = results_of_record(record, out_dir)?;
  Ok(match results.results {
    CommandResults::Timetree(results) => results
      .estimates
      .map(|estimates| RunHeadline {
        root_date: estimates.root_date,
        clock_rate: estimates.clock_rate.map(JsonFloat),
        r_squared: estimates.r_squared.map(JsonFloat),
        ..RunHeadline::default()
      })
      .unwrap_or_default(),
    CommandResults::Clock(results) => RunHeadline {
      root_date: None,
      clock_rate: results.estimates.clock_rate.map(JsonFloat),
      r_squared: results.estimates.r_squared.map(JsonFloat),
      ..RunHeadline::default()
    },
    CommandResults::Homoplasy(results) => RunHeadline {
      recurrent_substitutions: results.statistics.map(|statistics| statistics.recurrent_substitutions),
      ..RunHeadline::default()
    },
    CommandResults::Ancestral(_)
    | CommandResults::Mugration(_)
    | CommandResults::Optimize(_)
    | CommandResults::Prune(_) => RunHeadline::default(),
  })
}

/// Key results of a finished run, for run lists; the same values the run's results show.
#[skip_serializing_none]
#[derive(
  Clone, Debug, Default, PartialEq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize,
)]
#[deser(skip_serializing_optionals)]
pub struct RunHeadline {
  /// Date of the root of a time tree.
  pub root_date: Option<YearDate>,
  /// Clock rate in substitutions per site per year.
  pub clock_rate: Option<JsonFloat>,
  /// Coefficient of determination of the clock model.
  pub r_squared: Option<JsonFloat>,
  /// Number of distinct substitutions on two or more branches of a homoplasy run.
  pub recurrent_substitutions: Option<usize>,
}
