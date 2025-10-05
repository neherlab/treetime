use crate::commands::timetree::data::clock_model::ClockModel;
use eyre::Report;
use std::path::Path;

/// Write clock model parameters to file.
///
/// Core: Documents inferred molecular clock rate.
/// Why: Clock rate is key parameter for interpreting divergence times.
/// How: JSON file with rate, intercept, R², confidence intervals.
pub fn write_clock_model(_clock_model: &ClockModel, _out_base: &Path) -> Result<(), Report> {
  todo!("Write clock model to JSON file")
}
