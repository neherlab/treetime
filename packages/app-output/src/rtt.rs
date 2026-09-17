use eyre::Report;
use std::path::Path;
use treetime::clock::rtt::ClockRegressionResult;
use treetime_io::csv::CsvStructFileWriter;

/// Write the per-node root-to-tip regression results as delimited rows.
///
/// The core clock pipeline produces `ClockRegressionResult` values; this encoder projects them to the
/// `--output-clock-csv` file. Kept out of core so the core clock layer holds no output-format or
/// file-writing concern (PLAN W3.2).
pub fn write_clock_regression_result_csv(
  results: &[ClockRegressionResult],
  filepath: impl AsRef<Path>,
  delimiter: u8,
) -> Result<(), Report> {
  let mut rtt_writer = CsvStructFileWriter::new(filepath, delimiter)?;
  results.iter().try_for_each(|result| rtt_writer.write(result))
}
