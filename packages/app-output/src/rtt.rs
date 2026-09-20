use eyre::Report;
use std::path::Path;
use treetime::clock::rtt::ClockRegressionResult;
use treetime_io::csv::CsvStructFileWriter;

pub fn write_clock_regression_result_csv(
  results: &[ClockRegressionResult],
  filepath: impl AsRef<Path>,
  delimiter: u8,
) -> Result<(), Report> {
  let mut rtt_writer = CsvStructFileWriter::new(filepath, delimiter)?;
  results.iter().try_for_each(|result| rtt_writer.write(result))
}
