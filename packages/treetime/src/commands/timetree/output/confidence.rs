use crate::timetree::confidence::NodeConfidenceInterval;
use eyre::Report;
use std::path::Path;
use treetime_io::csv::CsvStructFileWriter;

/// Write confidence intervals for node dates to TSV file.
pub fn write_confidence_intervals(intervals: &[NodeConfidenceInterval], filepath: &Path) -> Result<(), Report> {
  let mut writer = CsvStructFileWriter::new(filepath, b'\t')?;
  intervals.iter().try_for_each(|ci| writer.write(ci))
}
