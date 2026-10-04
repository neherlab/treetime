use crate::output_plan::OutputSelection;
use crate::table_output::table_create;
use eyre::Report;
use std::path::Path;
use treetime::timetree::convergence::metrics::IterationRecord;
use treetime::timetree::convergence::optimizer::TraceSink;
use treetime_io::csv::CsvWriter;
use treetime_utils::io::file::FileWriter;

pub struct TraceCsvSink {
  writer: CsvWriter<FileWriter>,
}

impl TraceCsvSink {
  pub fn new(filepath: &Path) -> Result<Self, Report> {
    Ok(Self {
      writer: table_create(OutputSelection::Tracelog, filepath)?,
    })
  }

  pub fn finish(self) -> Result<(), Report> {
    self.writer.finish()
  }
}

impl TraceSink for TraceCsvSink {
  fn emit(&mut self, record: &IterationRecord) -> Result<(), Report> {
    self.writer.write_row(&record.metrics)
  }
}
