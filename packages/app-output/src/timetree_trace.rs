use eyre::Report;
use std::path::Path;
use treetime::timetree::convergence::metrics::IterationRecord;
use treetime::timetree::convergence::optimizer::TraceSink;
use treetime_io::csv::CsvStructFileWriter;

pub struct TraceCsvSink {
  writer: CsvStructFileWriter,
}

impl TraceCsvSink {
  pub fn new(filepath: impl AsRef<Path>) -> Result<Self, Report> {
    Ok(Self {
      writer: CsvStructFileWriter::new(filepath, b',')?,
    })
  }

  pub fn finish(self) -> Result<(), Report> {
    self.writer.finish()
  }
}

impl TraceSink for TraceCsvSink {
  fn emit(&mut self, record: &IterationRecord) -> Result<(), Report> {
    self.writer.write(&record.metrics)
  }
}
