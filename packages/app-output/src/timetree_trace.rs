use eyre::Report;
use std::io::Write;
use treetime::timetree::convergence::metrics::ConvergenceMetrics;
use treetime::timetree::convergence::optimizer::TraceSink;
use treetime_io::csv::CsvStructWriter;

pub struct TraceCsvSink {
  writer: CsvStructWriter<Box<dyn Write + Send>>,
}

impl TraceCsvSink {
  pub fn new(writer: Box<dyn Write + Send>) -> Result<Self, Report> {
    Ok(Self {
      writer: CsvStructWriter::new(writer, b',')?,
    })
  }
}

impl TraceSink for TraceCsvSink {
  fn emit(&mut self, metric: &ConvergenceMetrics) -> Result<(), Report> {
    self.writer.write(metric)
  }
}
