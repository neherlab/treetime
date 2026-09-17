use eyre::Report;
use std::io::Write;
use treetime::timetree::convergence::metrics::ConvergenceMetrics;
use treetime::timetree::convergence::optimizer::TraceSink;
use treetime_io::csv::CsvStructWriter;

/// Encodes timetree optimizer convergence metrics as CSV.
///
/// One comma-delimited row per optimization iteration, with the columns of `ConvergenceMetrics`. The
/// core optimizer emits metrics to this sink; the CSV format lives here rather than in the core.
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
