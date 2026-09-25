use app_output::timetree_trace::TraceCsvSink;
use eyre::Report;
use std::path::Path;
use treetime::progress::ProgressSink;
use treetime::timetree::convergence::metrics::IterationRecord;
use treetime::timetree::convergence::optimizer::TraceSink;
use treetime_utils::io::file::create_file_or_stdout;

pub(crate) fn timetree_trace_sink<'a>(
  tracelog: Option<&Path>,
  progress: &'a dyn ProgressSink,
) -> Result<Box<dyn TraceSink + 'a>, Report> {
  let mut sinks: Vec<Box<dyn TraceSink + 'a>> = vec![Box::new(ProgressTraceSink { progress })];
  if let Some(path) = tracelog {
    sinks.push(Box::new(TraceCsvSink::new(create_file_or_stdout(path)?)?));
  }
  Ok(Box::new(TraceSinks { sinks }))
}

struct ProgressTraceSink<'a> {
  progress: &'a dyn ProgressSink,
}

impl TraceSink for ProgressTraceSink<'_> {
  fn emit(&mut self, record: &IterationRecord) -> Result<(), Report> {
    self.progress.iteration(record);
    Ok(())
  }
}

struct TraceSinks<'a> {
  sinks: Vec<Box<dyn TraceSink + 'a>>,
}

impl TraceSink for TraceSinks<'_> {
  fn emit(&mut self, record: &IterationRecord) -> Result<(), Report> {
    self.sinks.iter_mut().try_for_each(|sink| sink.emit(record))
  }
}
