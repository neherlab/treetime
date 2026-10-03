use app_output::timetree_trace::TraceCsvSink;
use eyre::Report;
use std::path::Path;
use treetime::progress::StageSink;
use treetime::timetree::convergence::metrics::IterationRecord;
use treetime::timetree::convergence::optimizer::TraceSink;
use treetime_utils::io::file::create_file_or_stdout;

pub(crate) fn timetree_trace_sink<'a>(
  tracelog: Option<&Path>,
  stages: &'a dyn StageSink,
) -> Result<Box<dyn TraceSink + 'a>, Report> {
  let mut sinks: Vec<Box<dyn TraceSink + 'a>> = vec![Box::new(StageTraceSink { stages })];
  if let Some(path) = tracelog {
    sinks.push(Box::new(TraceCsvSink::new(Box::new(create_file_or_stdout(path)?))?));
  }
  Ok(Box::new(TraceSinks { sinks }))
}

struct StageTraceSink<'a> {
  stages: &'a dyn StageSink,
}

impl TraceSink for StageTraceSink<'_> {
  fn emit(&mut self, record: &IterationRecord) -> Result<(), Report> {
    self.stages.iteration(record);
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
