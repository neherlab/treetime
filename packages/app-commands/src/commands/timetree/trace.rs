use app_output::timetree_trace::TraceCsvSink;
use eyre::Report;
use std::path::Path;
use treetime::progress::StageSink;
use treetime::timetree::convergence::metrics::IterationRecord;
use treetime::timetree::convergence::optimizer::TraceSink;

pub(crate) struct TimetreeTraceSink<'a> {
  stages: &'a dyn StageSink,
  csv: Option<TraceCsvSink>,
}

impl<'a> TimetreeTraceSink<'a> {
  pub(crate) fn new(tracelog: Option<&Path>, stages: &'a dyn StageSink) -> Result<Self, Report> {
    Ok(Self {
      stages,
      csv: tracelog.map(TraceCsvSink::new).transpose()?,
    })
  }

  pub(crate) fn finish(self) -> Result<(), Report> {
    self.csv.map_or(Ok(()), TraceCsvSink::finish)
  }
}

impl TraceSink for TimetreeTraceSink<'_> {
  fn emit(&mut self, record: &IterationRecord) -> Result<(), Report> {
    self.stages.iteration(record);
    self.csv.as_mut().map_or(Ok(()), |csv| csv.emit(record))
  }
}
