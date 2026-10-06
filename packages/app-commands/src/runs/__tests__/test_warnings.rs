#[cfg(test)]
mod tests {
  use crate::runs::warnings::WarningCollector;
  use pretty_assertions::assert_eq;
  use treetime::o;
  use treetime::progress::{LogSink, NoopProgress, RunWarning, RunWarningKind};

  #[test]
  fn test_warnings_collector_keeps_warnings_when_the_inner_sink_logs_nothing() {
    let collector = WarningCollector::new(&NoopProgress);
    let warning = RunWarning {
      kind: RunWarningKind::DuplicateSequenceNames,
      message: o!("The alignment repeats A."),
      names: vec![o!("A")],
    };

    collector.warning(&warning);

    assert_eq!(vec![warning], collector.into_warnings());
  }
}
