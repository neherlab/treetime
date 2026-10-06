#[cfg(test)]
mod tests {
  use crate::progress::{LogLevel, LogSink, RunWarning, RunWarningKind};
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::warnings_enabled(  LogLevel::Warn,  vec![(LogLevel::Warn, o!("The tree repeats A."))])]
  #[case::warnings_disabled( LogLevel::Error, vec![])]
  #[trace]
  fn test_progress_default_warning_logs_the_message_at_warn_level(
    #[case] min_level: LogLevel,
    #[case] expected: Vec<(LogLevel, String)>,
  ) {
    let sink = helpers::RecordingSink {
      min_level,
      lines: Mutex::new(vec![]),
    };

    sink.warning(&RunWarning {
      kind: RunWarningKind::DuplicateNodeNames,
      message: o!("The tree repeats A."),
      names: vec![o!("A")],
    });

    assert_eq!(expected, sink.lines.into_inner());
  }

  mod helpers {
    use crate::progress::{LogLevel, LogSink};
    use parking_lot::Mutex;

    pub(super) struct RecordingSink {
      pub(super) min_level: LogLevel,
      pub(super) lines: Mutex<Vec<(LogLevel, String)>>,
    }

    impl LogSink for RecordingSink {
      fn log(&self, level: LogLevel, message: &str) {
        self.lines.lock().push((level, message.to_owned()));
      }

      fn log_enabled(&self, level: LogLevel) -> bool {
        level >= self.min_level
      }
    }
  }
}
