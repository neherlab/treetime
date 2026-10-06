#[cfg(test)]
mod tests {
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::given_seed_without_random_step(Some(7), None,                       (7, vec![]))]
  #[case::given_seed_with_random_step(   Some(7), Some("Polytomy resolution"), (7, vec!["Polytomy resolution is stochastic; seed 7 (pass --seed to reproduce this run)"]))]
  #[trace]
  fn test_seed_args_given_seed(
    #[case] seed: Option<u64>,
    #[case] random_step: Option<&str>,
    #[case] expected: (u64, Vec<&str>),
  ) {
    let (actual_seed, messages) = helpers::resolve(seed, random_step);
    let expected = (expected.0, expected.1.into_iter().map(str::to_owned).collect::<Vec<_>>());
    assert_eq!(expected, (actual_seed, messages));
  }

  #[test]
  fn test_seed_args_drawn_seed_is_logged_for_a_random_step() {
    let (seed, messages) = helpers::resolve(None, Some("Sampling from the profile"));
    let expected = vec![format!(
      "Sampling from the profile is stochastic; seed {seed} (pass --seed to reproduce this run)"
    )];
    assert_eq!(expected, messages);
  }

  #[test]
  fn test_seed_args_drawn_seed_is_silent_without_a_random_step() {
    let (_, messages) = helpers::resolve(None, None);
    assert_eq!(Vec::<String>::new(), messages);
  }

  mod helpers {
    use crate::commands::shared::seed::SeedArgs;
    use crate::job::{JobEvent, JobProgress};
    use parking_lot::Mutex;
    use pretty_assertions::assert_eq;
    use treetime::progress::{LogEvent, LogLevel};

    pub(super) fn resolve(seed: Option<u64>, random_step: Option<&str>) -> (u64, Vec<String>) {
      let events = Mutex::new(vec![]);
      let log = JobProgress::new(|event: JobEvent| {
        if let JobEvent::Log {
          data: LogEvent { level, message },
        } = event
        {
          assert_eq!(LogLevel::Info, level);
          events.lock().push(message);
        }
      });
      let seed = SeedArgs { seed }.resolve(random_step, &log);
      (seed, events.into_inner())
    }
  }
}
