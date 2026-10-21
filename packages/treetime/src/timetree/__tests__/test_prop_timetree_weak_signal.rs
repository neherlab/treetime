#[cfg(test)]
mod tests {
  use crate::cancel::NoopCancel;
  use crate::error::OperationError;
  use crate::progress::NoopProgress;
  use crate::test_utils::timetree_cell_params;
  use crate::timetree::params::TimetreeParams;
  use crate::timetree::pipeline;
  use generators::arb_weak_signal_case;
  use proptest::prelude::{ProptestConfig, TestCaseError, prop_assert, proptest};

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(8))]

    #[ignore = "kb/issues/M-timetree-weak-signal-errors-differ-by-code-path.md"]
    #[test]
    fn test_prop_timetree_weak_signal_succeeds_or_rejects_the_input(spec in arb_weak_signal_case()) {
      let case = spec.simulate().map_err(|error| TestCaseError::fail(format!("{error:?}")))?;
      let input = case.timetree_cell_input(spec.length).map_err(|error| TestCaseError::fail(format!("{error:?}")))?;
      let params = TimetreeParams {
        coalescent_opt: true,
        ..timetree_cell_params()
      };

      let result = pipeline::run(&params, input, None, None, &NoopCancel, &NoopProgress, &NoopProgress);

      prop_assert!(
        matches!(result, Ok(_) | Err(OperationError::InvalidInput(_))),
        "Timetree run of the case with seed {} ended with {:?}",
        spec.seed,
        result.err()
      );
    }
  }

  mod generators {
    use crate::test_utils::{CaseSpec, timetree_cell};
    use proptest::arbitrary::any;
    use proptest::sample::select;
    use proptest::strategy::Strategy;

    pub(super) fn arb_weak_signal_case() -> impl Strategy<Value = CaseSpec> {
      (any::<u64>(), select(vec![1e-4, 3e-4]), 200..=500_usize).prop_map(|(seed, mu, length)| CaseSpec {
        mu,
        length,
        ..timetree_cell(seed).expect("the standard timetree cell is valid")
      })
    }
  }
}
