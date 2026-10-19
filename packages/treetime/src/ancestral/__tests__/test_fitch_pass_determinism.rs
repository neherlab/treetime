#[cfg(test)]
mod tests {
  use super::super::test_fitch_pass_support::tests::helpers::{FitchInput, fitch_input, run_in_pool};
  use crate::partition::fitch::partition::PartitionFitch;
  use crate::partition::fitch::passes::compress_sequences;
  use eyre::Report;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_fitch_pass_determinism_results_are_equal_across_thread_counts_and_runs() -> Result<(), Report> {
    let compress = |threads| -> Result<PartitionFitch, Report> {
      let FitchInput {
        graph,
        partition,
        node_inputs,
        ..
      } = fitch_input()?;
      run_in_pool(threads, || compress_sequences(&graph, partition, node_inputs))
    };

    let expected = compress(1)?;

    assert_eq!(expected, compress(1)?);
    assert_eq!(expected, compress(4)?);
    assert_eq!(expected, compress(4)?);
    Ok(())
  }
}
