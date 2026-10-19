#[cfg(test)]
mod tests {
  use super::super::test_fitch_pass_support::tests::helpers::{FitchInput, fitch_input, run_in_pool};
  use crate::partition::fitch::passes::{attach_seqs_to_graph, fitch_backward};
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::first_leaf_one_thread(  "A", 1)]
  #[case::first_leaf_four_threads("A", 4)]
  #[case::middle_leaf_one_thread( "C", 1)]
  #[case::middle_leaf_four_threads("C", 4)]
  #[case::last_leaf_one_thread(   "D", 1)]
  #[case::last_leaf_four_threads( "D", 4)]
  #[trace]
  fn test_fitch_pass_failure_of_one_leaf_returns_its_error(
    #[case] failing_leaf: &str,
    #[case] threads: usize,
  ) -> Result<(), Report> {
    let FitchInput { graph, names, partition, node_inputs } = fitch_input()?;
    let mut partition = attach_seqs_to_graph(&graph, partition, node_inputs)?;
    let leaf = find_node_key_by_name(&graph, &names, failing_leaf).expect("leaf exists");
    partition.nodes.remove(&leaf);

    let result = run_in_pool(threads, || fitch_backward(&graph, partition));

    assert_error!(
      result,
      format!(
        "Leaf {leaf} has no sequence before the Fitch backward pass. This is an internal error. Please report it to developers."
      )
    );
    Ok(())
  }
}
