#[cfg(test)]
pub(super) mod tests {
  use crate::ancestral::__tests__::prop_generators::input::MarginalTestInput;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, SparseReconstruction};
  use crate::test_utils::{run_dense_marginal, run_sparse_marginal};
  use eyre::Report;
  use treetime_io::nwk::{NwkParse, nwk_read};

  pub(crate) fn run_dense_marginal_for_input(input: &MarginalTestInput) -> Result<(f64, DenseReconstruction), Report> {
    let NwkParse {
      graph,
      branch_lengths,
      names,
      ..
    } = nwk_read(input.newick.as_bytes())?;
    run_dense_marginal(&graph, &branch_lengths, &names, &input.alignment, input.gtr.clone())
  }

  pub(crate) fn run_sparse_marginal_for_input(
    input: &MarginalTestInput,
  ) -> Result<(f64, SparseReconstruction), Report> {
    let NwkParse {
      graph,
      branch_lengths,
      names,
      ..
    } = nwk_read(input.newick.as_bytes())?;
    run_sparse_marginal(&graph, &branch_lengths, &names, &input.alignment, input.gtr.clone())
  }
}
