#[cfg(test)]
mod tests {
  use super::super::test_marginal_analytical_support::tests::{analytical_star_tree_likelihood, state_index};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::GTR;
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::run_dense_marginal_with_newick;
  use eyre::Report;
  use helpers::nonuniform_pi_gtr;
  use itertools::Itertools;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::jc69_all_same(     jc69(JC69Params::default()).unwrap(), 0.1,  "AAAA")]
  #[case::jc69_mixed_states( jc69(JC69Params::default()).unwrap(), 0.2,  "ACGT")]
  #[case::nonuniform_pi(     nonuniform_pi_gtr(),                  0.15, "TTGC")]
  #[trace]
  fn test_star_tree_log_likelihood_matches_analytical(
    #[case] gtr: GTR,
    #[case] t: f64,
    #[case] tip_states: &str,
  ) -> Result<(), Report> {
    let observations = tip_states.chars().map(state_index).collect_vec();
    let expected_log_lh = analytical_star_tree_likelihood(&gtr, &observations, t).ln();

    let newick = format!("(A:{t},B:{t},C:{t},D:{t})root;");
    let aln = ["A", "B", "C", "D"]
      .iter()
      .zip(tip_states.chars())
      .map(|(name, state)| format!(">{name}\n{state}\n"))
      .join("");
    let (actual_log_lh, _) = run_dense_marginal_with_newick(&newick, &aln, &gtr)?;

    pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::gtr::gtr::{GTR, Mu};
    use ndarray::array;

    pub(super) fn nonuniform_pi_gtr() -> GTR {
      let n_states = Alphabet::new(AlphabetName::Nuc).unwrap().n_canonical();
      GTR::builder()
        .n_states(n_states)
        .mu(Mu::FromMatrix(1.0))
        .pi(array![0.1, 0.2, 0.3, 0.4])
        .build()
        .unwrap()
    }
  }
}
