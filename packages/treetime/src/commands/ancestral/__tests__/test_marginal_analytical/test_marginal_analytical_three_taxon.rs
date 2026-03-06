#![cfg(test)]

use super::test_marginal_analytical_support::{
  analytical_three_taxon_likelihood, run_dense_marginal_get_log_lh, state_index,
};
use crate::gtr::get_gtr::{JC69Params, jc69};
use crate::pretty_assert_ulps_eq;
use eyre::Report;

/// Three-taxon tree with a single position, verified against closed-form Felsenstein
/// site likelihood with explicit summation over all internal node states.
#[test]
fn test_three_taxon_single_position_exhaustive() -> Result<(), Report> {
  let gtr = jc69(JC69Params::default())?;

  let t_a = 0.1;
  let t_b = 0.2;
  let t_ab = 0.15;
  let t_c = 0.25;

  let obs_a = state_index('A');
  let obs_b = state_index('C');
  let obs_c = state_index('G');

  let expected_lh = analytical_three_taxon_likelihood(&gtr, obs_a, obs_b, obs_c, t_a, t_b, t_ab, t_c);
  let expected_log_lh = expected_lh.ln();

  let newick = "((A:0.1,B:0.2)AB:0.15,C:0.25)root;";
  let aln = ">A\nA\n>B\nC\n>C\nG\n";
  let actual_log_lh = run_dense_marginal_get_log_lh(newick, aln, gtr)?;

  pretty_assert_ulps_eq!(expected_log_lh, actual_log_lh, epsilon = 1e-7);
  Ok(())
}

/// Exhaustive verification of all 64 single-position state combinations on a three-taxon tree under JC69.
#[test]
fn test_three_taxon_all_combinations() -> Result<(), Report> {
  let gtr = jc69(JC69Params::default())?;

  let t_a = 0.1;
  let t_b = 0.2;
  let t_ab = 0.15;
  let t_c = 0.25;

  let states = ['A', 'C', 'G', 'T'];
  let mut mismatch = None;

  'state_space: for &state_a in &states {
    for &state_b in &states {
      for &state_c in &states {
        let obs_a = state_index(state_a);
        let obs_b = state_index(state_b);
        let obs_c = state_index(state_c);

        let expected_lh = analytical_three_taxon_likelihood(&gtr, obs_a, obs_b, obs_c, t_a, t_b, t_ab, t_c);
        let expected_log_lh = expected_lh.ln();

        let newick = "((A:0.1,B:0.2)AB:0.15,C:0.25)root;";
        let aln = format!(">A\n{state_a}\n>B\n{state_b}\n>C\n{state_c}\n");
        let actual_log_lh = run_dense_marginal_get_log_lh(newick, &aln, gtr.clone())?;

        if (expected_log_lh - actual_log_lh).abs() >= 1e-7 {
          mismatch = Some((state_a, state_b, state_c, expected_log_lh, actual_log_lh));
          break 'state_space;
        }
      }
    }
  }

  if let Some((state_a, state_b, state_c, expected_log_lh, actual_log_lh)) = mismatch {
    panic!("Mismatch for triplet ({state_a},{state_b},{state_c}): expected {expected_log_lh}, got {actual_log_lh}");
  }

  Ok(())
}
