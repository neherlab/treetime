#[cfg(test)]
mod tests {
  use crate::gtr::gtr_site_specific::GTRSiteSpecific;
  use ndarray::prelude::*;
  use rand::SeedableRng;
  use rand::rngs::StdRng;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_gtr_site_specific_random_zero_shapes_give_uniform_model() {
    let gtr = GTRSiteSpecific::random()
      .n_states(4)
      .seq_len(3)
      .avg_mu(0.5)
      .pi_dirichlet_alpha(0.0)
      .W_dirichlet_alpha(0.0)
      .mu_gamma_alpha(0.0)
      .rng(&mut StdRng::seed_from_u64(0))
      .call()
      .unwrap();

    let expected_pi = Array2::from_elem((4, 3), 0.25);
    let expected_W = Array2::from_shape_fn((4, 4), |(i, j)| if i == j { 0.0 } else { 4.0 / 3.0 });
    let expected_mu = array![0.5, 0.5, 0.5];

    pretty_assert_abs_diff_eq!(expected_pi, gtr.pi, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(expected_W, gtr.W, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(expected_mu, gtr.mu, epsilon = 1e-15);
  }
}
