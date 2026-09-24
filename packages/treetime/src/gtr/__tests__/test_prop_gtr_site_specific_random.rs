#[cfg(test)]
mod tests {
  use crate::gtr::gtr_site_specific::GTRSiteSpecific;
  use ndarray::prelude::*;
  use proptest::prelude::*;
  use rand::SeedableRng;
  use rand::rngs::StdRng;
  use treetime_utils::{prop_assert_array_abs_diff_eq, prop_assert_relative_eq};

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn test_prop_gtr_site_specific_random_pi_columns_normalized(
      (n_states, seq_len, seed) in generators::arb_shape_and_seed(),
      params in generators::arb_random_params(),
    ) {
      let gtr = helpers::random_gtr(n_states, seq_len, seed, params);
      prop_assert_array_abs_diff_eq!(gtr.pi.sum_axis(Axis(0)), Array1::ones(seq_len), epsilon = 1e-12);
    }

    #[test]
    fn test_prop_gtr_site_specific_random_w_symmetric_zero_diagonal(
      (n_states, seq_len, seed) in generators::arb_shape_and_seed(),
      params in generators::arb_random_params(),
    ) {
      let gtr = helpers::random_gtr(n_states, seq_len, seed, params);
      prop_assert_array_abs_diff_eq!(gtr.W.t(), gtr.W.view(), epsilon = 1e-15);
      prop_assert_array_abs_diff_eq!(gtr.W.diag().to_owned(), Array1::zeros(n_states), epsilon = 1e-15);
    }

    #[test]
    fn test_prop_gtr_site_specific_random_mean_rate_matches_avg_mu(
      (n_states, seq_len, seed) in generators::arb_shape_and_seed(),
      params in generators::arb_random_params(),
    ) {
      let gtr = helpers::random_gtr(n_states, seq_len, seed, params);
      let mean_rate = gtr.average_rate().mean().unwrap();
      prop_assert_relative_eq!(params.avg_mu, mean_rate, max_relative = 1e-12);
    }

    #[test]
    fn test_prop_gtr_site_specific_random_same_seed_same_model(
      (n_states, seq_len, seed) in generators::arb_shape_and_seed(),
      params in generators::arb_random_params(),
    ) {
      let first = helpers::random_gtr(n_states, seq_len, seed, params);
      let second = helpers::random_gtr(n_states, seq_len, seed, params);
      prop_assert_eq!(&first.mu, &second.mu);
      prop_assert_eq!(&first.W, &second.W);
      prop_assert_eq!(&first.pi, &second.pi);
    }
  }

  mod generators {
    use super::helpers::RandomParams;
    use proptest::prelude::*;

    pub(super) fn arb_shape_and_seed() -> impl Strategy<Value = (usize, usize, u64)> {
      (2_usize..=20, 1_usize..=8, any::<u64>())
    }

    pub(super) fn arb_random_params() -> impl Strategy<Value = RandomParams> {
      (0.01_f64..10.0, 0.5_f64..5.0, 0.5_f64..5.0, 0.5_f64..5.0).prop_map(
        |(avg_mu, pi_dirichlet_alpha, W_dirichlet_alpha, mu_gamma_alpha)| RandomParams {
          avg_mu,
          pi_dirichlet_alpha,
          W_dirichlet_alpha,
          mu_gamma_alpha,
        },
      )
    }
  }

  mod helpers {
    use super::*;

    #[derive(Clone, Copy, Debug)]
    pub(super) struct RandomParams {
      pub(crate) avg_mu: f64,
      pub(crate) pi_dirichlet_alpha: f64,
      pub(crate) W_dirichlet_alpha: f64,
      pub(crate) mu_gamma_alpha: f64,
    }

    pub(super) fn random_gtr(n_states: usize, seq_len: usize, seed: u64, params: RandomParams) -> GTRSiteSpecific {
      GTRSiteSpecific::random()
        .n_states(n_states)
        .seq_len(seq_len)
        .avg_mu(params.avg_mu)
        .pi_dirichlet_alpha(params.pi_dirichlet_alpha)
        .W_dirichlet_alpha(params.W_dirichlet_alpha)
        .mu_gamma_alpha(params.mu_gamma_alpha)
        .rng(&mut StdRng::seed_from_u64(seed))
        .call()
        .unwrap()
    }
  }
}
