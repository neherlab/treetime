#[cfg(test)]
mod tests {
  use crate::ancestral::sample::{SampleMode, resolve_profile, sample_from_profile};
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rand::SeedableRng;
  use rand::rngs::StdRng;

  #[test]
  fn test_sample_deterministic_profile() {
    let profile = array![0.0, 0.0, 1.0, 0.0];
    let mut rng = StdRng::seed_from_u64(42);
    let idx = sample_from_profile(profile.view(), &mut rng);
    assert_eq!(2, idx);
  }

  #[test]
  fn test_sample_reproducible_with_seed() {
    let profile = array![0.25, 0.25, 0.25, 0.25];
    let mut rng1 = StdRng::seed_from_u64(123);
    let mut rng2 = StdRng::seed_from_u64(123);

    let results1: Vec<usize> = (0..100)
      .map(|_| sample_from_profile(profile.view(), &mut rng1))
      .collect();
    let results2: Vec<usize> = (0..100)
      .map(|_| sample_from_profile(profile.view(), &mut rng2))
      .collect();

    assert_eq!(results1, results2);
  }

  #[test]
  fn test_sample_respects_distribution() {
    let profile = array![0.9, 0.1, 0.0, 0.0];
    let mut rng = StdRng::seed_from_u64(42);

    let mut counts = [0_usize; 4];
    for _ in 0..1000 {
      counts[sample_from_profile(profile.view(), &mut rng)] += 1;
    }

    assert!(counts[0] > 800, "expected state 0 to dominate, got {}", counts[0]);
    assert!(counts[2] == 0, "expected state 2 never sampled, got {}", counts[2]);
    assert!(counts[3] == 0, "expected state 3 never sampled, got {}", counts[3]);
  }

  #[test]
  fn test_sample_zero_profile_returns_zero() {
    let profile = array![0.0, 0.0, 0.0, 0.0];
    let mut rng = StdRng::seed_from_u64(42);
    let idx = sample_from_profile(profile.view(), &mut rng);
    assert_eq!(0, idx);
  }

  #[test]
  fn test_sample_mode_default_is_argmax() {
    assert_eq!(SampleMode::Argmax, SampleMode::default());
  }

  #[test]
  fn test_resolve_profile_argmax_when_not_sampling() {
    let profile = array![0.1, 0.3, 0.6, 0.0];
    let mut rng = StdRng::seed_from_u64(42);
    let idx = resolve_profile(profile.view(), false, &mut rng);
    assert_eq!(2, idx);
  }

  #[test]
  fn test_resolve_profile_samples_when_true() {
    let profile = array![0.5, 0.5, 0.0, 0.0];
    let mut rng = StdRng::seed_from_u64(42);

    let mut saw_zero = false;
    let mut saw_one = false;
    for _ in 0..100 {
      match resolve_profile(profile.view(), true, &mut rng) {
        0 => saw_zero = true,
        1 => saw_one = true,
        _ => {},
      }
    }
    assert!(saw_zero && saw_one, "expected both states sampled from 50/50 profile");
  }
}
