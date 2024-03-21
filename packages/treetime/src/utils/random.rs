use crate::make_internal_report;
use eyre::Report;
use ndarray_rand::rand::SeedableRng;
use rand::{seq::IteratorRandom, Rng};
use rand_isaac::Isaac64Rng;

pub fn get_random_number_generator(seed: Option<u64>) -> impl Rng {
  match seed {
    None => Isaac64Rng::from_entropy(),
    Some(seed) => Isaac64Rng::seed_from_u64(seed),
  }
}

pub fn random_choice_maybe<T>(iter: impl IntoIterator<Item = T>, rng: &mut impl Rng) -> Option<T> {
  iter.into_iter().choose(rng)
}

pub fn random_choice<T>(iter: impl IntoIterator<Item = T>, rng: &mut impl Rng) -> Result<T, Report> {
  random_choice_maybe(iter, rng)
    .ok_or_else(|| make_internal_report!("random_choice: expected at least one item, but none found"))
}
