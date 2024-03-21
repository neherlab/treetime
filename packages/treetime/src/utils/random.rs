use crate::make_internal_report;
use eyre::Report;
use ndarray_rand::rand::SeedableRng;
use rand::{seq::IteratorRandom, Rng};
use rand_isaac::Isaac64Rng;

pub fn get_random_number_generator(seed: Option<u64>) -> (impl Rng + Send + Sync + Clone) {
  match seed {
    None => Isaac64Rng::from_entropy(),
    Some(seed) => Isaac64Rng::seed_from_u64(seed),
  }
}

pub fn clone_random_number_generator(rng: &mut impl Rng) -> (impl Rng + Send + Sync + Clone) {
  Isaac64Rng::from_rng(rng).expect("Unable to clone random number generator")
}

pub fn random_choice_maybe<T>(iter: impl IntoIterator<Item = T>, rng: &mut impl Rng) -> Option<T> {
  iter.into_iter().choose(rng)
}

pub fn random_choice<T>(iter: impl IntoIterator<Item = T>, rng: &mut impl Rng) -> Result<T, Report> {
  random_choice_maybe(iter, rng)
    .ok_or_else(|| make_internal_report!("random_choice: expected at least one item, but none found"))
}

pub fn random_remove<T>(v: &mut Vec<T>, rng: &mut impl Rng) -> T {
  let index: usize = rng.gen_range(0..v.len());
  v.remove(index)
}
