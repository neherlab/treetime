use crate::make_internal_report;
use eyre::Report;
use rand::rngs::Xoshiro256PlusPlus;
use rand::seq::IteratorRandom;
use rand::{Rng, RngExt, SeedableRng};

pub fn get_random_number_generator(seed: u64) -> impl Rng + Send + Sync + Clone {
  Xoshiro256PlusPlus::seed_from_u64(seed)
}

pub fn random_choice<T>(iter: impl IntoIterator<Item = T>, rng: &mut impl Rng) -> Result<T, Report> {
  random_choice_maybe(iter, rng)
    .ok_or_else(|| make_internal_report!("random_choice: expected at least one item, but none found"))
}

pub fn random_choice_maybe<T>(iter: impl IntoIterator<Item = T>, rng: &mut impl Rng) -> Option<T> {
  iter.into_iter().choose(rng)
}

pub fn random_remove<T>(v: &mut Vec<T>, rng: &mut impl Rng) -> T {
  let index: usize = rng.random_range(0..v.len());
  v.remove(index)
}
