use crate::make_internal_report;
use eyre::Report;
use itertools::Itertools;
use ndarray_rand::rand::SeedableRng;
use rand::{seq::IteratorRandom, seq::SliceRandom, Rng};
use rand_isaac::Isaac64Rng;
use std::collections::BTreeSet;

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

pub fn random_pop<T: Clone + Ord>(s: &mut BTreeSet<T>, rng: &mut impl Rng) -> T {
  // TODO: inefficient. Try to avoid copying. Sets have no indexed access, so removing random element is tricky.
  // Might be possible with IndexSet?
  let mut v = s.iter().cloned().collect_vec();
  let item = random_remove(&mut v, rng);
  *s = BTreeSet::new();
  s.extend(v);
  item
}

pub fn random_sequence(length: usize, seed: Option<u64>) -> Vec<char> {
  let mut rng = get_random_number_generator(seed);
  let letters = ['A', 'C', 'G', 'T', 'N', '-'];
  letters.choose_multiple(&mut rng, length).copied().collect()
}
