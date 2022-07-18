use ndarray_rand::rand::SeedableRng;
use rand::Rng;
use rand_isaac::Isaac64Rng;

pub fn get_random_number_generator(seed: Option<u64>) -> impl Rng {
  match seed {
    None => Isaac64Rng::from_entropy(),
    Some(seed) => Isaac64Rng::seed_from_u64(seed),
  }
}
