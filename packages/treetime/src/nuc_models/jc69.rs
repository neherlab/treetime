use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::profile_map::ProfileMap;
use crate::gtr::gtr::{GTRParams, GTR};
use eyre::Report;
use ndarray::{Array1, Array2};

pub fn jc69() -> Result<GTR, Report> {
  let mu = 1.0;
  let alphabet_name = "nuc";

  let alphabet = Alphabet::new(alphabet_name)?;
  let profile_map = ProfileMap::from_alphabet(&alphabet);

  let num_chars = alphabet.len();
  let W = Array2::<f32>::ones((num_chars, num_chars));
  let pi = Array1::<f32>::ones(num_chars);

  let gtr = GTR::new(&GTRParams {
    alphabet,
    profile_map,
    mu,
    W,
    pi,
  })?;

  Ok(gtr)
}
