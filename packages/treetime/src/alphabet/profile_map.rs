use crate::alphabet::alphabet::Alphabet;
use ndarray::Array2;

#[derive(Clone, Debug)]
pub struct ProfileMap {}

impl ProfileMap {
  pub fn from_alphabet(alphabet: &Alphabet) -> Self {
    // {s:x for s,x in zip(self.alphabet, np.eye(len(self.alphabet)))}
    let mtx = Array2::<f32>::eye(alphabet.len());

    Self {}
  }

  #[inline]
  pub fn len(&self) -> usize {
    0
  }

  #[inline]
  pub fn is_empty(&self) -> bool {
    self.len() == 0
  }
}
