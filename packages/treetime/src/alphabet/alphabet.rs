use ndarray::Array1;
use num_traits::NumCast;
use std::ops::Index;

#[derive(Clone, Debug)]
pub struct Alphabet {}

const DUMMY: Array1<f32> = Array1::<f32>::from_iter([1.0, 0.0, 0.0, 0.0, 0.0]);

impl<I: NumCast> Index<I> for Alphabet {
  type Output = Array1<f32>;

  fn index(&self, index: I) -> &Self::Output {
    let index = index.to_usize().unwrap();
    &DUMMY
  }
}

impl Alphabet {
  pub fn new(name: &str) -> Self {
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

// const alphabets = {
// 'nuc':{
//     'A': np.array([1, 0, 0, 0, 0], dtype='float'),
//     'C': np.array([0, 1, 0, 0, 0], dtype='float'),
//     'G': np.array([0, 0, 1, 0, 0], dtype='float'),
//     'T': np.array([0, 0, 0, 1, 0], dtype='float'),
//     '-': np.array([0, 0, 0, 0, 1], dtype='float'),
//     'N': np.array([1, 1, 1, 1, 1], dtype='float'),
//     'X': np.array([1, 1, 1, 1, 1], dtype='float'),
//     'R': np.array([1, 0, 1, 0, 0], dtype='float'),
//     'Y': np.array([0, 1, 0, 1, 0], dtype='float'),
//     'S': np.array([0, 1, 1, 0, 0], dtype='float'),
//     'W': np.array([1, 0, 0, 1, 0], dtype='float'),
//     'K': np.array([0, 0, 1, 1, 0], dtype='float'),
//     'M': np.array([1, 1, 0, 0, 0], dtype='float'),
//     'D': np.array([1, 0, 1, 1, 0], dtype='float'),
//     'H': np.array([1, 1, 0, 1, 0], dtype='float'),
//     'B': np.array([0, 1, 1, 1, 0], dtype='float'),
//     'V': np.array([1, 1, 1, 0, 0], dtype='float')
//     },

// 'nuc_nogap':{
//     'A': np.array([1, 0, 0, 0], dtype='float'),
//     'C': np.array([0, 1, 0, 0], dtype='float'),
//     'G': np.array([0, 0, 1, 0], dtype='float'),
//     'T': np.array([0, 0, 0, 1], dtype='float'),
//     '-': np.array([1, 1, 1, 1], dtype='float'), # gaps are completely ignored in distance computations
//     'N': np.array([1, 1, 1, 1], dtype='float'),
//     'X': np.array([1, 1, 1, 1], dtype='float'),
//     'R': np.array([1, 0, 1, 0], dtype='float'),
//     'Y': np.array([0, 1, 0, 1], dtype='float'),
//     'S': np.array([0, 1, 1, 0], dtype='float'),
//     'W': np.array([1, 0, 0, 1], dtype='float'),
//     'K': np.array([0, 0, 1, 1], dtype='float'),
//     'M': np.array([1, 1, 0, 0], dtype='float'),
//     'D': np.array([1, 0, 1, 1], dtype='float'),
//     'H': np.array([1, 1, 0, 1], dtype='float'),
//     'B': np.array([0, 1, 1, 1], dtype='float'),
//     'V': np.array([1, 1, 1, 0], dtype='float')
//     },
//
// 'aa':{
//     'A': np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Alanine         Ala
//     'C': np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Cysteine        Cys
//     'D': np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Aspartic AciD   Asp
//     'E': np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamic Acid   Glu
//     'F': np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Phenylalanine   Phe
//     'G': np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glycine         Gly
//     'H': np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Histidine       His
//     'I': np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Isoleucine      Ile
//     'K': np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Lysine          Lys
//     'L': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Leucine         Leu
//     'M': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Methionine      Met
//     'N': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #AsparagiNe      Asn
//     'P': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Proline         Pro
//     'Q': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamine       Gln
//     'R': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #ARginine        Arg
//     'S': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype='float'), #Serine          Ser
//     'T': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], dtype='float'), #Threonine       Thr
//     'V': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype='float'), #Valine          Val
//     'W': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], dtype='float'), #Tryptophan      Trp
//     'Y': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], dtype='float'), #Tyrosine        Tyr
//     '*': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], dtype='float'), #stop
//     '-': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], dtype='float'), #gap
//     'X': np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype='float'), #not specified/any
//     'B': np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Asparagine/Aspartic Acid    Asx
//     'Z': np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamine/Glutamic Acid     Glx
//     },
//
// 'aa_nogap':{
//     'A': np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Alanine         Ala
//     'C': np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Cysteine        Cys
//     'D': np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Aspartic AciD   Asp
//     'E': np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamic Acid   Glu
//     'F': np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Phenylalanine   Phe
//     'G': np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Glycine         Gly
//     'H': np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Histidine       His
//     'I': np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Isoleucine      Ile
//     'K': np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Lysine          Lys
//     'L': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Leucine         Leu
//     'M': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Methionine      Met
//     'N': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #AsparagiNe      Asn
//     'P': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Proline         Pro
//     'Q': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamine       Gln
//     'R': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], dtype='float'), #ARginine        Arg
//     'S': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype='float'), #Serine          Ser
//     'T': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], dtype='float'), #Threonine       Thr
//     'V': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], dtype='float'), #Valine          Val
//     'W': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], dtype='float'), #Tryptophan      Trp
//     'Y': np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], dtype='float'), #Tyrosine        Tyr
//     'X': np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype='float'), #not specified/any
//     'B': np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], dtype='float'), #Asparagine/Aspartic Acid    Asx
//     'Z': np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], dtype='float'), #Glutamine/Glutamic Acid     Glx
//     }
// };
