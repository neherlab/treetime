use crate::alphabet::alphabet::Alphabet;
use eyre::Report;
use indexmap::IndexMap;
use itertools::Itertools;
use lazy_static::lazy_static;
use ndarray::array;
use ndarray::Array1;

lazy_static! {
  static ref PROFILE_MAP_NUC: IndexMap<char, Array1<f64>> = IndexMap::<char, Array1<f64>>::from([
    ('A', array![1.0, 0.0, 0.0, 0.0, 0.0]),
    ('C', array![0.0, 1.0, 0.0, 0.0, 0.0]),
    ('G', array![0.0, 0.0, 1.0, 0.0, 0.0]),
    ('T', array![0.0, 0.0, 0.0, 1.0, 0.0]),
    ('-', array![0.0, 0.0, 0.0, 0.0, 1.0]),
    ('N', array![1.0, 1.0, 1.0, 1.0, 1.0]),
    ('X', array![1.0, 1.0, 1.0, 1.0, 1.0]),
    ('R', array![1.0, 0.0, 1.0, 0.0, 0.0]),
    ('Y', array![0.0, 1.0, 0.0, 1.0, 0.0]),
    ('S', array![0.0, 1.0, 1.0, 0.0, 0.0]),
    ('W', array![1.0, 0.0, 0.0, 1.0, 0.0]),
    ('K', array![0.0, 0.0, 1.0, 1.0, 0.0]),
    ('M', array![1.0, 1.0, 0.0, 0.0, 0.0]),
    ('D', array![1.0, 0.0, 1.0, 1.0, 0.0]),
    ('H', array![1.0, 1.0, 0.0, 1.0, 0.0]),
    ('B', array![0.0, 1.0, 1.0, 1.0, 0.0]),
    ('V', array![1.0, 1.0, 1.0, 0.0, 0.0]),
  ]);
}

#[derive(Clone, Debug)]
pub struct ProfileMap {
  profile_map: IndexMap<char, Array1<f64>>,
}

impl ProfileMap {
  pub fn new(name: &str) -> Result<Self, Report> {
    // TODO: implement remaining profile maps as well as generation from alphabet
    let profile_map = match name {
      "nuc" => PROFILE_MAP_NUC.to_owned(),
      _ => unimplemented!(),
    };

    Ok(Self { profile_map })
  }

  pub fn from_alphabet(alphabet: &Alphabet) -> Result<Self, Report> {
    Self::new(&alphabet.name)
  }

  #[inline]
  pub fn get(&self, c: char) -> &Array1<f64> {
    self
      .profile_map
      .get(&c)
      .ok_or_else(|| {
        format!(
          "When accessing profile map: Unknown character: '{c}'. Known characters: {}",
          self.profile_map.keys().join(", ")
        )
      })
      .unwrap()
  }
}

//
// profile_maps = {
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
//
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
// }
