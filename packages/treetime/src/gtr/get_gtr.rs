use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::gtr::gtr::{GTRParams, GTR};
use crate::port::infer_gtr::{get_mutation_counts, infer_gtr, InferGtrOptions, InferGtrResult};
use crate::representation::graph_dense::DenseGraph;
use crate::representation::graph_sparse::SparseGraph;
use crate::{make_error, make_report};
use clap::ArgEnum;
use eyre::{Report, WrapErr};
use ndarray::{array, Array1, Array2};
use smart_default::SmartDefault;
use strum_macros::Display;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, SmartDefault, Display)]
pub enum GtrModelName {
  #[default]
  Infer,
  JC69,
  K80,
  F81,
  HKY85,
  T92,
}

pub fn get_gtr(name: &GtrModelName, alphabet: &Alphabet, graph: &SparseGraph) -> Result<GTR, Report> {
  match name {
    GtrModelName::Infer => {
      let counts = get_mutation_counts(graph, alphabet)?;
      let InferGtrResult { W, pi, mu } = infer_gtr(&counts, &InferGtrOptions::default())?;
      let alphabet = alphabet.to_owned();
      let W = Some(W);
      GTR::new(GTRParams { alphabet, mu, W, pi })
    }
    GtrModelName::JC69 => jc69(JC69Params::default()),
    GtrModelName::F81 => f81(F81Params::default()),
    GtrModelName::HKY85 => hky85(HKY85Params::default()),
    GtrModelName::K80 => k80(K80Params::default()),
    GtrModelName::T92 => t92(T92Params::default()),
  }
  .wrap_err_with(|| make_report!("When creating model '{name}'"))
}

pub fn get_gtr_dense(name: &GtrModelName, _alphabet: &Alphabet, _graph: &DenseGraph) -> Result<GTR, Report> {
  match name {
    GtrModelName::Infer => {
      unimplemented!("Model inference is not yet implemented for dense representation. Please set model explicitly.")
    }
    GtrModelName::JC69 => jc69(JC69Params::default()),
    GtrModelName::F81 => f81(F81Params::default()),
    GtrModelName::HKY85 => hky85(HKY85Params::default()),
    GtrModelName::K80 => k80(K80Params::default()),
    GtrModelName::T92 => t92(T92Params::default()),
  }
  .wrap_err_with(|| make_report!("When creating model '{name}'"))
}

#[derive(Copy, Clone, Debug, SmartDefault)]
pub struct JC69Params {
  /// Substitution rate
  #[default = 1.0]
  pub mu: f64,

  #[default(AlphabetName::Nuc)]
  pub alphabet: AlphabetName,

  pub treat_gap_as_unknown: bool,
}

/// Jukes-Cantor 1969 model.
///
/// This model assumes equal concentrations of the nucleotides and equal transition rates
/// between nucleotide states.
///
/// See: Jukes and Cantor (1969). Evolution of Protein Molecules. New York: Academic Press. pp. 21–132
pub fn jc69(
  JC69Params {
    mu,
    alphabet,
    treat_gap_as_unknown,
  }: JC69Params,
) -> Result<GTR, Report> {
  let alphabet = Alphabet::new(alphabet, treat_gap_as_unknown)?;
  let num_chars = alphabet.n_canonical();
  let W = Some(Array2::<f64>::ones((num_chars, num_chars)));
  let pi = Array1::<f64>::ones(num_chars);
  GTR::new(GTRParams { alphabet, mu, W, pi })
}

#[derive(Copy, Clone, Debug, SmartDefault)]
pub struct K80Params {
  /// Substitution rate
  #[default = 1.0]
  pub mu: f64,

  /// Ratio of transversion/transition rates
  #[default = 0.1]
  pub kappa: f64,

  #[default(AlphabetName::Nuc)]
  pub alphabet: AlphabetName,

  pub treat_gap_as_unknown: bool,
}

/// Kimura 1980 model.
///
/// Assumes equal concentrations across nucleotides, but
/// allows different rates between transitions and transversions. The ratio
/// of the transversion/transition rates is given by kappa parameter.
///
/// NOTE: Current implementation of the model does not account for the gaps.
///
/// See: Kimura (1980),  J. Mol. Evol. 16 (2): 111–120. doi:10.1007/BF01731581.
pub fn k80(
  K80Params {
    mu,
    kappa,
    alphabet,
    treat_gap_as_unknown,
  }: K80Params,
) -> Result<GTR, Report> {
  let alphabet = Alphabet::new(alphabet, treat_gap_as_unknown)?;
  let num_chars = alphabet.n_canonical();
  let W = Some(create_transversion_transition_W(&alphabet, kappa)?);
  let pi = Array1::<f64>::ones(num_chars) / (num_chars as f64);
  GTR::new(GTRParams { alphabet, mu, W, pi })
}

#[derive(Copy, Clone, Debug, SmartDefault)]
pub struct F81Params {
  /// Substitution rate
  #[default = 1.0]
  pub mu: f64,

  #[default(AlphabetName::Nuc)]
  pub alphabet: AlphabetName,

  pub treat_gap_as_unknown: bool,
}

/// Felsenstein 1981 model.
///
/// Assumes non-equal concentrations across nucleotides,
/// but the transition rate between all states is assumed to be equal.
///
/// See: Felsenstein (1981), J. Mol. Evol. 17  (6): 368–376. doi:10.1007/BF01734359
pub fn f81(
  F81Params {
    mu,
    alphabet,
    treat_gap_as_unknown,
  }: F81Params,
) -> Result<GTR, Report> {
  let alphabet = Alphabet::new(alphabet, treat_gap_as_unknown)?;
  let num_chars = alphabet.n_canonical();
  let W = Some(Array2::<f64>::ones((num_chars, num_chars)));
  let pi: Array1<f64> = {
    let pi = Array1::<f64>::ones(num_chars) / (num_chars as f64);
    let sum = pi.sum();
    pi / sum
  };
  GTR::new(GTRParams { alphabet, mu, W, pi })
}

#[derive(Copy, Clone, Debug, SmartDefault)]
pub struct HKY85Params {
  /// Substitution rate
  #[default = 1.0]
  pub mu: f64,

  /// Ratio of transversion/transition rates
  #[default = 0.1]
  pub kappa: f64,

  #[default(AlphabetName::Nuc)]
  pub alphabet: AlphabetName,

  pub treat_gap_as_unknown: bool,
}

/// Hasegawa, Kishino and Yano 1985 model.
///
/// Allows different concentrations of the nucleotides (as in F81) and distinguishes between transition/transversion
/// substitutions (similar to K80).
///
/// NOTE: Current implementation of the model does not account for the gaps
///
/// See: Hasegawa, Kishino, Yano (1985), J. Mol. Evol. 22 (2): 160–174. doi:10.1007/BF02101694
pub fn hky85(
  HKY85Params {
    mu,
    kappa,
    alphabet,
    treat_gap_as_unknown,
  }: HKY85Params,
) -> Result<GTR, Report> {
  let alphabet = Alphabet::new(alphabet, treat_gap_as_unknown)?;
  let num_chars = alphabet.n_canonical();
  let W = Some(create_transversion_transition_W(&alphabet, kappa)?);
  let pi: Array1<f64> = {
    let pi = Array1::<f64>::ones(num_chars) / (num_chars as f64);
    let sum = pi.sum();
    pi / sum
  };
  GTR::new(GTRParams { alphabet, mu, W, pi })
}

#[derive(Copy, Clone, Debug, SmartDefault)]
pub struct T92Params {
  /// Substitution rate
  #[default = 1.0]
  pub mu: f64,

  /// Ratio of transversion/transition rates
  #[default = 0.1]
  pub kappa: f64,

  /// Relative GC content
  #[default = 0.5]
  pub pi_GC: f64,

  #[default(AlphabetName::Nuc)]
  pub alphabet: AlphabetName,

  pub treat_gap_as_unknown: bool,
}

/// Tamura 1992 model.
///
/// Extending Kimura (1980) model for the case where a G+C-content bias exists.
///
/// NOTE: Current implementation of the model does not account for the gaps
///
/// See: Tamura K (1992),  Mol.  Biol. Evol. 9 (4): 678–687.  DOI: 10.1093/oxfordjournals.molbev.a040752
pub fn t92(
  T92Params {
    mu,
    kappa,
    pi_GC,
    alphabet,
    treat_gap_as_unknown,
  }: T92Params,
) -> Result<GTR, Report> {
  if !(0.0..=1.0).contains(&pi_GC) {
    return make_error!("The relative GC should be between 0 and 1, but found pi_GC={pi_GC}");
  }

  let alphabet = Alphabet::new(alphabet, treat_gap_as_unknown)?;
  let W = Some(create_transversion_transition_W(&alphabet, kappa)?);
  let pi = array![(1.0 - pi_GC) * 0.5, pi_GC * 0.5, pi_GC * 0.5, (1.0 - pi_GC) * 0.5];
  GTR::new(GTRParams { alphabet, mu, W, pi })
}

fn create_transversion_transition_W(alphabet: &Alphabet, kappa: f64) -> Result<Array2<f64>, Report> {
  let num_chars = alphabet.n_canonical();
  let mut W = Array2::<f64>::ones((num_chars, num_chars));
  W[[0, 2]] = kappa;
  W[[1, 3]] = kappa;
  W[[2, 0]] = kappa;
  W[[3, 1]] = kappa;
  Ok(W)
}
