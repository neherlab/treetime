use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::gtr::gtr::{GTRParams, GTR};
use crate::{make_error, make_report};
use clap::ArgEnum;
use eyre::{Report, WrapErr};
use ndarray::{array, Array1, Array2};
use smart_default::SmartDefault;
use strum_macros::Display;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, SmartDefault, Display)]
pub enum GtrModelName {
  #[default]
  JC69,
  K80,
  F81,
  HKY85,
  T92,
  Infer,
}

pub fn get_gtr(name: &GtrModelName, alphabet: &Option<AlphabetName>) -> Result<GTR, Report> {
  match name {
    GtrModelName::Infer => {
      unimplemented!("Not implemented: GTR inference is not yet implemented. Please provide the name of a concrete model with `--gtr=<model>` argument.")
    }
    GtrModelName::JC69 => jc69(JC69Params::default()),
    GtrModelName::F81 => f81(F81Params::default()),
    GtrModelName::HKY85 => hky85(HKY85Params::default()),
    GtrModelName::K80 => k80(K80Params::default()),
    GtrModelName::T92 => t92(T92Params::default()),
  }
  .wrap_err_with(|| {
    let alphabet_msg = match alphabet {
      None => "".to_owned(),
      Some(alphabet) => format!("with alphabet '{alphabet}'"),
    };
    make_report!("When creating model '{name}'{alphabet_msg}")
  })
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


//Trimmed down python version.

// def infer_gtr(nij, Ti, root_state, fixed_pi=None, alphabet, pc=1.0):
// r"""
// Infer a GTR model by specifying the number of transitions and time spent in each
// character. The basic equation that is being solved is

// :math:`n_{ij} = pi_i W_{ij} T_j`

// where :math:`n_{ij}` are the transitions, :math:`pi_i` are the equilibrium
// state frequencies, :math:`W_{ij}` is the "substitution attempt matrix",
// while :math:`T_i` is the time on the tree spent in character state
// :math:`i`. To regularize the process, we add pseudocounts and also need
// to account for the fact that the root of the tree is in a particular
// state. the modified equation is

// :math:`n_{ij} + pc = pi_i W_{ij} (T_j+pc+root\_state)`

// Parameters
// ----------

//  nij : nxn matrix
//     The number of times a change in character state is observed
//     between state j and i

//  Ti :n vector
//     The time spent in each character state

//  root_state : n vector
//     The number of characters in state i in the sequence
//     of the root node.

//  pc : float
//     Pseudocounts, this determines the lower cutoff on the rate when
//     no substitutions are observed
// fixed_pi : n vector of None.
//
// alphabet : Alphabet

// """
// //assert that size of n_ij is nxn and size of Ti is n, where n is the alphabet size
// dp = 1e-5 // convergence criterion -- could be made an argument
// Nit = 40  // maximum number of iterations -- could be made an argument
// pc_mat = pc*np.ones_like(nij)
// np.fill_diagonal(pc_mat, 0.0)
// np.fill_diagonal(nij, 0.0)
// pi_old = np.zeros_like(Ti)
// if fixed_pi is None:
//     pi = np.ones_like(Ti)
// else:
//     pi = np.copy(fixed_pi)
// pi/=pi.sum()
// W_ij = np.ones_like(nij)
// mu = (nij.sum()+pc)/(Ti.sum()+pc) // initial guess for the rate
// # if pi is fixed, this will immediately converge
// iteration_counter = 0
// while distance(pi_old, pi) > dp and interaction_counter < Nit:  // distance is sum([(x_i - y_i)**2 for x_i, y_i in zip(pi_old, pi)])**0.5
//     iteration_counter += 1
//     pi_old = np.copy(pi)
//     W_ij = (nij+nij.T+2*pc_mat)/mu/(np.outer(pi,Ti) + np.outer(Ti,pi) + 2*pc_mat) // np.outer(x,y) of two vectors x and y is a matrix where the i,j element is x_i*y_j

//     np.fill_diagonal(W_ij, 0)
//     scale_factor = avg_transition(W_ij,pi) // this function already exists

//     W_ij = W_ij/scale_factor
//     if fixed_pi is None:
//         pi = (np.sum(nij+pc_mat,axis=1)+root_state)/(mu*np.dot(W_ij,Ti)+root_state.sum()+np.sum(pc_mat, axis=1))
//         pi /= pi.sum()
//         mu = (nij.sum() + pc)/(np.sum(pi * (W_ij.dot(Ti)))+pc)
//     else:
//         mu = (nij.sum() + pc)/(np.sum(pi * (W_ij.dot(pi)))*Ti.sum() + pc)

// if count >= Nit:
//     if distance(pi_old,pi) > dp:
//         gtr.logger('the iterative scheme has not converged',3,warn=True)
//     elif np.abs(1-np.max(pi.sum(axis=0))) > dp:
//         gtr.logger('the iterative scheme has converged, but proper normalization was not reached',3,warn=True)

// gtr.assign_rates(mu=mu, W=W_ij, pi=pi)
// return gtr
