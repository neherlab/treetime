use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::gtr::gtr::{GTR, GTRParams};
use crate::gtr::infer_gtr::{infer_gtr_dense, infer_gtr_sparse};
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::payload::ancestral::GraphAncestral;
use crate::{make_error, make_report};
use clap::ValueEnum;
use eyre::{Report, WrapErr};
use ndarray::{Array1, Array2, array};
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::sync::Arc;
use strum_macros::Display;

#[derive(
  Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ValueEnum, SmartDefault, Display, Serialize, Deserialize,
)]
pub enum GtrModelName {
  /// Infer GTR parameters from data.
  ///
  /// Dense mode requires pre-populated profiles, causing two reconstruction passes.
  #[default]
  Infer,
  JC69,
  K80,
  F81,
  HKY85,
  T92,
  TN93,
  #[value(name = "jtt92")]
  Jtt92,
}

pub fn get_gtr_sparse(
  name: &GtrModelName,
  partition: &Arc<RwLock<PartitionMarginalSparse>>,
  graph: &GraphAncestral,
) -> Result<GTR, Report> {
  match name {
    GtrModelName::Infer => infer_gtr_sparse(partition, graph),
    _ => get_gtr_by_name(*name),
  }
  .wrap_err_with(|| make_report!("When creating model '{name}'"))
}

/// Get GTR model for dense representation.
pub fn get_gtr_dense(
  name: &GtrModelName,
  partition: &Arc<RwLock<PartitionMarginalDense>>,
  graph: &GraphAncestral,
) -> Result<GTR, Report> {
  match name {
    GtrModelName::Infer => infer_gtr_dense(partition, graph),
    _ => get_gtr_by_name(*name),
  }
  .wrap_err_with(|| make_report!("When creating model '{name}'"))
}

fn get_gtr_by_name(name: GtrModelName) -> Result<GTR, Report> {
  match name {
    GtrModelName::Infer => make_error!("Cannot get GTR by name for 'Infer'"),
    GtrModelName::JC69 => jc69(JC69Params::default()),
    GtrModelName::F81 => f81(F81Params::default()),
    GtrModelName::HKY85 => hky85(HKY85Params::default()),
    GtrModelName::K80 => k80(K80Params::default()),
    GtrModelName::T92 => t92(T92Params::default()),
    GtrModelName::TN93 => tn93(TN93Params::default()),
    GtrModelName::Jtt92 => jtt92(Jtt92Params::default()),
  }
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
/// See: Jukes and Cantor (1969). Evolution of Protein Molecules. New York: Academic Press. pp. 21-132
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
/// See: Kimura (1980),  J. Mol. Evol. 16 (2): 111-120. doi:10.1007/BF01731581.
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
/// See: Felsenstein (1981), J. Mol. Evol. 17  (6): 368-376. doi:10.1007/BF01734359
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
/// See: Hasegawa, Kishino, Yano (1985), J. Mol. Evol. 22 (2): 160-174. doi:10.1007/BF02101694
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
/// See: Tamura K (1992),  Mol.  Biol. Evol. 9 (4): 678-687.  DOI: 10.1093/oxfordjournals.molbev.a040752
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

#[derive(Copy, Clone, Debug, SmartDefault)]
pub struct Jtt92Params {
  #[default = 1.0]
  pub mu: f64,

  #[default(AlphabetName::AaNoStop)]
  pub alphabet: AlphabetName,

  pub treat_gap_as_unknown: bool,
}

/// Jones-Taylor-Thornton 1992 amino acid substitution model.
///
/// Empirical model derived from protein sequence alignments. Uses a 20x20 rate matrix
/// for the 20 standard amino acids.
///
/// See: Jones, Taylor, Thornton (1992). CABIOS 8(3):275-282
#[allow(clippy::excessive_precision)]
pub fn jtt92(
  Jtt92Params {
    mu,
    alphabet,
    treat_gap_as_unknown,
  }: Jtt92Params,
) -> Result<GTR, Report> {
  let alphabet = Alphabet::new(alphabet, treat_gap_as_unknown)?;

  // Stationary frequencies (pi) from JTT92 model
  // Order: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
  #[rustfmt::skip]
  let pi = array![
    0.07674789, 0.05169087, 0.04264509, 0.05154407, 0.01980301,
    0.04075195, 0.06182989, 0.07315199, 0.02294399, 0.05376110,
    0.09190390, 0.05867583, 0.02382594, 0.04012589, 0.05090097,
    0.06876503, 0.05856501, 0.01426057, 0.03210196, 0.06600504
  ];

  // Rate matrix Q from JTT92 model (20x20)
  #[rustfmt::skip]
  let q = array![
    [-1.247831, 0.044229, 0.041179, 0.061769, 0.042704, 0.043467, 0.08007, 0.136501, 0.02059, 0.027453, 0.022877, 0.02669, 0.041179, 0.011439, 0.14794, 0.288253, 0.362223, 0.006863, 0.008388, 0.227247],
    [0.029789, -1.025965, 0.023112, 0.008218, 0.058038, 0.159218, 0.014895, 0.070364, 0.168463, 0.011299, 0.019517, 0.33179, 0.022599, 0.002568, 0.038007, 0.051874, 0.032871, 0.064714, 0.010272, 0.008731],
    [0.022881, 0.019068, -1.280568, 0.223727, 0.014407, 0.03644, 0.024576, 0.034322, 0.165676, 0.019915, 0.005085, 0.11144, 0.012712, 0.004237, 0.006356, 0.213134, 0.098304, 0.00339, 0.029661, 0.00678],
    [0.041484, 0.008194, 0.270413, -1.044903, 0.005121, 0.025095, 0.392816, 0.066579, 0.05736, 0.005634, 0.003585, 0.013316, 0.007682, 0.002049, 0.007682, 0.030217, 0.019462, 0.002049, 0.023559, 0.015877],
    [0.011019, 0.022234, 0.00669, 0.001968, -0.56571, 0.001771, 0.000984, 0.011609, 0.013577, 0.003345, 0.004526, 0.001377, 0.0061, 0.015348, 0.002755, 0.043878, 0.008264, 0.022628, 0.041124, 0.012199],
    [0.02308, 0.125524, 0.034823, 0.019841, 0.003644, -1.04415, 0.130788, 0.010528, 0.241735, 0.003644, 0.029154, 0.118235, 0.017411, 0.00162, 0.066406, 0.021461, 0.020651, 0.007288, 0.009718, 0.008098],
    [0.064507, 0.017816, 0.035632, 0.471205, 0.003072, 0.198435, -0.944343, 0.073107, 0.015973, 0.007372, 0.005529, 0.111197, 0.011058, 0.003072, 0.011058, 0.01843, 0.019659, 0.006143, 0.0043, 0.027646],
    [0.130105, 0.099578, 0.058874, 0.09449, 0.042884, 0.018898, 0.086495, -0.647831, 0.016717, 0.004361, 0.004361, 0.019625, 0.010176, 0.003634, 0.017444, 0.146096, 0.023986, 0.039976, 0.005815, 0.034162],
    [0.006155, 0.074775, 0.089138, 0.025533, 0.01573, 0.1361, 0.005927, 0.005243, -1.135695, 0.003648, 0.012767, 0.010259, 0.007523, 0.009119, 0.026217, 0.016642, 0.010487, 0.001824, 0.130629, 0.002508],
    [0.01923, 0.011752, 0.025106, 0.005876, 0.009081, 0.004808, 0.00641, 0.003205, 0.008547, -1.273602, 0.122326, 0.011218, 0.25587, 0.047542, 0.005342, 0.021367, 0.130873, 0.004808, 0.017094, 0.513342],
    [0.027395, 0.0347, 0.010958, 0.006392, 0.021003, 0.065748, 0.008219, 0.005479, 0.051137, 0.209115, -0.668139, 0.012784, 0.354309, 0.226465, 0.093143, 0.053877, 0.022829, 0.047485, 0.021916, 0.16437],
    [0.020405, 0.376625, 0.153332, 0.015158, 0.004081, 0.170239, 0.105525, 0.015741, 0.026235, 0.012243, 0.008162, -0.900734, 0.037896, 0.002332, 0.012243, 0.027401, 0.06005, 0.00583, 0.004664, 0.008162],
    [0.012784, 0.010416, 0.007102, 0.003551, 0.007339, 0.01018, 0.004261, 0.003314, 0.007812, 0.113397, 0.091854, 0.015388, -1.182051, 0.01018, 0.003788, 0.006865, 0.053503, 0.005682, 0.004261, 0.076466],
    [0.00598, 0.001993, 0.003987, 0.001595, 0.031098, 0.001595, 0.001993, 0.001993, 0.015948, 0.035484, 0.098877, 0.001595, 0.017144, -0.637182, 0.006778, 0.03668, 0.004784, 0.021131, 0.213701, 0.024719],
    [0.098117, 0.037426, 0.007586, 0.007586, 0.007081, 0.082944, 0.009104, 0.012138, 0.058162, 0.005058, 0.051587, 0.010621, 0.008092, 0.008598, -0.727675, 0.144141, 0.059679, 0.003035, 0.005058, 0.011632],
    [0.258271, 0.069009, 0.343678, 0.040312, 0.152366, 0.036213, 0.020498, 0.137334, 0.049878, 0.02733, 0.040312, 0.032113, 0.019814, 0.06286, 0.194728, -1.447863, 0.325913, 0.023914, 0.043045, 0.025964],
    [0.276406, 0.037242, 0.135003, 0.022112, 0.02444, 0.029677, 0.018621, 0.019203, 0.026768, 0.142567, 0.014548, 0.059936, 0.131511, 0.006983, 0.068665, 0.27757, -1.335389, 0.006983, 0.01222, 0.065174],
    [0.001275, 0.017854, 0.001134, 0.000567, 0.016295, 0.002551, 0.001417, 0.007793, 0.001134, 0.001275, 0.007368, 0.001417, 0.003401, 0.00751, 0.00085, 0.004959, 0.0017, -0.312785, 0.010061, 0.003542],
    [0.003509, 0.006379, 0.022328, 0.014673, 0.066664, 0.007655, 0.002233, 0.002552, 0.182769, 0.010207, 0.007655, 0.002552, 0.005741, 0.170967, 0.00319, 0.020095, 0.006698, 0.022647, -0.605978, 0.005103],
    [0.195438, 0.011149, 0.010493, 0.020331, 0.040662, 0.013117, 0.029512, 0.030824, 0.007214, 0.630254, 0.11805, 0.009182, 0.211834, 0.040662, 0.015084, 0.024922, 0.073453, 0.016396, 0.010493, -1.241722]
  ];

  // Compute W from Q: W_ij = Q_ij * sqrt(pi_j / pi_i)
  // This is the symmetric exchangeability matrix
  let n = 20;
  let mut w = Array2::<f64>::zeros((n, n));
  for i in 0..n {
    for j in 0..n {
      if i != j {
        let ratio: f64 = pi[j] / pi[i];
        w[[i, j]] = q[[i, j]] * ratio.sqrt();
      }
    }
  }

  GTR::new(GTRParams {
    alphabet,
    mu,
    W: Some(w),
    pi,
  })
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

#[derive(Clone, Debug, SmartDefault)]
pub struct TN93Params {
  /// Substitution rate
  #[default = 1.0]
  pub mu: f64,

  /// Transversion rate (A<->C, A<->T, G<->C, G<->T) relative to A<->G = 1
  #[default = 1.0]
  pub kappa1: f64,

  /// Pyrimidine transition rate (C<->T) relative to A<->G = 1
  #[default = 1.0]
  pub kappa2: f64,

  /// Equilibrium frequencies [piA, piC, piG, piT]
  #[default(None)]
  pub pi: Option<Array1<f64>>,

  #[default(AlphabetName::Nuc)]
  pub alphabet: AlphabetName,

  pub treat_gap_as_unknown: bool,
}

/// Tamura-Nei 1993 model.
///
/// Distinguishes between the two types of transitions: A<->G has rate 1 (reference),
/// C<->T has rate kappa2. All transversions have rate kappa1. Equilibrium frequencies
/// can be non-uniform.
///
/// See: Tamura, Nei (1993), Mol Biol Evol. 10(3): 512-526. DOI: 10.1093/oxfordjournals.molbev.a040023
pub fn tn93(
  TN93Params {
    mu,
    kappa1,
    kappa2,
    pi,
    alphabet,
    treat_gap_as_unknown,
  }: TN93Params,
) -> Result<GTR, Report> {
  let alphabet = Alphabet::new(alphabet, treat_gap_as_unknown)?;
  let num_chars = alphabet.n_canonical();

  // W matrix for alphabet [A, C, G, T]:
  // - A<->G (purine transition): rate = 1 (reference)
  // - C<->T (pyrimidine transition): rate = kappa2
  // - All transversions: rate = kappa1
  #[rustfmt::skip]
  let W = Some(array![
    [1.0,    kappa1, 1.0,    kappa1],  // A
    [kappa1, 1.0,    kappa1, kappa2],  // C
    [1.0,    kappa1, 1.0,    kappa1],  // G
    [kappa1, kappa2, kappa1, 1.0   ]   // T
  ]);

  let pi = pi.unwrap_or_else(|| Array1::ones(num_chars) / (num_chars as f64));
  let pi = &pi / pi.sum();

  GTR::new(GTRParams { alphabet, mu, W, pi })
}
