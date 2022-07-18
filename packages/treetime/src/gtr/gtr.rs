#![allow(clippy::type_complexity)]
#![allow(non_snake_case)]

use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::profile_map::ProfileMap;
use crate::make_report;
use crate::utils::ndarray::{clamp_min, outer};
use eyre::Report;
use ndarray::prelude::*;
use ndarray_einsum_beta::einsum;
use ndarray_linalg::UPLO::Lower;
use ndarray_linalg::{Eig, Eigh};
use num_traits::abs;
use num_traits::real::Real;

pub fn avg_transition(W: &Array2<f64>, pi: &Array1<f64>, gap_index: Option<usize>) -> Result<f64, Report> {
  let result = einsum("i,ij,j", &[pi, W, pi])
    .map_err(|err| make_report!("einsum: {err}"))?
    .into_dimensionality::<Ix0>()?
    .into_scalar();

  Ok(if let Some(gap_index) = gap_index {
    // np.sum(pi*W[:,gap_index]) *pi[gap_index])/(1-pi[gap_index])
    let W_slice = W.slice(s!(.., gap_index));
    let pi_mul_W_slice = pi * &W_slice;
    let pi_mul_W_slice_sum = pi_mul_W_slice.sum();
    (result - pi_mul_W_slice_sum * pi[gap_index]) / (1.0 - pi[gap_index])
  } else {
    result
  })
}

/// Performs eigendecomposition of the rate matrix and stores the left- and right-
/// matrices to convert the sequence profiles to the GTR matrix eigenspace
/// and hence to speed-up the computations.
/// NOTE: this assumes the diagonal of W is all zeros
fn eig_single_site(W: &Array2<f64>, pi: &Array1<f64>) -> Result<(Array1<f64>, Array2<f64>, Array2<f64>), Report> {
  assert!(abs(W.diag().sum()) < 1e-10);

  let sqrt_pi: Array1<f64> = pi.mapv(f64::sqrt);
  let mut sym_Q: Array2<f64> = W * outer(&sqrt_pi, &sqrt_pi)?;

  let diag = -(W * pi).sum_axis(Axis(1));
  sym_Q.diag_mut().assign(&diag);

  let (eigvals, eigvecs) = sym_Q.eigh(Lower)?;

  let tmp_v: Array2<f64> = eigvecs.t().to_owned() * sqrt_pi.to_owned();
  let one_norm: Array1<f64> = tmp_v.mapv(f64::abs).sum_axis(Axis(1));

  let v = tmp_v.t().to_owned() / &one_norm;
  let v_inv = (eigvecs * one_norm).t().to_owned() / sqrt_pi;

  Ok((eigvals, v, v_inv))
}

#[derive(Clone, Debug)]
pub struct GTRParams {
  pub alphabet: Alphabet,
  pub profile_map: ProfileMap,
  pub mu: f64,
  pub W: Array2<f64>,
  pub pi: Array1<f64>,
}

/// Defines General-Time-Reversible model of character evolution.
#[derive(Clone, Debug)]
pub struct GTR {
  pub debug: bool,
  pub is_site_specific: bool,
  pub alphabet: Alphabet,
  pub profile_map: ProfileMap,
  pub mu: f64,
  pub W: Array2<f64>,
  pub pi: Array1<f64>,
  pub eigvals: Array1<f64>,
  pub v: Array2<f64>,
  pub v_inv: Array2<f64>,
}

impl GTR {
  pub fn new(
    GTRParams {
      alphabet,
      profile_map,
      mu,
      W,
      pi,
    }: &GTRParams,
  ) -> Result<Self, Report> {
    assert!(!alphabet.is_empty(), "Alphabet should not be empty");
    assert_eq!(
      pi.shape().to_vec(),
      [alphabet.len()],
      "Length of equilibrium frequency vector (`pi`) does not match the alphabet length"
    );
    assert_eq!(
      W.shape().to_vec(),
      [alphabet.len(), alphabet.len()],
      "Dimensions of substitution matrix (`W`) don't match the alphabet size"
    );

    // self.state_index= {s:si for si,s in enumerate(self.alphabet)}
    // self.state_index.update({s:si for si,s in enumerate(self.alphabet)})

    // self.ambiguous = None
    // let gap_index = Self::assign_gap_and_ambiguous(alphabet);

    let pi = {
      let pi_sum = pi.sum();
      pi / pi_sum
    };

    let average_rate = avg_transition(W, &pi, alphabet.gap_index())?;

    let W = {
      let mut W = 0.5 * (&W.view() + &W.t());
      W.diag_mut().fill(0.0);
      W / average_rate
    };

    let mu = mu * average_rate;

    let (eigvals, v, v_inv) = eig_single_site(&W, &pi)?;

    Ok(Self {
      debug: false,
      is_site_specific: false,
      alphabet: alphabet.clone(),
      profile_map: profile_map.clone(),
      mu,
      W,
      pi,
      eigvals,
      v,
      v_inv,
    })
  }

  #[allow(clippy::missing_const_for_fn)]
  fn assign_gap_and_ambiguous(alphabet: &Alphabet) -> Option<usize> {
    // let n_states = self.alphabet.len();

    // // determine if a character exists that corresponds to no info, i.e. all one profile
    // if any([x.sum()==n_states for x in self.profile_map.values()]):
    //     amb_states = [c for c,x in self.profile_map.items() if x.sum()==n_states]
    //     self.ambiguous = 'N' if 'N' in amb_states else amb_states[0]
    // else:
    //     self.ambiguous=None
    //
    // // check for a gap symbol
    // try:
    //     self.gap_index = self.state_index['-']
    // except:
    //     self.gap_index=None

    None
  }

  /// Compute the probability of the sequence state of the child
  /// at time t later, given the parent profile.
  ///
  /// Parameters
  /// ----------
  ///
  ///  profile : numpy.array
  ///     Sequence profile. Shape = (L, a),
  ///     where L - sequence length, a - alphabet size.
  ///
  ///  t : double
  ///     Time to propagate
  ///
  ///  return_log: bool
  ///     If true, return log-probability
  ///
  /// Returns
  /// -------
  ///
  ///  res : np.array
  ///     Profile of the sequence after time t in the future.
  ///     Shape = (L, a), where L - sequence length, a - alphabet size.
  pub fn evolve(&self, profile: &Array2<f64>, t: f64, return_log: bool) -> Array2<f64> {
    let Qt = self.expQt(t);
    let res = profile.dot(&Qt.t());
    if return_log {
      res.mapv(f64::ln)
    } else {
      res
    }
  }

  /// Compute the probability of the sequence state of the parent
  /// at time (t+t0, backwards), given the sequence state of the
  /// child (profile) at time t0.
  ///
  /// Parameters
  /// ----------
  ///
  ///  profile : numpy.array
  ///     Sequence profile. Shape = (L, a),
  ///     where L - sequence length, a - alphabet size.
  ///
  ///  t : double
  ///     Time to propagate
  ///
  ///  return_log: bool
  ///     If True, return log-probability
  ///
  /// Returns
  /// -------
  ///
  ///  res : np.array
  ///     Profile of the sequence after time t in the past.
  ///     Shape = (L, a), where L - sequence length, a - alphabet size.
  pub fn propagate_profile(&self, profile: &Array2<f64>, t: f64, return_log: bool) -> Array2<f64> {
    let Qt = self.expQt(t);
    let res = profile.dot(&Qt);

    if return_log {
      res.mapv(f64::ln)
    } else {
      res
    }
  }

  /// Matrix exponential of exo(Qt)
  pub fn expQt(&self, t: f64) -> Array2<f64> {
    let eLambdaT: Array2<f64> = Array2::from_diag(&self.exp_lt(t)); // vector length = a

    let eLambdaT_dot_v_inv: Array2<f64> = eLambdaT.dot(&self.v_inv);

    let Qt: Array2<f64> = self.v.dot(&eLambdaT_dot_v_inv); // This is P(nuc1 | given nuc_2)

    clamp_min(&Qt, 0.0)
  }

  fn exp_lt(&self, t: f64) -> Array1<f64> {
    (self.mu * t * &self.eigvals).mapv(f64::exp)
  }
}

#[allow(clippy::excessive_precision)]
#[cfg(test)]
mod test {
  use super::*;
  use crate::nuc_models::jc69::jc69;
  use approx::{assert_abs_diff_eq, assert_ulps_eq};
  use eyre::Report;
  use ndarray::{array, Array1, Array2};
  use rstest::rstest;

  #[rstest]
  fn computes_eig_single_site() -> Result<(), Report> {
    let W: Array2<f64> = array![
      [0.00, 1.25, 1.25, 1.25, 1.25],
      [1.25, 0.00, 1.25, 1.25, 1.25],
      [1.25, 1.25, 0.00, 1.25, 1.25],
      [1.25, 1.25, 1.25, 0.00, 1.25],
      [1.25, 1.25, 1.25, 1.25, 0.00],
    ];

    let pi: Array1<f64> = array![0.2, 0.2, 0.2, 0.2, 0.2];

    let (eigvals, v, v_inv) = eig_single_site(&W, &pi)?;

    #[rustfmt::skip]
    assert_ulps_eq!(
      eigvals,
      array![ -1.2500000000000002e+00, -1.2500000000000000e+00, -1.2499999999999998e+00, -1.2499999999999996e+00, 1.6728544325582780e-30 ],
    );

    #[rustfmt::skip]
    assert_ulps_eq!(
      v,
      array![
        [  1.2561935999126067e-01,  0.0000000000000000e+00,  4.4842166780616644e-01,  1.3740991709998629e-01, 1.9999999999999998e-01 ],
        [  3.6748767484131067e-01,  0.0000000000000000e+00, -7.9161468922062522e-02, -3.0952717143492964e-01, 2.0000000000000001e-01 ],
        [ -5.0000000000000011e-01, -4.9452390820531411e-17,  5.1578332193833711e-02, -1.9047282856507031e-01, 1.9999999999999998e-01 ],
        [  3.4464825837142365e-03, -4.9999999999999989e-01, -2.1041926553896867e-01,  1.8129504145000691e-01, 1.9999999999999998e-01 ],
        [  3.4464825837142365e-03,  5.0000000000000000e-01, -2.1041926553896867e-01,  1.8129504145000691e-01, 1.9999999999999998e-01 ],
      ],
    );

    #[rustfmt::skip]
    assert_ulps_eq!(
      v_inv,
      array![
        [ 3.1338154657190392e-01,  9.1676836990647648e-01, -1.2473457379248944e+00,  8.5979107232566597e-03,  8.5979107232566597e-03 ],
        [ 0.0000000000000000e+00,  0.0000000000000000e+00, -9.8904781641062834e-17, -1.0000000000000000e+00,  1.0000000000000002e+00 ],
        [ 1.5019412552636964e+00, -2.6514302170767157e-01,  1.7275620372809566e-01, -7.0477721864205978e-01, -7.0477721864205978e-01 ],
        [ 6.3408966320044879e-01, -1.4283392642158992e+00, -8.7895294795774470e-01,  8.3660127448659782e-01,  8.3660127448659782e-01 ],
        [ 9.9999999999999933e-01,  9.9999999999999967e-01,  9.9999999999999933e-01,  9.9999999999999956e-01,  9.9999999999999956e-01 ],
      ],
    );

    Ok(())
  }

  #[rstest]
  fn avg_transition_test() -> Result<(), Report> {
    let pi = array![1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];

    let Wi: Array2<f64> = array![
      [0.0, 4.0 / 3.0, 4.0 / 3.0],
      [4.0 / 3.0, 0.0, 4.0 / 3.0],
      [4.0 / 3.0, 4.0 / 3.0, 0.0],
    ];
    // test without gap index
    assert_ulps_eq!(avg_transition(&Wi, &pi, None)?, 8.0 / 9.0);

    // test with gap index - the index is wrong
    assert_ulps_eq!(avg_transition(&Wi, &pi, Some(1),)?, 8.0 / 9.0);

    Ok(())
  }

  #[rstest]
  fn jc69_creates() -> Result<(), Report> {
    let gtr = jc69()?;

    assert_eq!(gtr.pi, array![0.2, 0.2, 0.2, 0.2, 0.2]);

    assert_ulps_eq!(gtr.mu, 0.8);

    let diag = Array2::from_diag(&gtr.eigvals);
    assert_abs_diff_eq!(
      gtr.v.dot(&diag).dot(&gtr.v_inv),
      array![
        [-1.00, 0.25, 0.25, 0.25, 0.25],
        [0.25, -1.00, 0.25, 0.25, 0.25],
        [0.25, 0.25, -1.00, 0.25, 0.25],
        [0.25, 0.25, 0.25, -1.00, 0.25],
        [0.25, 0.25, 0.25, 0.25, -1.00],
      ],
      epsilon = 1e-14
    );

    #[rustfmt::skip]
    assert_ulps_eq!(
      gtr.W,
      array![
        [ 0.0000000000000000,  1.2499999999999993,  1.2499999999999993,  1.2499999999999993,  1.2499999999999993],
        [ 1.2499999999999993,  0.0000000000000000,  1.2499999999999993,  1.2499999999999993,  1.2499999999999993],
        [ 1.2499999999999993,  1.2499999999999993,  0.0000000000000000,  1.2499999999999993,  1.2499999999999993],
        [ 1.2499999999999993,  1.2499999999999993,  1.2499999999999993,  0.0000000000000000,  1.2499999999999993],
        [ 1.2499999999999993,  1.2499999999999993,  1.2499999999999993,  1.2499999999999993,  0.0000000000000000]
      ]
    );

    #[rustfmt::skip]
    assert_abs_diff_eq!(
      gtr.eigvals,
      array![-1.2499999999999995559107901499374, -1.2499999999999995559107901499374, -1.2499999999999993338661852249061, -1.2499999999999991118215802998748, -0.0000000000000006661338147750936],
      epsilon = 1e-14
    );

    #[rustfmt::skip]
    assert_abs_diff_eq!(
      gtr.v,

      // actual:
      // [
      //   [  0.0,                   -0.2000000000000003,       -0.4285714285714284,  -0.0,                   -0.2                 ],
      //   [ -2.912503921549619e-17,  0.49999999999999994,      -0.07142857142857173, -7.648537249746667e-17, -0.2                 ],
      //   [ -0.19039613892487,      -0.1, 0.16666666666666669, -0.49999999999999983,                         -0.2                 ],
      //   [ -0.30960386107512994,   -0.09999999999999988,       0.16666666666666666,  0.4263789812020191,    -0.19999999999999998 ],
      //   [  0.5000000000000001,    -0.09999999999999988,       0.16666666666666666,  0.07362101879798083,   -0.19999999999999998 ],
      // ]

      // numpy:
      array![
        [ 0.000000000000000000000000000000000,  0.000000000000000000000000000000000,  0.499999999999999833466546306226519, -0.043155100467879077552613154011851, -0.199999999999999955591079014993738],
        [-0.000000000000000025495124165822226,  0.000000000000000000000000000000000, -0.168258119423281476967346748097043, -0.456844899532120984897431981153204, -0.199999999999999983346654630622652],
        [ 0.500000000000000000000000000000000, -0.000000000000000022262105026958554, -0.110580626858906211351651904806204,  0.166666666666666657414808128123695, -0.199999999999999983346654630622652],
        [-0.249999999999999916733273153113259, -0.500000000000000000000000000000000, -0.110580626858906239107227520435117,  0.166666666666666601903656896865868, -0.199999999999999983346654630622652],
        [-0.250000000000000055511151231257827,  0.500000000000000000000000000000000, -0.110580626858906239107227520435117,  0.166666666666666601903656896865868, -0.199999999999999983346654630622652],
      ],
      epsilon = 1e-14
    );

    #[rustfmt::skip]
    assert_abs_diff_eq!(
      gtr.v_inv,

      // actual:
      // [
      //   [  0.0,                -7.622255893850872e-17,  -0.49828193581088565, -0.8102580866511949,  1.308540022462081   ],
      //   [ -0.6250000000000011,  1.5625000000000004,     -0.31250000000000017, -0.3124999999999998, -0.3124999999999998  ],
      //   [ -1.5750000000000004, -0.2625000000000013,      0.6125000000000005,   0.6125000000000004,  0.6125000000000004  ],
      //   [ -0.0,                -1.7493603148484694e-16, -1.1435914199845278,   0.9752066891287469,  0.16838473085578107 ],
      //   [ -0.9999999999999999, -0.9999999999999999,     -0.9999999999999999,  -0.9999999999999998, -0.9999999999999998  ],
      // ]

      // numpy:
      array![
        [ 0.000000000000000000000000000000000, -0.000000000000000067986997775525899,  1.333333333333332815229255174926948, -0.666666666666666185570022662432166, -0.666666666666666629659232512494782],
        [ 0.000000000000000000000000000000000,  0.000000000000000000000000000000000, -0.000000000000000044524210053917096, -0.999999999999999777955395074968692,  0.999999999999999777955395074968692],
        [ 1.587326682879073302956385305151343, -0.534161205143256712091215376858599, -0.351055159245272418999661567795556, -0.351055159245272474510812799053383, -0.351055159245272474510812799053383],
        [-0.146834522264184169637601939939486, -1.554407285450700948459257233480457,  0.567080602571628178409923748404253,  0.567080602571627956365318823372945,  0.567080602571627956365318823372945],
        [-1.000000000000000666133814775093924, -1.000000000000000888178419700125232, -1.000000000000000888178419700125232, -1.000000000000000888178419700125232, -1.000000000000000888178419700125232],
      ],
      epsilon = 1e-14
    );

    assert_eq!(gtr.alphabet, Alphabet::new("nuc")?);

    // assert_eq!(gtr.profile_map, {
    //   '-': array([0.0, 0.0, 0.0, 0.0, 1.0]),
    //   'A': array([1.0, 0.0, 0.0, 0.0, 0.0]),
    //   'B': array([0.0, 1.0, 1.0, 1.0, 0.0]),
    //   'C': array([0.0, 1.0, 0.0, 0.0, 0.0]),
    //   'D': array([1.0, 0.0, 1.0, 1.0, 0.0]),
    //   'G': array([0.0, 0.0, 1.0, 0.0, 0.0]),
    //   'H': array([1.0, 1.0, 0.0, 1.0, 0.0]),
    //   'K': array([0.0, 0.0, 1.0, 1.0, 0.0]),
    //   'M': array([1.0, 1.0, 0.0, 0.0, 0.0]),
    //   'N': array([1.0, 1.0, 1.0, 1.0, 1.0]),
    //   'R': array([1.0, 0.0, 1.0, 0.0, 0.0]),
    //   'S': array([0.0, 1.0, 1.0, 0.0, 0.0]),
    //   'T': array([0.0, 0.0, 0.0, 1.0, 0.0]),
    //   'V': array([1.0, 1.0, 1.0, 0.0, 0.0]),
    //   'W': array([1.0, 0.0, 0.0, 1.0, 0.0]),
    //   'X': array([1.0, 1.0, 1.0, 1.0, 1.0]),
    //   'Y': array([0.0, 1.0, 0.0, 1.0, 0.0])
    // });

    Ok(())
  }

  #[rstest]
  fn jc69_calculates_exp_qt() -> Result<(), Report> {
    let gtr = jc69()?;

    let t = (1.0 / 5.0).ln() / gtr.mu;
    let Qs = gtr.expQt(t);

    assert_abs_diff_eq!(
      Qs,
      array![
        [6.18139512, 0., 0., 0., 0.],
        [0., 6.18139512, 0., 0., 0.],
        [0., 0., 6.18139512, 0., 0.],
        [0., 0., 0., 6.18139512, 0.],
        [0., 0., 0., 0., 6.18139512]
      ],
      epsilon = 1e-8
    );

    Ok(())
  }

  #[rstest]
  fn gtr_test_theoretical_limits() -> Result<(), Report> {
    // symmetric rate matrix with some variation in entries (test doesn't depend on precise values)
    let W: Array2<f64> = array![
      [0.00, 1.25, 2.25, 1.25, 1.25],
      [1.25, 0.00, 1.25, 3.25, 1.25],
      [2.25, 1.25, 0.00, 1.25, 1.25],
      [1.25, 3.25, 1.25, 0.00, 1.25],
      [1.25, 1.25, 1.25, 1.25, 0.00],
    ];

    // pi vector of equilibrium probability, some variation to be general, doesn't depend on details
    let pi: Array1<f64> = array![0.18, 0.35, 0.25, 0.18, 0.04];
    let alphabet = Alphabet::new("nuc")?;
    let profile_map = ProfileMap::from_alphabet(&alphabet)?;
    let mu = 1.0;

    let params = GTRParams {
      alphabet,
      profile_map,
      mu,
      W,
      pi,
    };

    let gtr = GTR::new(&params)?;

    // initial profile to be back_propagated or evolved
    let profile: Array2<f64> = array![[0.00, 0.8, 0.0, 0.2, 0.0],];
    // propagate forward and backward in time by a large amount (as long as exp(-mu*large_t) is tiny, value of large_t doesn't matter)
    let large_t = 100.0;
    let distant_past = gtr.propagate_profile(&profile, large_t, false);
    let distant_future = gtr.evolve(&profile, large_t, false);

    // the "distant past profile" is the product of a vector [1,1,1,1,1] times the dot product of pi and initial profile
    let mut weight = 0.0;
    for i in 0..5 {
      weight += params.pi[i] * profile[[0, i]];
    }

    assert_ulps_eq!(
      distant_past,
      array![[1.0, 1.0, 1.0, 1.0, 1.0]] * weight,
      epsilon = 1e-14
    );

    // propagating the profile far into the future gives the equilibrium probabilities pi
    assert_ulps_eq!(distant_future.slice(s![0, ..]), params.pi, epsilon = 1e-14);

    Ok(())
  }
}
