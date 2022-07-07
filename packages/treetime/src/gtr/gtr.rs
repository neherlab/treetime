#![allow(non_snake_case)]

use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::profile_map::ProfileMap;
use crate::make_report;
use crate::utils::ndarray::{clamp_min, outer, to_col};
use eyre::Report;
use ndarray::prelude::*;
use ndarray_einsum_beta::einsum;
use ndarray_linalg::UPLO::Lower;
use ndarray_linalg::{Eig, Eigh};
use num_traits::abs;
use num_traits::real::Real;

pub fn avg_transition(W: &Array2<f32>, pi: &Array1<f32>, gap_index: Option<usize>) -> Result<f32, Report> {
  let result = einsum("i,ij,j", &[pi, W, pi])
    .map_err(|err| make_report!("einsum: {err}"))?
    .into_dimensionality::<Ix0>()?
    .into_scalar();

  Ok(if let Some(gap_index) = gap_index {
    // np.sum(pi*W[:,gap_index]) *pi[gap_index])/(1-pi[gap_index])
    let W_slice = W.slice(s!(.., gap_index));
    let pi_mul_W_slice = pi * &W_slice;
    let pi_mul_W_slice_sum = pi_mul_W_slice.sum();
    let result = result - pi_mul_W_slice_sum * &pi[gap_index] / (1.0 - pi[gap_index]);
    result
  } else {
    result
  })
}

/// Performs eigendecomposition of the rate matrix and stores the left- and right-
/// matrices to convert the sequence profiles to the GTR matrix eigenspace
/// and hence to speed-up the computations.
/// NOTE: this assumes the diagonal of W is all zeros
fn eig_single_site(W: &Array2<f32>, pi: &Array1<f32>) -> Result<(Array1<f32>, Array2<f32>, Array2<f32>), Report> {
  assert!(abs(W.diag().sum()) < 1e-10);

  let tmp_pi: Array1<f32> = pi.mapv(f32::sqrt);
  let mut sym_Q: Array2<f32> = W * outer(&tmp_pi, &tmp_pi)?;

  let diag = -(W * pi).sum_axis(Axis(1));
  sym_Q.diag_mut().assign(&diag);

  let (eigvals, eigvecs) = sym_Q.eigh(Lower)?;

  let tmp_v: Array2<f32> = eigvecs.t().to_owned().dot(&to_col(&tmp_pi)?);
  let one_norm: Array1<f32> = tmp_v.mapv(f32::abs).sum_axis(Axis(1));

  let v = tmp_v / &one_norm;
  let v_inv = (eigvecs * one_norm) / tmp_pi;

  Ok((eigvals, v, v_inv))
}

#[derive(Clone, Debug)]
pub struct GTRParams {
  pub alphabet: Alphabet,
  pub profile_map: ProfileMap,
  pub mu: f32,
  pub W: Array2<f32>,
  pub pi: Array1<f32>,
}

/// Defines General-Time-Reversible model of character evolution.
#[derive(Clone, Debug)]
pub struct GTR {
  pub debug: bool,
  pub is_site_specific: bool,
  pub alphabet: Alphabet,
  pub profile_map: ProfileMap,
  pub mu: f32,
  pub W: Array2<f32>,
  pub pi: Array1<f32>,
  pub eigvals: Array1<f32>,
  pub v: Array2<f32>,
  pub v_inv: Array2<f32>,
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

    // self.state_index= {s:si for si,s in enumerate(self.alphabet)}
    // self.state_index.update({s:si for si,s in enumerate(self.alphabet)})

    // self.ambiguous = None
    // self.gap_index = None
    // self.assign_gap_and_ambiguous();

    let n = alphabet.len();

    let mut pi = if pi.len() == n { pi.clone() } else { Array1::ones(n) };
    let pi_sum = pi.sum();
    pi = pi / pi_sum;

    let mut W: Array2<f32> = 0.5 * W + W.t();
    W.diag_mut().fill(0.0);

    let average_rate = avg_transition(&W, &pi, None)?;
    W = W / average_rate;
    let mu = mu * average_rate;

    let (eigvals, v, v_inv) = eig_single_site(&W, &pi)?;

    Ok(Self {
      debug: false,
      is_site_specific: false,
      alphabet: alphabet.to_owned(),
      profile_map: profile_map.to_owned(),
      mu,
      W,
      pi,
      eigvals,
      v,
      v_inv,
    })
  }

  fn assign_gap_and_ambiguous(&self) {
    let n_states = self.alphabet.len();

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
  pub fn evolve(&self, profile: &Array2<f32>, t: f32, return_log: bool) -> Array2<f32> {
    let Qt = self.expQt(t);
    let res = profile.dot(&Qt.t());
    if return_log {
      res.mapv(f32::ln)
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
  pub fn propagate_profile(&self, profile: &Array2<f32>, t: f32, return_log: bool) -> Array2<f32> {
    let Qt = self.expQt(t);
    let res = profile.dot(&Qt);
    if return_log {
      res.mapv(f32::ln)
    } else {
      res
    }
  }

  /// Matrix exponential of exo(Qt)
  fn expQt(&self, t: f32) -> Array2<f32> {
    let eLambdaT: Array2<f32> = Array2::from_diag(&self.exp_lt(t)); // vector length = a

    let eLambdaT_dot_v_inv: Array2<f32> = eLambdaT.dot(&self.v_inv);

    let Qs: Array2<f32> = self.v.dot(&eLambdaT_dot_v_inv); // This is P(nuc1 | given nuc_2)

    clamp_min(&Qs, 0.0)
  }

  fn exp_lt(&self, t: f32) -> Array1<f32> {
    (self.mu * t * &self.eigvals).mapv(f32::exp)
  }
}

#[allow(clippy::excessive_precision)]
#[cfg(test)]
mod test {
  use super::*;
  use crate::nuc_models::jc69::jc69;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use ndarray::{array, Array1, Array2};
  use rstest::rstest;

  #[rstest]
  fn computes_eig_single_site() -> Result<(), Report> {
    let W: Array2<f32> = array![
      [0.00, 1.25, 1.25, 1.25, 1.25],
      [1.25, 0.00, 1.25, 1.25, 1.25],
      [1.25, 1.25, 0.00, 1.25, 1.25],
      [1.25, 1.25, 1.25, 0.00, 1.25],
      [1.25, 1.25, 1.25, 1.25, 0.00],
    ];

    let pi: Array1<f32> = array![0.2, 0.2, 0.2, 0.2, 0.2];

    let (eigvals, v, v_inv) = eig_single_site(&W, &pi)?;

    #[rustfmt::skip]
    assert_ulps_eq!(eigvals, array![-1.25000000e+00, -1.25000000e+00, -1.25000000e+00, -1.25000000e+00, -6.66133815e-16]);

    #[rustfmt::skip]
    assert_ulps_eq!(v, array![
      [0.00000000e+00,0.00000000e+00,5.00000000e-01,-4.31551005e-02,-2.00000000e-01],
      [-2.54951242e-17,0.00000000e+00,-1.68258119e-01,-4.56844900e-01,-2.00000000e-01],
      [5.00000000e-01,-2.22621050e-17,-1.10580627e-01,1.66666667e-01,-2.00000000e-01],
      [-2.50000000e-01,-5.00000000e-01,-1.10580627e-01,1.66666667e-01,-2.00000000e-01],
      [-2.50000000e-01,5.00000000e-01,-1.10580627e-01,1.66666667e-01,-2.00000000e-01]
    ]);

    #[rustfmt::skip]
    assert_ulps_eq!(v_inv, array![
      [0.00000000e+00,-6.79869978e-17,1.33333333e+00,-6.66666667e-01,-6.66666667e-01],
      [0.00000000e+00,0.00000000e+00,-4.45242101e-17,-1.00000000e+00,1.00000000e+00],
      [1.58732668e+00,-5.34161205e-01,-3.51055159e-01,-3.51055159e-01,-3.51055159e-01],
      [-1.46834522e-01,-1.55440729e+00,5.67080603e-01,5.67080603e-01,5.67080603e-01],
      [-1.00000000e+00,-1.00000000e+00,-1.00000000e+00,-1.00000000e+00,-1.00000000e+00]
    ]);

    Ok(())
  }

  #[rstest]
  fn avg_transition_test() -> Result<(), Report> {
    let pi = array![1.0/3.0, 1.0/3.0, 1.0/3.0];

    let Wi: Array2<f32> = array![
                        [0.0, 4.0 / 3.0, 4.0 / 3.0], 
                        [4.0 / 3.0, 0.0, 4.0 / 3.0], 
                        [4.0 / 3.0, 4.0 / 3.0, 0.0],
                        ];
    // test without gap index
    assert_ulps_eq!(avg_transition(
      &Wi,
      &pi,
      None
    )?, 8.0/9.0);

    // // test with gap index - the index is wrong
    // let i: usize = 0;
    // assert_eq!(avg_transition(
    //   &Wi,
    //   &pi,
    //   &i,
    // )?, (8.0/9.0 - 8.0/27.0)*(3.0/2.0));

    Ok(())
  }

  #[rstest]
  fn jc69_creates() -> Result<(), Report> {

    let gtr = jc69()?;

    assert_eq!(gtr.pi, array![0.2, 0.2, 0.2, 0.2, 0.2]);

    // need to add gap index for avg_transition with 'nuc' alphabet as this contains '-'
    // assert_ulps_eq!(gtr.mu, 0.8);
    println!("mu = {}", gtr.mu);

    assert_eq!(
      gtr.W,
      array![
        [0.00, 1.25, 1.25, 1.25, 1.25],
        [1.25, 0.00, 1.25, 1.25, 1.25],
        [1.25, 1.25, 0.00, 1.25, 1.25],
        [1.25, 1.25, 1.25, 0.00, 1.25],
        [1.25, 1.25, 1.25, 1.25, 0.00],
      ]
    );

    #[rustfmt::skip]
    assert_ulps_eq!(
      gtr.eigvals,
      array![-1.25000000e+00,-1.25000000e+00,-1.25000000e+00,-1.25000000e+00,-2.22044605e-16]
    );

    #[rustfmt::skip]
    assert_ulps_eq!(
      gtr.v,
      array![
        [-1.90476190e-01,0.00000000e+00,0.00000000e+00,-5.00000000e-01,-2.00000000e-01],
        [-3.09523810e-01,0.00000000e+00,1.31972474e-17,3.75000000e-01,-2.00000000e-01],
        [1.66666667e-01,-3.66025404e-01,-3.66025404e-01,4.16666667e-02,-2.00000000e-01],
        [1.66666667e-01,5.00000000e-01,-1.33974596e-01,4.16666667e-02,-2.00000000e-01],
        [1.66666667e-01,-1.33974596e-01,5.00000000e-01,4.16666667e-02,-2.00000000e-01]
      ]
    );

    #[rustfmt::skip]
    assert_ulps_eq!(
      gtr.v_inv,
      array![
        [-8.84210526e-01,-1.43684211e+00,7.73684211e-01,7.73684211e-01,7.73684211e-01],
        [0.00000000e+00,0.00000000e+00,-9.10683603e-01,1.24401694e+00,-3.33333333e-01],
        [0.00000000e+00,3.28351985e-17,-9.10683603e-01,-3.33333333e-01,1.24401694e+00],
        [-1.26315789e+00,9.47368421e-01,1.05263158e-01,1.05263158e-01,1.05263158e-01],
        [-1.00000000e+00,-1.00000000e+00,-1.00000000e+00,-1.00000000e+00,-1.00000000e+00]
      ]
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

    let Qs = gtr.expQt(0.1);

    assert_ulps_eq!(
      Qs,
      array![
        [0.90599752, 0.02350062, 0.02350062, 0.02350062, 0.02350062],
        [0.02350062, 0.90599752, 0.02350062, 0.02350062, 0.02350062],
        [0.02350062, 0.02350062, 0.90599752, 0.02350062, 0.02350062],
        [0.02350062, 0.02350062, 0.02350062, 0.90599752, 0.02350062],
        [0.02350062, 0.02350062, 0.02350062, 0.02350062, 0.90599752]
      ]
    );

    Ok(())
  }

  #[rstest]
  fn jc69_evolves() -> Result<(), Report> {
    let gtr = jc69()?;

    let norm_prof: Array2<f32> = array![
      [0.19356424, 0.25224431, 0.21259213, 0.19217803, 0.14942128],
      [0.19440831, 0.13170981, 0.26841564, 0.29005381, 0.11541244],
      [0.27439982, 0.18330691, 0.19687558, 0.32079767, 0.02462001],
      [0.03366488, 0.00781195, 0.32170632, 0.30066296, 0.33615390],
      [0.31185458, 0.25466645, 0.14705881, 0.24872985, 0.03769030],
      [0.24016971, 0.05380214, 0.35454510, 0.19585567, 0.15562739],
      [0.12705805, 0.37184099, 0.21907519, 0.27300161, 0.00902417],
    ];

    let res = gtr.evolve(&norm_prof, 0.1, false);

    assert_ulps_eq!(
      res,
      array![
        [0.19432046, 0.24610544, 0.21111252, 0.19309714, 0.15536444],
        [0.19506535, 0.13973412, 0.26037659, 0.27947221, 0.12535174],
        [0.26565761, 0.18526840, 0.19724271, 0.30660357, 0.04522771],
        [0.05320977, 0.03039464, 0.30740545, 0.28883475, 0.32015539],
        [0.29871132, 0.24824297, 0.15327957, 0.24300395, 0.05676219],
        [0.23544964, 0.07098084, 0.33638557, 0.19634264, 0.16084131],
        [0.13562895, 0.35164914, 0.21683379, 0.26442369, 0.03146442],
      ]
    );

    Ok(())
  }
}
