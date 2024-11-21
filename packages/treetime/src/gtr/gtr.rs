use crate::alphabet::alphabet::Alphabet;
use crate::utils::ndarray::{clamp_min, outer};
use eyre::Report;
use ndarray::prelude::*;
use ndarray_linalg::Eigh;
use ndarray_linalg::UPLO::Lower;
use num_traits::abs;

pub fn avg_transition(W: &Array2<f64>, pi: &Array1<f64>) -> Result<f64, Report> {
  Ok(pi.dot(W).dot(pi))
}

/// Performs eigendecomposition of the rate matrix and stores the left- and right-
/// matrices to convert the sequence profiles to the GTR matrix eigenspace
/// and hence to speed up the computations.
/// NOTE: this assumes the diagonal of W is all zeros
#[allow(clippy::type_complexity)]
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
  pub mu: f64,
  pub W: Option<Array2<f64>>,
  pub pi: Array1<f64>,
}

/// Defines General-Time-Reversible model of character evolution.
#[derive(Clone, Debug)]
pub struct GTR {
  pub debug: bool,
  pub is_site_specific: bool,
  pub average_rate: f64,
  pub mu: f64,
  pub W: Array2<f64>,
  pub pi: Array1<f64>,
  pub eigvals: Array1<f64>,
  pub v: Array2<f64>,
  pub v_inv: Array2<f64>,
}

impl GTR {
  pub fn new(GTRParams { alphabet, mu, W, pi }: GTRParams) -> Result<Self, Report> {
    let n = alphabet.n_canonical();

    assert_eq!(
      pi.shape().to_vec(),
      [n],
      "Length of equilibrium frequency vector (`pi`) does not match the alphabet length"
    );

    if let Some(W) = &W {
      assert_eq!(
        W.shape().to_vec(),
        [n, n],
        "Dimensions of substitution matrix (`W`) don't match the alphabet size"
      );
    }

    // self.state_index= {s:si for si,s in enumerate(self.alphabet)}
    // self.state_index.update({s:si for si,s in enumerate(self.alphabet)})

    // self.ambiguous = None
    // let gap_index = Self::assign_gap_and_ambiguous(alphabet);

    let W = {
      let W = W.unwrap_or_else(|| {
        let mut W = Array2::<f64>::ones([n, n]);
        W.diag_mut().fill(0.0);
        let s = -W.sum_axis(Axis(0));
        W.diag_mut().assign(&s);
        W
      });
      let mut W = 0.5 * (&W.view() + &W.t());
      W.diag_mut().fill(0.0);
      W
    };

    let pi = {
      let pi_sum = pi.sum();
      pi / pi_sum
    };

    // FIXME: alphabet.gap_index() does not exist anymore
    let average_rate = avg_transition(&W, &pi)?;
    let mu = mu * average_rate;
    let W = W / average_rate;

    let (eigvals, v, v_inv) = eig_single_site(&W, &pi)?;

    Ok(Self {
      debug: false,
      is_site_specific: false,
      average_rate,
      mu,
      W,
      pi,
      eigvals,
      v,
      v_inv,
    })
  }

  pub const fn average_rate(&self) -> f64 {
    self.average_rate
  }

  // fn assign_gap_and_ambiguous(alphabet: &Alphabet) -> Option<usize> {
  //   // let n_states = self.alphabet.len();
  //
  //   // // determine if a character exists that corresponds to no info, i.e. all one profile
  //   // if any([x.sum()==n_states for x in self.profile_map.values()]):
  //   //     amb_states = [c for c,x in self.profile_map.items() if x.sum()==n_states]
  //   //     self.ambiguous = 'N' if 'N' in amb_states else amb_states[0]
  //   // else:
  //   //     self.ambiguous=None
  //   //
  //   // // check for a gap symbol
  //   // try:
  //   //     self.gap_index = self.state_index['-']
  //   // except:
  //   //     self.gap_index=None
  //
  //   None
  // }

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

  /// Matrix exponential of exp(Qt)
  pub fn expQt(&self, t: f64) -> Array2<f64> {
    let eLambdaT: Array2<f64> = Array2::from_diag(&self.exp_lt(t)); // vector length = a

    let eLambdaT_dot_v_inv: Array2<f64> = eLambdaT.dot(&self.v_inv);

    let Qt: Array2<f64> = self.v.dot(&eLambdaT_dot_v_inv); // This is P(nuc1 | given nuc_2)

    clamp_min(&Qt, 0.0)
  }

  fn exp_lt(&self, t: f64) -> Array1<f64> {
    (self.mu * t * &self.eigvals).mapv(f64::exp)
  }

  pub fn is_multi_site(&self) -> bool {
    self.pi.shape().len() > 1
  }

  /// Product of the transition matrix and the equilibrium frequencies to obtain the rate matrix
  /// of the GTR model
  pub fn Q(&self) -> Array2<f64> {
    let mut Q = (&self.W * &self.pi).t().to_owned();
    let diag = -Q.sum_axis(Axis(0));
    Q.diag_mut().assign(&diag);
    Q
  }

  // pub fn print<W: Write>(&self, w: &mut W) -> Result<(), Report> {
  //   if self.is_multi_site() {
  //     writeln!(w, "Average substitution rate (mu): {:.6}", self.average_rate)?;
  //   } else {
  //     writeln!(w, "Substitution rate (mu): {:.6}", self.mu)?;
  //     writeln!(w, "\nEquilibrium frequencies (pi_i):")?;
  //     for (a, p) in zip(self.alphabet.canonical(), &self.pi) {
  //       writeln!(w, "{a}:\t{p:.4}")?;
  //     }
  //   }
  //
  //   writeln!(w, "\nSymmetrized rates from j->i (W_ij):")?;
  //   writeln!(w, "\t{}", self.alphabet.canonical().join("\t"))?;
  //   for (a, Wi) in zip(self.alphabet.canonical(), self.W.rows()) {
  //     writeln!(
  //       w,
  //       "{a}\t{}",
  //       Wi.iter().map(|Wij| format!("{:.4}", Wij.max(0.0))).join("\t")
  //     )?;
  //   }
  //
  //   if !self.is_multi_site() {
  //     writeln!(w, "\nActual rates from j->i (Q_ij):")?;
  //     writeln!(w, "\t{}", self.alphabet.canonical().join("\t"))?;
  //     for (a, Qi) in zip(self.alphabet.canonical(), self.Q().rows()) {
  //       writeln!(
  //         w,
  //         "{a}\t{}",
  //         Qi.iter().map(|Qij| format!("{:.4}", Qij.max(0.0))).join("\t")
  //       )?;
  //     }
  //   }
  //   writeln!(w)?;
  //   Ok(())
  // }
}

// impl Display for GTR {
//   fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//     let mut buf = vec![];
//     self.print(&mut buf).unwrap();
//     write!(f, "{}", String::from_utf8(buf).unwrap())
//   }
// }

#[cfg(test)]
mod tests {
  #![allow(clippy::excessive_precision)]

  use super::*;
  use crate::alphabet::alphabet::AlphabetName;
  use crate::gtr::get_gtr::{jc69, JC69Params};
  use crate::pretty_assert_ulps_eq;
  use approx::{assert_abs_diff_eq, assert_ulps_eq};
  use eyre::Report;
  use lazy_static::lazy_static;
  use ndarray::{array, Array1, Array2};
  use rstest::rstest;

  lazy_static! {
    static ref ALPHABET: Alphabet = Alphabet::new(AlphabetName::Nuc, false).unwrap();
  }

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
    pretty_assert_ulps_eq!(
      eigvals,
      array![-1.25, -1.25, -1.25, -1.25, 0.0],
      epsilon = 1e-12,
    );

    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      v,
      array![
        [0.3529411764705881,     0.3157894736842106,   0.0,                  0.0,                 0.20000000000000007],
        [-0.44117647058823545,   0.1842105263157894,  -0.07026264899630773,  0.15000000000000005, 0.2                ],
        [-0.029411764705882307, -0.12280701754385968, -0.2657911700123075,  -0.41515307716504657, 0.2                ],
        [-0.029411764705882335, -0.12280701754385966,  0.49999999999999994, -0.08484692283495342, 0.2                ],
        [0.14705882352941183,   -0.2543859649122807,  -0.16394618099138472,  0.3499999999999999,  0.20000000000000007]
      ],
      epsilon = 1e-12,
    );

    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      v_inv,
      array![
        [1.03030303030303,  -1.2878787878787885,  -0.08585858585858573, -0.08585858585858581,  0.4292929292929295 ],
        [1.381818181818182,  0.8060606060606058,  -0.5373737373737375,  -0.5373737373737374,  -1.113131313131313  ],
        [0.0,               -0.19934920800230116, -0.7541027842366552,   1.4186001442443257,  -0.46514815200536946],
        [0.0,                0.4621768660251716,  -1.2791609874989829,  -0.26142856591825553,  1.0784126873920665 ],
        [0.9999999999999999, 0.9999999999999997,   0.9999999999999997,   0.9999999999999997,   0.9999999999999999 ],
      ],
      epsilon = 1e-12,
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

    assert_ulps_eq!(avg_transition(&Wi, &pi)?, 8.0 / 9.0);

    Ok(())
  }

  #[rstest]
  fn avg_transition_from_alphabet() -> Result<(), Report> {
    let num_chars = ALPHABET.n_canonical();

    let W = array![
      [0., 1., 1., 1., 1.],
      [1., 0., 1., 1., 1.],
      [1., 1., 0., 1., 1.],
      [1., 1., 1., 0., 1.],
      [1., 1., 1., 1., 0.]
    ];

    let pi = array![0.2, 0.2, 0.2, 0.2, 0.2];

    assert_ulps_eq!(avg_transition(&W, &pi)?, 0.8000000000000005);

    Ok(())
  }

  #[rstest]
  fn creates_gtr_general() {
    // equilibrium distribution
    let pi = array![0.1, 0.15, 0.35, 0.4];

    // symmetric rate matrix
    let W = array![
      [0.0, 0.2, 0.5, 0.2],
      [0.0, 0.0, 0.3, 0.5],
      [0.0, 0.0, 0.0, 0.1],
      [0.0, 0.0, 0.0, 0.0],
    ];

    let W = &W + &W.t();

    let mu = 1.0;

    let gtr = GTR::new(GTRParams {
      alphabet: ALPHABET.clone(),
      W: Some(W),
      pi: pi.clone(),
      mu,
    })
    .unwrap();

    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      gtr.Q(),
      array![
        [-1.6147308781869687, 0.11331444759206799, 0.28328611898017, 0.11331444759206799],
        [0.16997167138810199, -1.8413597733711047, 0.2549575070821529, 0.4249291784702549],
        [0.9915014164305948, 0.5949008498583568, -0.7648725212464589, 0.19830028328611898],
        [0.45325779036827196, 1.13314447592068, 0.22662889518413598, -0.7365439093484418]],
      epsilon = 1e-12
    );

    pretty_assert_ulps_eq!(
      gtr.eigvals,
      array![
        -2.233008645590864,
        -1.855684299231632,
        -0.8688141373304785,
        6.071532165918825e-17
      ],
      epsilon = 1e-12
    );

    // NOTE: If this is failing, then note that it is allowed that columns can have numbers with different signs, depending on conventions of the underlying linear algebra package.
    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      gtr.v,
      array![
        [0.05873974053212546, 0.4712453851141975, 0.0825744262968125, 0.09999999999999999],
        [0.4412602594678747, -0.14612175247186746, -0.06583552802849169, 0.14999999999999986],
        [-0.1745079879591264, -0.3538782475281327, 0.4174255737031872, 0.35000000000000014],
        [-0.3254920120408734, 0.02875461488580234, -0.4341644719715086, 0.39999999999999997]],
      epsilon = 1e-12
    );

    // NOTE: If this is failing, then similarly to the above, certain differences are tolerable, depending on the underlying math implementation used
    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      gtr.v_inv,
      array![
        [0.34871852609386533, 1.7464113836330555, -0.2959990133787373, -0.48308476367946657],
        [1.7306539716673766, -0.3577558518352421, -0.3713205257613536, 0.026400411562555803],
        [0.7744972596093181, -0.4116644690913425, 1.118627080433171, -1.018048834372102],
        [1.0000000000000009, 0.9999999999999999, 1.0000000000000013, 1.0000000000000009]],
      epsilon = 1e-12
    );

    // Orthogonality of Eigenvectors
    pretty_assert_ulps_eq!(gtr.v_inv.dot(&gtr.v), Array2::<f64>::eye(4), epsilon = 1e-14);

    // Equilibrium distribution
    pretty_assert_ulps_eq!(gtr.v.slice(s![.., -1]), pi, epsilon = 1e-14);

    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      gtr.expQt(0.1),
      array![
        [0.9719558927286881, 0.0019852130818534012, 0.004904636725461245, 0.0019750147773543642],
        [0.0029778196227800645, 0.9681216438807416, 0.004420048636871237, 0.0073423860817647905],
        [0.01716622853911446, 0.010313446819366263, 0.986663935318014, 0.0035099569046974977],
        [0.007900059109417457, 0.019579696218039466, 0.0040113793196542066, 0.9871726422361846]],
      epsilon = 1e-12
    );

    // Propagation backwards in time
    let profile = array![
      [1.0, 0.0, 0.0, 0.0],
      [0.0, 1.0, 0.0, 0.0],
      [0.0, 0.0, 1.0, 0.0],
      [0.0, 0.0, 0.0, 1.0],
      [0.0, 1.0, 0.0, 1.0],
      [0.0, 1.0, 1.0, 0.0],
      [1.0, 1.0, 1.0, 1.0],
    ];

    // Propagate for short time
    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      gtr.propagate_profile(&profile, 0.1, false),
      array![
        [0.97195589, 0.00198521, 0.00490464, 0.00197501],
        [0.00297782, 0.96812164, 0.00442005, 0.00734239],
        [0.01716623, 0.01031345, 0.98666394, 0.00350996],
        [0.00790006, 0.0195797 , 0.00401138, 0.98717264],
        [0.01087788, 0.98770134, 0.00843143, 0.99451503],
        [0.02014405, 0.97843509, 0.99108398, 0.01085234],
        [1.        , 1.        , 1.        , 1.        ]],
      epsilon = 1e-8
    );

    // Propagate for long time
    #[rustfmt::skip]
    assert_ulps_eq!(
      gtr.propagate_profile(&profile, 1000.0, false),
      array![
        [0.1 , 0.1 , 0.1 , 0.1 ],
        [0.15, 0.15, 0.15, 0.15],
        [0.35, 0.35, 0.35, 0.35],
        [0.4 , 0.4 , 0.4 , 0.4 ],
        [0.55, 0.55, 0.55, 0.55],
        [0.5 , 0.5 , 0.5 , 0.5 ],
        [1.  , 1.  , 1.  , 1.  ]],
      epsilon = 1e-12
    );

    // Propagation forward in time
    let profile = array![
      [1.00, 0.00, 0.00, 0.00],
      [0.00, 1.00, 0.00, 0.00],
      [0.00, 0.00, 1.00, 0.00],
      [0.00, 0.00, 0.00, 1.00],
      [0.00, 0.30, 0.20, 0.50],
      [0.00, 0.80, 0.20, 0.00],
      [0.10, 0.15, 0.35, 0.40], // pi
    ];

    // Evolve for short time
    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      gtr.evolve(&profile, 0.1, false),
      array![
        [0.97195589, 0.00297782, 0.01716623, 0.00790006],
        [0.00198521, 0.96812164, 0.01031345, 0.0195797 ],
        [0.00490464, 0.00442005, 0.98666394, 0.00401138],
        [0.00197501, 0.00734239, 0.00350996, 0.98717264],
        [0.002564  , 0.2949917 , 0.2021818 , 0.50026251],
        [0.0025691 , 0.77538132, 0.20558354, 0.01646603],
        [0.1       , 0.15      , 0.35      , 0.4       ]],
      epsilon = 1e-8
    );

    // Evolve for long time
    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      gtr.evolve(&profile, 1000.0, false),
      array![
        [0.10, 0.15, 0.35, 0.40], // pi
        [0.10, 0.15, 0.35, 0.40], // pi
        [0.10, 0.15, 0.35, 0.40], // pi
        [0.10, 0.15, 0.35, 0.40], // pi
        [0.10, 0.15, 0.35, 0.40], // pi
        [0.10, 0.15, 0.35, 0.40], // pi
        [0.10, 0.15, 0.35, 0.40], // pi
      ],
      epsilon = 1e-12
    );
  }

  #[rstest]
  fn jc69_creates() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;

    assert_eq!(gtr.pi, array![0.25, 0.25, 0.25, 0.25]);

    pretty_assert_ulps_eq!(gtr.mu, 0.75);

    let diag = Array2::from_diag(&gtr.eigvals);
    let od: f64 = 1.0 / 3.0;
    pretty_assert_ulps_eq!(
      gtr.v.dot(&diag).dot(&gtr.v_inv),
      array![
        [-1.00, od, od, od],
        [od, -1.00, od, od],
        [od, od, -1.00, od],
        [od, od, od, -1.00],
      ],
      epsilon = 1e-14
    );

    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      gtr.W,
      array![
        [ 0.0000000000000000,  1.0+od,  1.0+od,  1.0+od],
        [ 1.0+od,  0.0000000000000000,  1.0+od,  1.0+od],
        [ 1.0+od,  1.0+od,  0.0000000000000000,  1.0+od],
        [ 1.0+od,  1.0+od,  1.0+od,  0.0000000000000000],
      ]
    );

    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      gtr.eigvals,
      array![-1.0 - od, -1.0 - od, -1.0 - od, 0.0]
    );

    // NOTE: These are failing likely due to the different conventions in linear algebra packages used in Python version
    // (the source of expected test values) and in Rust version. And these differences are hard to test. But this
    // also tests internals of the model implementation, which is a kind of antipattern. The most important thing is
    // to ensure that the model usage produces correct results. While the differences in the internals might be
    // tolerable.

    // #[rustfmt::skip]
    // pretty_assert_ulps_eq!(
    //   gtr.v,
    //
    //   // Rust result
    //   // [
    //   //   [0.38888888888888884,   0.0,                    -0.25000000000000006,  0.0,                0.2],
    //   //   [0.11111111111111116,   1.1518476660775412e-16,  0.4999999999999998,   0.0,                0.2],
    //   //   [-0.16666666666666663, -0.5000000000000001,     -0.08333333333333331, -0.1744576301870094, 0.19999999999999998],
    //   //   [-0.16666666666666669,  0.09150635094610969,    -0.08333333333333341,  0.5,                0.2],
    //   //   [-0.16666666666666669,  0.40849364905389013,    -0.08333333333333341, -0.3255423698129907, 0.2]
    //   // ]
    //
    //   // Python result (openblas backend)
    //   array![
    //     [-0.046874760469663719,  0.500000000000000111, -0.001426887227439782,  0.007485879491300664, -0.199999999999999956],
    //     [-0.453125239530336232, -0.172109713611864223,  0.022158582410371695, -0.004910409593215660, -0.200000000000000094],
    //     [ 0.176921578205358310, -0.102772856485764411,  0.229538218058494475, -0.495089590406784297, -0.200000000000000011],
    //     [ 0.174649639361912162, -0.112629155793391236,  0.248303199531133917,  0.481091398682775395, -0.199999999999999983],
    //     [ 0.148428782432729611, -0.112488274108980033, -0.498573112772560278,  0.011422721825924000, -0.199999999999999983],
    //   ]
    // );

    // #[rustfmt::skip]
    // pretty_assert_ulps_eq!(
    //   gtr.v_inv,
    //
    //   // Rust result
    //   // [
    //   //   [1.5750000000000002,  0.45000000000000023,   -0.6749999999999999,  -0.6750000000000002,  -0.6750000000000002],
    //   //   [0.0,                 2.708697166989168e-16, -1.1758052938602825,   0.21518730372854528,  0.9606179901317368],
    //   //   [-0.7500000000000001, 1.4999999999999996,    -0.25,                -0.2500000000000002,  -0.2500000000000002],
    //   //   [0.0,                 0.0,                   -0.4514793629381211,   1.2939513234650695,  -0.8424719605269488],
    //   //   [0.9999999999999997,  0.9999999999999997,     0.9999999999999993,   0.9999999999999997,   0.9999999999999997],
    //   // ]
    //
    //   // Python result (openblas backend)
    //   array![
    //    [-0.160885619055681939, -1.555236420221760563,  0.607238039163926158,  0.599440190521669192,  0.509443809591847541],
    //    [ 1.584670771738339923, -0.545474465385955010, -0.325722283602099816, -0.356960262462701905, -0.356513760287582582],
    //    [-0.003926379079478472,  0.060973980798110977,  0.631622485641511799,  0.683258262642032932, -1.371928350002177277],
    //    [ 0.015701130983909538, -0.010299255324283767, -1.038417265036203574,  1.009056995203920337,  0.023958394172657615],
    //    [-0.999999999999999556, -1.000000000000000222, -0.999999999999999889, -0.999999999999999778, -0.999999999999999778],
    // ]);

    Ok(())
  }

  #[rstest]
  fn jc69_calculates_exp_qt() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;

    let t = 0.2 / gtr.mu;
    let Qs = gtr.expQt(t);

    assert_abs_diff_eq!(
      Qs,
      array![
        [0.82444625, 0.05851792, 0.05851792, 0.05851792],
        [0.05851792, 0.82444625, 0.05851792, 0.05851792],
        [0.05851792, 0.05851792, 0.82444625, 0.05851792],
        [0.05851792, 0.05851792, 0.05851792, 0.82444625]
      ],
      epsilon = 1e-8
    );

    Ok(())
  }

  #[rstest]
  fn gtr_test_theoretical_limits() -> Result<(), Report> {
    // symmetric rate matrix with some variation in entries (test doesn't depend on precise values)
    let W: Array2<f64> = array![
      [0.00, 1.25, 2.25, 1.25],
      [1.25, 0.00, 1.25, 3.25],
      [2.25, 1.25, 0.00, 1.25],
      [1.25, 3.25, 1.25, 0.00],
    ];

    // pi vector of equilibrium probability, some variation to be general, doesn't depend on details
    let pi: Array1<f64> = array![0.18, 0.35, 0.25, 0.22];
    let mu = 1.0;

    // initial profile to be back_propagated or evolved
    let profile: Array2<f64> = array![[0.00, 0.8, 0.0, 0.2],];

    let params = GTRParams {
      alphabet: ALPHABET.clone(),
      mu,
      W: Some(W),
      pi: pi.clone(),
    };

    let gtr = GTR::new(params)?;

    // propagate forward and backward in time by a large amount (as long as exp(-mu*large_t) is tiny, value of large_t doesn't matter)
    let large_t = 100.0;
    let distant_past = gtr.propagate_profile(&profile, large_t, false);
    let distant_future = gtr.evolve(&profile, large_t, false);

    // the "distant past profile" is the product of a vector [1,1,1,1,1] times the dot product of pi and initial profile
    let mut weight = 0.0;
    for i in 0..4 {
      weight += pi[i] * profile[[0, i]];
    }
    let distant_past_expected = array![[1.0, 1.0, 1.0, 1.0]] * weight;
    assert_ulps_eq!(distant_past, distant_past_expected, epsilon = 1e-10);

    // propagating the profile far into the future gives the equilibrium probabilities pi
    assert_ulps_eq!(distant_future.slice(s![0, ..]), pi, epsilon = 1e-10);

    Ok(())
  }
}
