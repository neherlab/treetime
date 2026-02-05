use crate::alphabet::alphabet::Alphabet;
use eyre::Report;
use ndarray::prelude::*;
use ndarray_linalg::Eigh;
use ndarray_linalg::UPLO::Lower;
use num_traits::abs;
use treetime_utils::ndarray::{clamp_min, outer};

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
    if return_log { res.mapv(f64::ln) } else { res }
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

    if return_log { res.mapv(f64::ln) } else { res }
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

