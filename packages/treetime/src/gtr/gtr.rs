use eyre::Report;
use ndarray::prelude::*;
use ndarray_linalg::Eigh;
use ndarray_linalg::UPLO::Lower;
use num_traits::abs;
use treetime_utils::array::ndarray::{clamp_min, outer};

/// Compute the average substitution rate for normalization.
///
/// Returns pi' * W * pi, which represents the expected number of substitutions
/// per unit time when the system is at equilibrium. Used to normalize W so that
/// branch lengths are measured in expected substitutions per site.
pub fn avg_transition(W: &Array2<f64>, pi: &Array1<f64>) -> Result<f64, Report> {
  Ok(pi.dot(W).dot(pi))
}

/// Eigendecomposition of the GTR rate matrix for efficient matrix exponentiation.
///
/// # Mathematical Background
///
/// The GTR rate matrix Q governs continuous-time Markov evolution: dP/dt = Q * P.
/// The solution P(t) = exp(Q * t) requires matrix exponentiation, which is expensive
/// to compute directly. Eigendecomposition enables efficient computation:
///
///   P(t) = V * exp(Lambda * t) * V^{-1}
///
/// where Lambda is diagonal (eigenvalues) and exp(Lambda * t) is trivial to compute.
///
/// # Symmetrization for Numerical Stability
///
/// GTR satisfies detailed balance: pi[i] * Q[i,j] = pi[j] * Q[j,i], which means
/// Q is similar to a symmetric matrix. We exploit this via the similarity transform:
///
///   S = D^{1/2} * Q * D^{-1/2}    where D = diag(pi)
///
/// For the rate matrix Q with Q[i,j] = W[i,j] * pi[j] (row-stochastic convention):
///
///   S[i,j] = sqrt(pi[i]) * W[i,j] * pi[j] / sqrt(pi[j])
///          = W[i,j] * sqrt(pi[i] * pi[j])
///
/// This S is symmetric (since W is symmetric), guaranteeing:
/// - Real eigenvalues (essential for exp(lambda * t) to be real)
/// - Orthogonal eigenvectors (numerically stable via symmetric eigensolvers)
///
/// # Convention Choice
///
/// Internally we use row-stochastic Q where Q[i,j] = W[i,j] * pi[j] and rows sum to 0.
/// This matches the similarity transform derivation above. The returned v and v_inv
/// matrices incorporate D^{-1/2} and D^{1/2} factors so that expQt() produces
/// column-stochastic P(t) matrices (columns sum to 1) for compatibility with the
/// rest of the codebase.
///
/// Note: The `Q()` display method uses column-stochastic convention (Q[i,j] = W[i,j] * pi[i])
/// for human readability. This is a different mathematical object, not an inconsistency.
///
/// # Arguments
///
/// * `W` - Symmetric exchangeability matrix with zero diagonal. W[i,j] = W[j,i] encodes
///   the relative rate of exchange between states i and j, independent of frequencies.
/// * `pi` - Equilibrium frequency vector summing to 1. pi[i] is the stationary probability
///   of state i.
///
/// # Returns
///
/// * `eigvals` - Eigenvalues of Q (same as eigenvalues of S). All non-positive; exactly one zero.
/// * `v` - Right transformation matrix incorporating D^{-1/2} factor.
/// * `v_inv` - Left transformation matrix incorporating D^{1/2} factor.
///
/// These satisfy: P(t) = v * diag(exp(eigvals * mu * t)) * v_inv
#[allow(clippy::type_complexity)]
pub(super) fn eig_single_site(
  W: &Array2<f64>,
  pi: &Array1<f64>,
) -> Result<(Array1<f64>, Array2<f64>, Array2<f64>), Report> {
  // W must have zero diagonal (off-diagonal rates only)
  assert!(abs(W.diag().sum()) < 1e-10);

  // Build symmetric matrix S = W * outer(sqrt(pi), sqrt(pi))
  // This is the similarity transform of Q: S[i,j] = W[i,j] * sqrt(pi[i] * pi[j])
  let sqrt_pi: Array1<f64> = pi.mapv(f64::sqrt);
  let mut sym_Q: Array2<f64> = W * outer(&sqrt_pi, &sqrt_pi)?;

  // Set diagonal so rows sum to zero: S[i,i] = -sum_{j != i} W[i,j] * pi[j]
  // This corresponds to the row-stochastic Q convention internally
  let diag = -(W * pi).sum_axis(Axis(1));
  sym_Q.diag_mut().assign(&diag);

  // Symmetric eigendecomposition: S = U * Lambda * U'  where U is orthogonal
  let (eigvals, eigvecs) = sym_Q.eigh(Lower)?;

  // Transform eigenvectors to get v and v_inv for the original (non-symmetric) Q.
  // Since S = D^{1/2} Q D^{-1/2}, we have Q = D^{-1/2} S D^{1/2}.
  // If S = U Lambda U', then Q = (D^{-1/2} U) Lambda (U' D^{1/2}).
  // We normalize by one_norm to ensure numerical stability and proper scaling.
  let tmp_v: Array2<f64> = eigvecs.t().to_owned() * sqrt_pi.to_owned();
  let one_norm: Array1<f64> = tmp_v.mapv(f64::abs).sum_axis(Axis(1));

  let v = tmp_v.t().to_owned() / &one_norm;
  let v_inv = (eigvecs * one_norm).t().to_owned() / sqrt_pi;

  Ok((eigvals, v, v_inv))
}

/// Parameters for constructing a GTR model.
#[derive(Clone, Debug)]
pub struct GTRParams {
  /// Number of states (e.g., 4 for nucleotides, 20 for amino acids).
  pub n_states: usize,
  /// Base substitution rate (will be scaled by average transition rate).
  pub mu: f64,
  /// Exchangeability matrix. If None, defaults to all-ones (equal rates).
  /// Will be symmetrized and have diagonal zeroed during construction.
  pub W: Option<Array2<f64>>,
  /// Equilibrium frequencies (will be normalized to sum to 1).
  pub pi: Array1<f64>,
}

/// General Time-Reversible (GTR) model of character evolution.
///
/// GTR is the most general neutral, independent, reversible Markov model of sequence
/// evolution. It describes how characters (nucleotides, amino acids) change over time
/// along phylogenetic branches.
///
/// # Model Components
///
/// The model is parameterized by:
///
/// * **W** - Symmetric exchangeability matrix. W[i,j] = W[j,i] represents the intrinsic
///   rate of exchange between states i and j, independent of their frequencies. For
///   nucleotides, this captures transition/transversion biases.
///
/// * **pi** - Equilibrium frequencies. pi[i] is the stationary probability of state i.
///   The process converges to this distribution as t -> infinity.
///
/// * **mu** - Overall substitution rate, scaling branch lengths. After normalization,
///   mu incorporates the average rate so that branch length t represents expected
///   substitutions per site.
///
/// # Rate Matrix Q
///
/// The instantaneous rate matrix Q determines evolution: dP/dt = Q * P.
///
/// Two equivalent conventions exist for Q (both satisfy detailed balance when W is symmetric):
///
/// * **Row-stochastic**: Q[i,j] = W[i,j] * pi[j], rows sum to 0. Used internally for
///   eigendecomposition because it yields a clean similarity transform.
///
/// * **Column-stochastic**: Q[i,j] = W[i,j] * pi[i], columns sum to 0. Common in textbooks.
///   The `Q()` method returns this form for display/debugging.
///
/// # Transition Probability Matrix P(t)
///
/// P(t) = exp(Q * mu * t) gives transition probabilities over time t.
///
/// This implementation produces column-stochastic P(t): P[i,j] = Prob(state j -> state i),
/// and columns sum to 1. This means column j represents "starting from state j" and
/// entry P[i,j] is the probability of ending in state i.
///
/// # Eigendecomposition
///
/// Computing exp(Q * t) via eigendecomposition: P(t) = v * diag(exp(eigvals * mu * t)) * v_inv.
/// The eigvals, v, and v_inv are precomputed in the constructor for efficiency.
///
/// See `eig_single_site()` for details on the symmetrization trick that ensures
/// numerical stability.
#[derive(Clone, Debug)]
pub struct GTR {
  pub debug: bool,
  /// Average substitution rate before normalization, used for scaling.
  pub average_rate: f64,
  /// Overall substitution rate (incorporates average_rate after normalization).
  pub mu: f64,
  /// Symmetric exchangeability matrix (zero diagonal, normalized by average_rate).
  pub W: Array2<f64>,
  /// Equilibrium frequencies (sum to 1).
  pub pi: Array1<f64>,
  /// Eigenvalues of the rate matrix (all <= 0, exactly one zero).
  pub eigvals: Array1<f64>,
  /// Right transformation matrix for P(t) = v * exp(Lambda * t) * v_inv.
  pub v: Array2<f64>,
  /// Left transformation matrix (inverse of v, adjusted for non-symmetric Q).
  pub v_inv: Array2<f64>,
  /// Per-site rate multipliers for among-site rate variation.
  ///
  /// When `Some`, each element `site_rates[a]` is a relative rate multiplier for
  /// alignment position `a`. Sites with rate > 1 evolve faster than average; sites
  /// with rate < 1 evolve slower. The eigendecomposition (W, pi, eigvals, v, v_inv)
  /// is shared across all sites - only the rate scaling changes per site.
  ///
  /// The effective rate at site `a` is `mu * site_rates[a]`. The matrix exponential
  /// becomes `P_a(t) = V * diag(exp(eigvals * mu * site_rates[a] * t)) * V_inv`.
  pub site_rates: Option<Array1<f64>>,
  /// Whether the one-dimensional branch-length likelihood $L(t)$ is guaranteed
  /// unimodal on $(0, \infty)$ for this model.
  ///
  /// Dinh & Matsen (2017) prove that models with a single distinct nonzero
  /// eigenvalue (JC69, F81, binary symmetric) have at most one stationary point
  /// (Corollary 3.1). With one distinct nonzero eigenvalue, the per-site
  /// characteristic polynomial factors into linear terms with real roots,
  /// satisfying the condition of Theorem 3.1. For these models, a negative
  /// derivative at $t = 0$ guarantees zero is the global maximum on
  /// $[0, \infty)$, enabling the zero-branch shortcut.
  ///
  /// Models with multiple distinct nonzero eigenvalues (K80, HKY85, TN93, GTR)
  /// can have two local maxima. The shortcut is not valid for these models.
  pub unimodal_branch_likelihood: bool,
}

impl GTR {
  /// Construct a new GTR model from parameters.
  ///
  /// # Processing Steps
  ///
  /// 1. **Symmetrize W**: Computes W = (W + W') / 2 to ensure symmetry, then zeros the diagonal.
  ///    Symmetry is required for detailed balance and stable eigendecomposition.
  ///
  /// 2. **Normalize pi**: Scales pi to sum to 1.
  ///
  /// 3. **Normalize W by average rate**: Computes avg_rate = pi' * W * pi (expected substitutions
  ///    per unit time at equilibrium), then divides W by avg_rate and multiplies mu by avg_rate.
  ///    After normalization, branch length t corresponds to t expected substitutions per site.
  ///
  /// 4. **Eigendecomposition**: Calls `eig_single_site()` to precompute eigvals, v, v_inv
  ///    for efficient matrix exponentiation.
  ///
  /// # Panics
  ///
  /// Panics if dimensions of W or pi don't match n_states.
  pub fn new(GTRParams { n_states, mu, W, pi }: GTRParams) -> Result<Self, Report> {
    let n = n_states;

    assert_eq!(
      pi.shape().to_vec(),
      [n],
      "Length of equilibrium frequency vector (`pi`) does not match n_states"
    );

    if let Some(W) = &W {
      assert_eq!(
        W.shape().to_vec(),
        [n, n],
        "Dimensions of substitution matrix (`W`) don't match n_states"
      );
    }

    // Symmetrize W and zero diagonal
    let W = {
      let W = W.unwrap_or_else(|| {
        // Default: equal exchange rates between all states
        let mut W = Array2::<f64>::ones([n, n]);
        W.diag_mut().fill(0.0);
        W
      });
      // Enforce symmetry: W = (W + W') / 2
      let mut W = 0.5 * (&W.view() + &W.t());
      // Diagonal must be zero (off-diagonal rates only; diagonal is set by row/column sum constraint)
      W.diag_mut().fill(0.0);
      W
    };

    // Normalize pi to sum to 1
    let pi = {
      let pi_sum = pi.sum();
      pi / pi_sum
    };

    // Normalize W so that average substitution rate = 1, absorb scaling into mu
    let average_rate = avg_transition(&W, &pi)?;
    let mu = mu * average_rate;
    let W = W / average_rate;

    // Precompute eigendecomposition for efficient exp(Q*t)
    let (eigvals, v, v_inv) = eig_single_site(&W, &pi)?;

    // A 2-state rate matrix has exactly one nonzero eigenvalue, so L(t) is
    // unimodal (Dinh & Matsen 2017, Corollary 3.1). JC69/F81 constructors
    // override this to true for 4-state nucleotide models separately.
    let unimodal_branch_likelihood = n == 2;

    Ok(Self {
      debug: false,
      average_rate,
      mu,
      W,
      pi,
      eigvals,
      v,
      v_inv,
      site_rates: None,
      unimodal_branch_likelihood,
    })
  }

  pub const fn average_rate(&self) -> f64 {
    self.average_rate
  }

  /// Whether per-site rate variation is active.
  pub fn has_site_rates(&self) -> bool {
    self.site_rates.is_some()
  }

  /// Set per-site rate multipliers for among-site rate variation.
  pub fn set_site_rates(&mut self, rates: Array1<f64>) {
    self.site_rates = Some(rates);
  }

  /// Clear per-site rates, reverting to uniform scalar mu.
  pub fn clear_site_rates(&mut self) {
    self.site_rates = None;
  }

  /// Compute P(t) for a single site with a custom rate multiplier.
  ///
  /// Returns exp(Q * mu * rate * t). The eigendecomposition is shared;
  /// only the eigenvalue scaling changes.
  pub fn expQt_with_rate(&self, t: f64, rate: f64) -> Array2<f64> {
    let exp_lt = self.exp_lt_scaled(t, rate);
    let scaled_v_inv = &self.v_inv * &exp_lt.view().insert_axis(Axis(1));
    let Qt = self.v.dot(&scaled_v_inv);
    clamp_min(&Qt, 0.0)
  }

  /// Evolve a sequence profile forward in time (parent -> child).
  ///
  /// Given a probability distribution over states at the parent node, compute the
  /// distribution at a child node after time t of evolution.
  ///
  /// # Arguments
  ///
  /// * `profile` - Parent state probabilities. Shape (L, a) where L is sequence length
  ///   and a is alphabet size. Each row sums to 1 (probability distribution over states).
  /// * `t` - Branch length (time) to evolve forward.
  /// * `return_log` - If true, return log-probabilities.
  ///
  /// # Returns
  ///
  /// Child state probabilities. Shape (L, a), rows sum to 1.
  ///
  /// # Computation
  ///
  /// Uses profile * P(t)' where P(t)' is the transpose of the transition matrix.
  ///
  /// Since P(t) is column-stochastic (P[i,j] = Prob(j -> i)), we need P(t)' to get
  /// row-stochastic behavior for the dot product: result[i] = sum_j profile[j] * P[i,j].
  /// This computes the probability of arriving at state i from any starting state j,
  /// weighted by the parent's probability of being in state j.
  pub fn evolve(&self, profile: &Array2<f64>, t: f64, return_log: bool) -> Array2<f64> {
    let res = match &self.site_rates {
      None => {
        let Qt = self.expQt(t);
        profile.dot(&Qt.t())
      },
      Some(rates) => {
        // Efficient batched per-site computation avoiding L separate matrix multiplies.
        // P_l(t)^T = V_inv^T * diag(exp(λ * μ * r_l * t)) * V^T
        // result_l = profile_l @ P_l^T
        //          = (profile_l @ V_inv^T) .* exp_lt_l @ V^T
        let transformed = profile.dot(&self.v_inv.t());
        let scaled_eigvals = &self.eigvals * (self.mu * t);
        // exp_lt[l, k] = exp(rates[l] * scaled_eigvals[k]) via broadcasting
        let exp_lt = (&rates.view().insert_axis(Axis(1)) * &scaled_eigvals.view().insert_axis(Axis(0))).mapv(f64::exp);
        let scaled = &transformed * &exp_lt;
        clamp_min(&scaled.dot(&self.v.t()), 0.0)
      },
    };
    if return_log { res.mapv(f64::ln) } else { res }
  }

  /// Propagate a sequence profile backward in time (child -> parent).
  ///
  /// Given observed state probabilities at a child node, compute the likelihood
  /// contribution to the parent node. This is the "upward" pass in Felsenstein's
  /// pruning algorithm.
  ///
  /// # Arguments
  ///
  /// * `profile` - Child state likelihoods. Shape (L, a). These are partial likelihoods,
  ///   not necessarily normalized to sum to 1.
  /// * `t` - Branch length (time) from parent to child.
  /// * `return_log` - If true, return log-likelihoods.
  ///
  /// # Returns
  ///
  /// Parent partial likelihoods. Shape (L, a).
  ///
  /// # Computation
  ///
  /// Uses profile * P(t) directly (no transpose).
  ///
  /// Since P(t) is column-stochastic (P[i,j] = Prob(j -> i)), the dot product
  /// result[j] = sum_i profile[i] * P[i,j] computes: for each parent state j,
  /// sum over all child states i of (likelihood of observing child in state i) *
  /// (probability that parent j produces child i).
  ///
  /// This is the likelihood of the observed data given parent state j.
  pub fn propagate_profile(&self, profile: &Array2<f64>, t: f64, return_log: bool) -> Array2<f64> {
    let res = match &self.site_rates {
      None => {
        let Qt = self.expQt(t);
        profile.dot(&Qt)
      },
      Some(rates) => {
        // Efficient batched per-site computation avoiding L separate matrix multiplies.
        // P_l(t) = V * diag(exp(λ * μ * r_l * t)) * V_inv
        // result_l = profile_l @ P_l
        //          = (profile_l @ V) .* exp_lt_l @ V_inv
        let transformed = profile.dot(&self.v);
        let scaled_eigvals = &self.eigvals * (self.mu * t);
        let exp_lt = (&rates.view().insert_axis(Axis(1)) * &scaled_eigvals.view().insert_axis(Axis(0))).mapv(f64::exp);
        let scaled = &transformed * &exp_lt;
        clamp_min(&scaled.dot(&self.v_inv), 0.0)
      },
    };
    if return_log { res.mapv(f64::ln) } else { res }
  }

  /// Compute the transition probability matrix P(t) = exp(Q * mu * t).
  ///
  /// Returns a column-stochastic matrix where P[i,j] = Prob(j -> i | time t).
  /// Column j sums to 1, representing all possible destinations from state j.
  ///
  /// # Computation
  ///
  /// Uses precomputed eigendecomposition for efficiency:
  ///
  ///   P(t) = v * diag(exp(eigvals * mu * t)) * v_inv
  ///
  /// This avoids expensive matrix exponentiation on each call. The eigenvalues are
  /// all non-positive (rates of decay to equilibrium), so exp(eigval * t) is bounded
  /// in [0, 1] and numerically stable.
  ///
  /// # Boundary Behavior
  ///
  /// * t = 0: P(0) = I (identity matrix, no evolution)
  /// * t -> infinity: P[i,j] -> pi[i] for all j (convergence to equilibrium)
  ///
  /// Results are clamped to [0, infinity) to handle floating-point errors that might
  /// produce tiny negative values.
  pub fn expQt(&self, t: f64) -> Array2<f64> {
    // P(t) = V * diag(exp(Lambda * t)) * V_inv
    // Row-scale V_inv by the diagonal vector instead of constructing a full
    // diagonal matrix: O(n^2) broadcast vs O(n^3) matrix multiply.
    let exp_lt = self.exp_lt(t);
    let scaled_v_inv = &self.v_inv * &exp_lt.view().insert_axis(Axis(1));
    let Qt = self.v.dot(&scaled_v_inv);

    // Clamp to handle floating-point errors (should never be significantly negative)
    clamp_min(&Qt, 0.0)
  }

  /// Compute exp(eigenvalue * mu * t) for each eigenvalue.
  fn exp_lt(&self, t: f64) -> Array1<f64> {
    (self.mu * t * &self.eigvals).mapv(f64::exp)
  }

  /// Compute exp(eigenvalue * mu * rate * t) for each eigenvalue, with a custom rate multiplier.
  fn exp_lt_scaled(&self, t: f64, rate: f64) -> Array1<f64> {
    (self.mu * rate * t * &self.eigvals).mapv(f64::exp)
  }

  /// Exponentiated eigenvalues in branch-length space: `exp(eigvals * d)`.
  ///
  /// Branch length `d` is in substitutions per site (= `mu * t`). The `mu`
  /// scaling is already absorbed into `d`, so this method uses raw eigenvalues
  /// without the `mu` factor. This is the correct formulation for the per-edge
  /// likelihood in branch-length optimization.
  ///
  /// Contrast with `exp_lt(t)` which operates in time space and includes `mu`:
  /// `exp(eigvals * mu * t)`. The two are equivalent because `d = mu * t`.
  pub fn exp_eigvals_branch_length(&self, branch_length: f64) -> Array1<f64> {
    (&self.eigvals * branch_length).mapv(f64::exp)
  }

  /// Construct the rate matrix Q in column-stochastic form for display.
  ///
  /// Returns Q where Q[i,j] = W[i,j] * pi[i] for i != j, and columns sum to 0.
  /// This is the standard textbook convention where Q[i,j] represents the rate
  /// of transition INTO state i FROM state j.
  ///
  /// # Convention Note
  ///
  /// This method returns a **different Q** than what `eig_single_site()` uses internally.
  /// The eigendecomposition uses row-stochastic Q (Q[i,j] = W[i,j] * pi[j], rows sum to 0)
  /// because it yields a cleaner similarity transform for symmetrization.
  ///
  /// Both conventions are mathematically valid and satisfy detailed balance. The
  /// eigendecomposition's v and v_inv matrices handle the conversion so that `expQt()`
  /// produces column-stochastic P(t) regardless of internal convention.
  ///
  /// This method exists for debugging, display, and verifying detailed balance in tests.
  /// It is NOT used in the actual evolution computations.
  pub fn Q(&self) -> Array2<f64> {
    // (W * pi) broadcasts pi across columns: result[i,j] = W[i,j] * pi[j]
    // Transpose gives: Q[i,j] = W[j,i] * pi[i] = W[i,j] * pi[i] (W is symmetric)
    let mut Q = (&self.W * &self.pi).t().to_owned();
    // Set diagonal so columns sum to 0
    let diag = -Q.sum_axis(Axis(0));
    Q.diag_mut().assign(&diag);
    Q
  }
}
