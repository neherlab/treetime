#![allow(non_snake_case)]

use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::profile_map::ProfileMap;
use crate::make_report;
use crate::utils::ndarray::outer;
use eyre::Report;
use ndarray::prelude::*;
use ndarray::Zip;
use ndarray_einsum_beta::einsum;
use ndarray_linalg::Eigh;
use ndarray_linalg::UPLO::Lower;
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

/// Performs eigendecompositon of the rate matrix and stores the left- and right-
/// matrices to convert the sequence profiles to the GTR matrix eigenspace
/// and hence to speed-up the computations.
/// NOTE: this assumes the diagonal of W is all zeros
fn eig_single_site(W: &Array2<f32>, pi: &Array1<f32>) -> Result<(Array1<f32>, Array1<f32>, Array2<f32>), Report> {
  assert!(abs(W.diag().sum()) < 1e-10);

  let tmp_pi = pi.mapv(|x| x.sqrt());
  let mut sym_Q: Array2<f32> = W * outer(&tmp_pi, &tmp_pi)?;

  let diag = -(W * pi).sum_axis(Axis(1));
  sym_Q.diag_mut().assign(&diag);

  let (eigvals, eigvecs) = sym_Q.eigh(Lower)?;
  let tmp_v = eigvecs.t().dot(&tmp_pi);
  let one_norm: f32 = tmp_v.mapv(|x| x.abs()).sum();

  let v = tmp_v / one_norm;
  let v_inv = (eigvecs * one_norm) / tmp_pi;

  Ok((eigvals, v, v_inv))
}

pub struct GTRParams {
  pub alphabet: Alphabet,
  pub profile_map: ProfileMap,
  pub mu: f32,
  pub W: Array2<f32>,
  pub pi: Array1<f32>,
}

/// Defines General-Time-Reversible model of character evolution.
pub struct GTR {
  pub debug: bool,
  pub is_site_specific: bool,
  pub alphabet: Alphabet,
  pub profile_map: ProfileMap,
  pub mu: f32,
  pub W: Array2<f32>,
  pub pi: Array1<f32>,
  pub eigvals: Array1<f32>,
  pub v: Array1<f32>,
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
    let Qt = self.expQt(t).t();
    let res = profile.dot(&Qt);
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
  pub fn propagate_profile(self, profile: &Array2<f32>, t: f32, return_log: bool) -> Array2<f32> {
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
    let eLambdaT = self.exp_lt(t).diag(); // vector length = a

    let dot1 = eLambdaT.dot(&self.v_inv);
    let Qs: Array2<f32> = self.v.dot(&dot1); // This is P(nuc1 | given nuc_2)

    // np.maximum(0, Qs)
    let zeros: Array2<f32> = Array2::<f32>::zeros(Qs.dim());
    Zip::from(&zeros).and(&Qs).map_collect(|a, &b| a.max(b))
  }

  fn exp_lt(&self, t: f32) -> Array1<f32> {
    (self.mu * t * &self.eigvals).mapv(f32::exp)
  }

  //     @classmethod
  //     def custom(cls, mu=1.0, pi=None, W=None, **kwargs):
  //         """
  //         Create a GTR model by specifying the matrix explicitly
  //
  //         Parameters
  //         ----------
  //
  //          mu : float
  //             Substitution rate
  //
  //          W : nxn matrix
  //             Substitution matrix
  //
  //          pi : n vector
  //             Equilibrium frequencies
  //
  //          **kwargs:
  //             Key word arguments to be passed
  //
  //         Keyword Args
  //         ------------
  //
  //          alphabet : str
  //             Specify alphabet when applicable. If the alphabet specification is
  //             required, but no alphabet is specified, the nucleotide alphabet will be used as
  //             default.
  //
  //         """
  //         gtr = cls(**kwargs)
  //         gtr.assign_rates(mu=mu, pi=pi, W=W)
  //         return gtr
  //
  //     @staticmethod
  //     def standard(model, **kwargs):
  //         """
  //         Create standard model of molecular evolution.
  //
  //         Parameters
  //         ----------
  //
  //          model : str
  //             Model to create. See list of available models below
  //          **kwargs:
  //             Key word arguments to be passed to the model
  //
  //
  //         **Available models**
  //
  //         - JC69:
  //
  //           Jukes-Cantor 1969 model. This model assumes equal frequencies
  //           of the nucleotides and equal transition rates between nucleotide states.
  //           For more info, see: Jukes and Cantor (1969).
  //           Evolution of Protein Molecules. New York: Academic Press. pp. 21-132.
  //           To create this model, use:
  //
  //           :code:`mygtr = GTR.standard(model='jc69', mu=<my_mu>, alphabet=<my_alph>)`
  //
  //           :code:`my_mu` - substitution rate (float)
  //
  //           :code:`my_alph` - alphabet (str: :code:`'nuc'` or  :code:`'nuc_nogap'`)
  //
  //
  //
  //         - K80:
  //
  //           Kimura 1980 model. Assumes equal concentrations across nucleotides, but
  //           allows different rates between transitions and transversions. The ratio
  //           of the transversion/transition rates is given by kappa parameter.
  //           For more info, see
  //           Kimura (1980),  J. Mol. Evol. 16 (2): 111-120. doi:10.1007/BF01731581.
  //           Current implementation of the model does not account for the gaps.
  //
  //
  //           :code:`mygtr = GTR.standard(model='k80', mu=<my_mu>, kappa=<my_kappa>)`
  //
  //           :code:`mu` - overall substitution rate (float)
  //
  //           :code:`kappa` - ratio of transversion/transition rates (float)
  //
  //
  //         - F81:
  //
  //           Felsenstein 1981 model. Assumes non-equal concentrations across nucleotides,
  //           but the transition rate between all states is assumed to be equal. See
  //           Felsenstein (1981), J. Mol. Evol. 17  (6): 368-376. doi:10.1007/BF01734359
  //           for details.
  //
  //           :code:`mygtr = GTR.standard(model='F81', mu=<mu>, pi=<pi>, alphabet=<alph>)`
  //
  //           :code:`mu` -  substitution rate  (float)
  //
  //           :code:`pi`  - : nucleotide concentrations (numpy.array)
  //
  //           :code:`alphabet' -  alphabet to use. (:code:`'nuc'` or  :code:`'nuc_nogap'`)
  //
  //
  //         - HKY85:
  //
  //           Hasegawa, Kishino and Yano 1985 model. Allows different concentrations of the
  //           nucleotides (as in F81) + distinguishes between transition/transversion substitutions
  //           (similar to K80). Link:
  //           Hasegawa, Kishino, Yano (1985), J. Mol. Evol. 22 (2): 160-174. doi:10.1007/BF02101694
  //
  //           Current implementation of the model does not account for the gaps
  //
  //           :code:`mygtr = GTR.standard(model='HKY85', mu=<mu>, pi=<pi>, kappa=<kappa>)`
  //
  //           :code:`mu` -  substitution rate  (float)
  //
  //           :code:`pi`  - : nucleotide concentrations (numpy.array)
  //
  //           :code:`kappa` - ratio of transversion/transition rates (float)
  //
  //
  //
  //         - T92:
  //
  //           Tamura 1992 model. Extending Kimura  (1980) model for the case where a
  //           G+C-content bias exists. Link:
  //           Tamura K (1992),  Mol.  Biol. Evol. 9 (4): 678-687.  DOI: 10.1093/oxfordjournals.molbev.a040752
  //
  //           Current implementation of the model does not account for the gaps
  //
  //           :code:`mygtr = GTR.standard(model='T92', mu=<mu>, pi_GC=<pi_gc>, kappa=<kappa>)`
  //
  //           :code:`mu` -  substitution rate  (float)
  //
  //           :code:`pi_GC`  - : relative GC content
  //
  //           :code:`kappa` - ratio of transversion/transition rates (float)
  //
  //
  //         - TN93:
  //
  //           Tamura and Nei 1993. The model distinguishes between the two different types of
  //           transition: (A <-> G) is allowed to have a different rate to (C<->T).
  //           Transversions have the same rate. The frequencies of the nucleotides are allowed
  //           to be different. Link: Tamura, Nei (1993), MolBiol Evol. 10 (3): 512-526.
  //           DOI:10.1093/oxfordjournals.molbev.a040023
  //
  //           :code:`mygtr = GTR.standard(model='TN93', mu=<mu>, kappa1=<k1>, kappa2=<k2>)`
  //
  //           :code:`mu` -  substitution rate  (float)
  //
  //           :code:`kappa1`  - relative A<-->C, A<-->T, T<-->G and G<-->C rates (float)
  //
  //           :code:`kappa` - relative C<-->T rate (float)
  //
  //           .. Note::
  //                 Rate of A<-->G substitution is set to one. All other rates
  //                 (kappa1, kappa2) are specified relative to this rate
  //
  //         """
  //         from .nuc_models import JC69, K80, F81, HKY85, T92, TN93
  //         from .aa_models  import JTT92
  //
  //         if model.lower() in ['jc', 'jc69', 'jukes-cantor', 'jukes-cantor69', 'jukescantor', 'jukescantor69']:
  //             model = JC69(**kwargs)
  //         elif model.lower() in ['k80', 'kimura80', 'kimura1980']:
  //             model = K80(**kwargs)
  //         elif model.lower() in ['f81', 'felsenstein81', 'felsenstein1981']:
  //             model = F81(**kwargs)
  //         elif model.lower() in ['hky', 'hky85', 'hky1985']:
  //             model = HKY85(**kwargs)
  //         elif model.lower() in ['t92', 'tamura92', 'tamura1992']:
  //             model = T92(**kwargs)
  //         elif model.lower() in ['tn93', 'tamura_nei_93', 'tamuranei93']:
  //             model = TN93(**kwargs)
  //         elif model.lower() in ['jtt', 'jtt92']:
  //             model = JTT92(**kwargs)
  //         else:
  //             raise KeyError("The GTR model '{}' is not in the list of available models."
  //                 "".format(model))
  //
  //         model.mu = kwargs['mu'] if 'mu' in kwargs else 1.0
  //         return model
  //
  //
  //     @classmethod
  //     def random(cls, mu=1.0, alphabet='nuc'):
  //         """
  //         Creates a random GTR model
  //
  //         Parameters
  //         ----------
  //
  //          mu : float
  //             Substitution rate
  //
  //          alphabet : str
  //             Alphabet name (should be standard: 'nuc', 'nuc_gap', 'aa', 'aa_gap')
  //
  //
  //         """
  //
  //         alphabet=alphabets[alphabet]
  //         gtr = cls(alphabet)
  //         n = gtr.alphabet.shape[0]
  //         pi = 1.0*np.random.randint(0,100,size=(n))
  //         W = 1.0*np.random.randint(0,100,size=(n,n)) // with gaps
  //
  //         gtr.assign_rates(mu=mu, pi=pi, W=W)
  //         return gtr
  //
  //
  //     @classmethod
  //     def infer(cls, nij, Ti, root_state, fixed_pi=None, pc=1.0, gap_limit=0.01, **kwargs):
  //         r"""
  //         Infer a GTR model by specifying the number of transitions and time spent in each
  //         character. The basic equation that is being solved is
  //
  //         :math:`n_{ij} = pi_i W_{ij} T_j`
  //
  //         where :math:`n_{ij}` are the transitions, :math:`pi_i` are the equilibrium
  //         state frequencies, :math:`W_{ij}` is the "substitution attempt matrix",
  //         while :math:`T_i` is the time on the tree spent in character state
  //         :math:`i`. To regularize the process, we add pseudocounts and also need
  //         to account for the fact that the root of the tree is in a particular
  //         state. the modified equation is
  //
  //         :math:`n_{ij} + pc = pi_i W_{ij} (T_j+pc+root\_state)`
  //
  //         Parameters
  //         ----------
  //
  //          nij : nxn matrix
  //             The number of times a change in character state is observed
  //             between state j and i
  //
  //          Ti :n vector
  //             The time spent in each character state
  //
  //          root_state : n vector
  //             The number of characters in state i in the sequence
  //             of the root node.
  //
  //          pc : float
  //             Pseudocounts, this determines the lower cutoff on the rate when
  //             no substitutions are observed
  //
  //          **kwargs:
  //             Key word arguments to be passed
  //
  //         Keyword Args
  //         ------------
  //
  //          alphabet : str
  //             Specify alphabet when applicable. If the alphabet specification
  //             is required, but no alphabet is specified, the nucleotide alphabet will be used as default.
  //
  //         """
  //         from scipy import linalg as LA
  //         gtr = cls(**kwargs)
  //         dp = 1e-5
  //         Nit = 40
  //         pc_mat = pc*np.ones_like(nij)
  //         np.fill_diagonal(pc_mat, 0.0)
  //         np.fill_diagonal(nij, 0.0)
  //         count = 0
  //         pi_old = np.zeros_like(Ti)
  //         if fixed_pi is None:
  //             pi = np.ones_like(Ti)
  //         else:
  //             pi = np.copy(fixed_pi)
  //         pi/=pi.sum()
  //         W_ij = np.ones_like(nij)
  //         mu = nij.sum()/Ti.sum()
  //         // if pi is fixed, this will immediately converge
  //         while LA.norm(pi_old-pi) > dp and count < Nit:
  //             count += 1
  //             pi_old = np.copy(pi)
  //             W_ij = (nij+nij.T+2*pc_mat)/mu/(np.outer(pi,Ti) + np.outer(Ti,pi)
  //                                                     + ttconf.TINY_NUMBER + 2*pc_mat)
  //
  //             np.fill_diagonal(W_ij, 0)
  //             scale_factor = avg_transition(W_ij,pi, gap_index=gtr.gap_index)
  //
  //             W_ij = W_ij/scale_factor
  //             if fixed_pi is None:
  //                 pi = (np.sum(nij+pc_mat,axis=1)+root_state)/(ttconf.TINY_NUMBER + mu*np.dot(W_ij,Ti)+root_state.sum()+np.sum(pc_mat, axis=1))
  //                 pi /= pi.sum()
  //                 mu = nij.sum()/(ttconf.TINY_NUMBER + np.sum(pi * (W_ij.dot(Ti))))
  //             else:
  //                 mu = nij.sum()/(ttconf.TINY_NUMBER + np.sum(pi * (W_ij.dot(pi)))*Ti.sum())
  //
  //         if count >= Nit:
  //             if LA.norm(pi_old-pi) > dp:
  //             elif np.abs(1-np.max(pi.sum(axis=0))) > dp:
  //         if gtr.gap_index is not None:
  //             if pi[gtr.gap_index]<gap_limit:
  //                        ' this can potentially result in artificats.'+
  //                        ' gap fraction will be set to %1.4f'%gap_limit,2,warn=true)
  //             pi[gtr.gap_index] = gap_limit
  //             pi /= pi.sum()
  //
  //         gtr.assign_rates(mu=mu, W=W_ij, pi=pi)
  //         return gtr
  //
  // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ////// prepare model
  // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //
  //     def _eig_single_site(self, W, p):
  //         """
  //         Perform eigendecompositon of the rate matrix and stores the left- and right-
  //         matrices to convert the sequence profiles to the GTR matrix eigenspace
  //         and hence to speed-up the computations.
  //         NOTE: this assumes the diagonal of W is all zeros
  //         """
  //         // eigendecomposition of the rate matrix
  //         assert np.abs(np.diag(W).sum())<1e-10
  //
  //         tmpp = np.sqrt(p)
  //         symQ = W*np.outer(tmpp, tmpp)
  //         np.fill_diagonal(symQ, -np.sum(W*p, axis=1))
  //
  //         eigvals, eigvecs = np.linalg.eigh(symQ)
  //         tmp_v = eigvecs.T*tmpp
  //         one_norm = np.sum(np.abs(tmp_v), axis=1)
  //         return eigvals, tmp_v.T/one_norm, (eigvecs*one_norm).T/tmpp
  //
  //
  //     def state_pair(self, seq_p, seq_ch, pattern_multiplicity=None,
  //                                ignore_gaps=false):
  //         '''
  //         Make a compressed representation of a pair of sequences, only counting
  //         the number of times a particular pair of states (e.g. (A,T)) is observed
  //         in the aligned sequences of parent and child.
  //
  //         Parameters
  //         ----------
  //
  //          seq_p: numpy array
  //             Parent sequence as numpy array of chars
  //
  //          seq_ch: numpy array
  //             Child sequence as numpy array of chars
  //
  //          pattern_multiplicity : numpy array
  //             If sequences are reduced by combining identical alignment patterns,
  //             these multplicities need to be accounted for when counting the number
  //             of mutations across a branch. If None, all pattern are assumed to
  //             occur exactly once.
  //
  //          ignore_gap: bool
  //             Whether or not to include gapped positions of the alignment
  //             in the multiplicity count
  //
  //         Returns
  //         -------
  //           seq_pair : list
  //             :code:`[(0,1), (2,2), (3,4)]` list of parent_child state pairs
  //             as indices in the alphabet
  //
  //           multiplicity : numpy array
  //             Number of times a particular pair is observed
  //
  //         '''
  //         if pattern_multiplicity is None:
  //             pattern_multiplicity = np.ones_like(seq_p, dtype=float)
  //
  //         from collections import Counter
  //         if seq_ch.shape != seq_p.shape:
  //             raise ValueError("GTR.state_pair: Sequence lengths do not match!")
  //
  //         if len(self.alphabet)<10: // for small alphabet, repeatedly check array for all state pairs
  //             pair_count = []
  //             bool_seqs_p = []
  //             bool_seqs_ch = []
  //             for seq, bs in [(seq_p,bool_seqs_p), (seq_ch, bool_seqs_ch)]:
  //                 for ni,nuc in enumerate(self.alphabet):
  //                     bs.append(seq==nuc)
  //
  //             for n1,nuc1 in enumerate(self.alphabet):
  //                 if (self.gap_index is None) or (not ignore_gaps) or (n1!=self.gap_index):
  //                     for n2,nuc2 in enumerate(self.alphabet):
  //                         if (self.gap_index is None) or (not ignore_gaps) or (n2!=self.gap_index):
  //                             count = ((bool_seqs_p[n1]&bool_seqs_ch[n2])*pattern_multiplicity).sum()
  //                             if count: pair_count.append(((n1,n2), count))
  //         else: // enumerate state pairs of the sequence for large alphabets
  //             num_seqs = []
  //             for seq in [seq_p, seq_ch]: // for each sequence (parent and child) construct a numerical sequence [0,5,3,1,2,3...]
  //                 tmp = np.ones_like(seq, dtype=int)
  //                 for ni,nuc in enumerate(self.alphabet):
  //                     tmp[seq==nuc] = ni  // set each position corresponding to a state to the corresponding index
  //                 num_seqs.append(tmp)
  //             pair_count = defaultdict(int)
  //             if ignore_gaps:  // if gaps are ignored skip positions where one or the other sequence is gapped
  //                 for i in range(len(seq_p)):
  //                     if self.gap_index!=num_seqs[0][i] and self.gap_index!=num_seqs[1][i]:
  //                         pair_count[(num_seqs[0][i],num_seqs[1][i])]+=pattern_multiplicity[i]
  //             else: // otherwise, just count
  //                 for i in range(len(seq_p)):
  //                     pair_count[(num_seqs[0][i],num_seqs[1][i])]+=pattern_multiplicity[i]
  //             pair_count = pair_count.items()
  //
  //         return (np.array([x[0] for x in pair_count], dtype=int),    // [(child_nuc, parent_nuc),()...]
  //                 np.array([x[1] for x in pair_count], dtype=int))    // multiplicity of each parent/child nuc pair
  //
  //
  // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ////// evolution functions
  // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //     def prob_t_compressed(self, seq_pair, multiplicity, t, return_log=false):
  //         '''
  //         Calculate the probability of observing a sequence pair at a distance t,
  //         for compressed sequences
  //
  //         Parameters
  //         ----------
  //
  //           seq_pair : numpy array
  //             :code:`np.array([(0,1), (2,2), ()..])` as indicies of
  //             pairs of aligned positions. (e.g. 'A'==0, 'C'==1 etc).
  //             This only lists all occuring parent-child state pairs, order is irrelevant
  //
  //           multiplicity : numpy array
  //             The number of times a parent-child state pair is observed.
  //             This allows compression of the sequence representation
  //
  //           t : float
  //             Length of the branch separating parent and child
  //
  //           return_log : bool
  //             Whether or not to exponentiate the result
  //
  //         '''
  //         if t<0:
  //             logP = -ttconf.BIG_NUMBER
  //         else:
  //             tmp_eQT = self.expQt(t)
  //             logQt = np.log(np.maximum(tmp_eQT, ttconf.SUPERTINY_NUMBER))
  //             // logQt[~np.isfinite(logQt)] = -ttconf.BIG_NUMBER
  //             logP = np.sum(logQt[seq_pair[:,1], seq_pair[:,0]]*multiplicity)
  //
  //         return logP if return_log else np.exp(logP)
  //
  //
  //     def prob_t(self, seq_p, seq_ch, t, pattern_multiplicity = None,
  //                return_log=false, ignore_gaps=true):
  //         """
  //         Compute the probability to observe seq_ch (child sequence) after time t starting from seq_p
  //         (parent sequence).
  //
  //         Parameters
  //         ----------
  //
  //          seq_p : character array
  //             Parent sequence
  //
  //          seq_c : character array
  //             Child sequence
  //
  //          t : double
  //             Time (branch len) separating the profiles.
  //
  //          pattern_multiplicity : numpy array
  //             If sequences are reduced by combining identical alignment patterns,
  //             these multplicities need to be accounted for when counting the number
  //             of mutations across a branch. If None, all pattern are assumed to
  //             occur exactly once.
  //
  //          return_log : bool
  //             It true, return log-probability.
  //
  //         Returns
  //         -------
  //          prob : np.array
  //             Resulting probability
  //
  //         """
  //         seq_pair, multiplicity = self.state_pair(seq_p, seq_ch,
  //                                         pattern_multiplicity=pattern_multiplicity, ignore_gaps=ignore_gaps)
  //         return self.prob_t_compressed(seq_pair, multiplicity, t, return_log=return_log)
  //
  //
  //     def optimal_t(self, seq_p, seq_ch, pattern_multiplicity=None, ignore_gaps=false):
  //         '''
  //         Find the optimal distance between the two sequences
  //
  //         Parameters
  //         ----------
  //
  //          seq_p : character array
  //             Parent sequence
  //
  //          seq_c : character array
  //             Child sequence
  //
  //          pattern_multiplicity : numpy array
  //             If sequences are reduced by combining identical alignment patterns,
  //             these multplicities need to be accounted for when counting the number
  //             of mutations across a branch. If None, all pattern are assumed to
  //             occur exactly once.
  //
  //          ignore_gaps : bool
  //             If true, ignore gaps in distance calculations
  //
  //         '''
  //         seq_pair, multiplicity = self.state_pair(seq_p, seq_ch,
  //                                         pattern_multiplicity = pattern_multiplicity,
  //                                         ignore_gaps=ignore_gaps)
  //         return self.optimal_t_compressed(seq_pair, multiplicity)
  //
  //
  //     def optimal_t_compressed(self, seq_pair, multiplicity, profiles=false, tol=1e-10):
  //         """
  //         Find the optimal distance between the two sequences represented as state_pairs
  //         or as pair of profiles
  //
  //         Parameters
  //         ----------
  //
  //          seq_pair : state_pair, tuple
  //             Compressed representation of sequences along a branch, either
  //             as tuple of state pairs or as tuple of profiles.
  //
  //          multiplicity : array
  //             Number of times each state pair in seq_pair appears (if profile==false)
  //
  //             Number of times an alignment pattern is observed (if profiles==true)
  //
  //          profiles : bool, default false
  //             The standard branch length optimization assumes fixed sequences at
  //             either end of the branch. With profiles==true, optimization is performed
  //             while summing over all possible states of the nodes at either end of the
  //             branch. Note that the meaning/format of seq_pair and multiplicity
  //             depend on the value of :profiles:.
  //
  //         """
  //
  //         def _neg_prob(t, seq_pair, multiplicity):
  //             """
  //             Probability to observe a child given the the parent state, transition
  //             matrix, and the time of evolution (branch length).
  //
  //             Parameters
  //             ----------
  //
  //              t : double
  //                 Branch length (time between sequences)
  //
  //              seq_pair : tuple of profiles
  //                 parent and child sequences
  //
  //              multiplicity : vector containing the number of times each alignment pattern is observed
  //
  //             Returns
  //             -------
  //
  //              prob : double
  //                 Negative probability of the two given sequences
  //                 to be separated by the time t.
  //             """
  //             if profiles:
  //                 res = -1.0*self.prob_t_profiles(seq_pair, multiplicity,t**2, return_log=true)
  //                 return res + np.exp(t**4/10000)
  //             else:
  //                 return -1.0*self.prob_t_compressed(seq_pair, multiplicity,t**2, return_log=true)
  //
  //         try:
  //             from scipy.optimize import minimize_scalar
  //             opt = minimize_scalar(_neg_prob,
  //                     bounds=[-np.sqrt(ttconf.MAX_BRANCH_LENGTH),np.sqrt(ttconf.MAX_BRANCH_LENGTH)],
  //                     args=(seq_pair, multiplicity), tol=tol)
  //             new_len = opt["x"]**2
  //             if 'success' not in opt:
  //                 opt['success'] = true
  //         except ImportError:
  //             import scipy
  //             print('legacy scipy', scipy.__version__)
  //             from scipy.optimize import fminbound
  //             new_len = fminbound(_neg_prob,
  //                       -np.sqrt(ttconf.MAX_BRANCH_LENGTH),np.sqrt(ttconf.MAX_BRANCH_LENGTH),
  //                        args=(seq_pair, multiplicity))
  //             new_len = new_len**2
  //             opt={'success':true}
  //
  //         if new_len > .9 * ttconf.MAX_BRANCH_LENGTH:
  //
  //         if opt["success"] != true:
  //             // return hamming distance: number of state pairs where state differs/all pairs
  //             new_len =  np.sum(multiplicity[seq_pair[:,1]!=seq_pair[:,0]])/np.sum(multiplicity)
  //
  //         return new_len
  //
  //
  //     def prob_t_profiles(self, profile_pair, multiplicity, t,
  //                         return_log=false, ignore_gaps=true):
  //         '''
  //         Calculate the probability of observing a node pair at a distance t
  //
  //         Parameters
  //         ----------
  //
  //           profile_pair: numpy arrays
  //             Probability distributions of the nucleotides at either
  //             end of the branch. pp[0] = parent, pp[1] = child
  //
  //           multiplicity : numpy array
  //             The number of times an alignment pattern is observed
  //
  //           t : float
  //             Length of the branch separating parent and child
  //
  //           ignore_gaps: bool
  //             If true, ignore mutations to and from gaps in distance calculations
  //
  //           return_log : bool
  //             Whether or not to exponentiate the result
  //
  //         '''
  //         if t<0:
  //             logP = -ttconf.BIG_NUMBER
  //         else:
  //             Qt = self.expQt(t)
  //             if len(Qt.shape)==3: // site specific GTR model
  //                 res = np.einsum('ai,ija,aj->a', profile_pair[1], Qt, profile_pair[0])
  //             else:
  //                 res = np.einsum('ai,ij,aj->a', profile_pair[1], Qt, profile_pair[0])
  //             if ignore_gaps and (self.gap_index is not None): // calculate the probability that neither outgroup/node has a gap
  //                 non_gap_frac = (1-profile_pair[0][:,self.gap_index])*(1-profile_pair[1][:,self.gap_index])
  //                 // weigh log LH by the non-gap probability
  //                 logP = np.sum(multiplicity*np.log(res+ttconf.SUPERTINY_NUMBER)*non_gap_frac)
  //             else:
  //                 logP = np.sum(multiplicity*np.log(res+ttconf.SUPERTINY_NUMBER))
  //
  //         return logP if return_log else np.exp(logP)
  //
  //
  //     def propagate_profile(self, profile, t, return_log=false):
  //         """
  //         Compute the probability of the sequence state of the parent
  //         at time (t+t0, backwards), given the sequence state of the
  //         child (profile) at time t0.
  //
  //         Parameters
  //         ----------
  //
  //          profile : numpy.array
  //             Sequence profile. Shape = (L, a),
  //             where L - sequence length, a - alphabet size.
  //
  //          t : double
  //             Time to propagate
  //
  //          return_log: bool
  //             If true, return log-probability
  //
  //         Returns
  //         -------
  //
  //          res : np.array
  //             Profile of the sequence after time t in the past.
  //             Shape = (L, a), where L - sequence length, a - alphabet size.
  //
  //         """
  //         Qt = self.expQt(t)
  //         res = profile.dot(Qt)
  //
  //         return np.log(res) if return_log else res
  //
  //

  //
  //
  //     def expQs(self, s):
  //         return self.expQt(s**2)
  //
  //
  //     def expQsds(self, s):
  //         r'''
  //         Returns
  //         -------
  //         Qtds :  Returns 2 V_{ij} \lambda_j s e^{\lambda_j s**2 } V^{-1}_{jk}
  //                 This is the derivative of the branch probability with respect to s=\sqrt(t)
  //         '''
  //         lambda_eLambdaT = np.diag(2.0*self._exp_lt(s**2)*self.eigenvals*s)
  //         return self.v.dot(lambda_eLambdaT.dot(self.v_inv))
  //
  //
  //     def sequence_logLH(self,seq, pattern_multiplicity=None):
  //         """
  //         Returns the log-likelihood of sampling a sequence from equilibrium frequency.
  //         Expects a sequence as numpy array
  //
  //         Parameters
  //         ----------
  //
  //          seq : numpy array
  //             Compressed sequence as an array of chars
  //
  //          pattern_multiplicity : numpy_array
  //             The number of times each position in sequence is observed in the
  //             initial alignment. If None, sequence is assumed to be not compressed
  //
  //         """
  //         if pattern_multiplicity is None:
  //             pattern_multiplicity = np.ones_like(seq, dtype=float)
  //         return np.sum([np.sum((seq==state)*pattern_multiplicity*np.log(self.Pi[si]))
  //                       for si,state in enumerate(self.alphabet)])
  //
  //     def average_rate(self):
  //         return self.mu*avg_transition(self.W, self.Pi, gap_index=self.gap_index)
  //
  //     def save_to_npz(self, outfile):
  //         full_gtr = self.mu * np.dot(self.Pi, self.W)
  //         desc=np.array(["GTR matrix description\n", "Substitution rate: " + str(self.mu)])
  //         np.savez(outfile,   description=desc,
  //                             full_gtr=full_gtr,
  //                             char_dist=self.Pi,
  //                             flow_matrix=self.W)
  //
  //
  // if __name__ == "__main__":
  //     pass

  //     @property
  //     def mu(self):
  //         return self._mu
  //
  //     @property
  //     def Pi(self):
  //         return self._Pi
  //
  //     @property
  //     def W(self):
  //         return self._W
  //
  //     @W.setter
  //     def W(self, value):
  //         self.assign_rates(mu=self.mu, pi=self.Pi, W=value)
  //
  //     @Pi.setter
  //     def Pi(self, value):
  //         self.assign_rates(mu=self.mu, pi=value, W=self.W)
  //
  //     @mu.setter
  //     def mu(self, value):
  //         self.assign_rates(mu=value, pi=self.Pi, W=self.W)

  //
  //     @property
  //     def Q(self):
  //         """function that return the product of the transition matrix
  //            and the equilibrium frequencies to obtain the rate matrix
  //            of the GTR model
  //         """
  //         Q_tmp = (self.W*self.Pi).T
  //         Q_diag = -np.sum(Q_tmp, axis=0)
  //         np.fill_diagonal(Q_tmp, Q_diag)
  //         return Q_tmp

  // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // //// constructor methods
  // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //     def __str__(self):
  //         '''
  //         String representation of the GTR model for pretty printing
  //         '''
  //         multi_site = len(self.Pi.shape)==2
  //
  //         if multi_site:
  //             eq_freq_str = "Average substitution rate (mu): "+str(np.round(self.average_rate,6))+'\n'
  //         else:
  //             eq_freq_str = "Substitution rate (mu): "+str(np.round(self.mu,6))+'\n'
  //
  //         if not multi_site:
  //             eq_freq_str += "\nEquilibrium frequencies (pi_i):\n"
  //             for a,p in zip(self.alphabet, self.Pi):
  //                 eq_freq_str+='  '+a+': '+str(np.round(p,4))+'\n'
  //
  //         W_str = "\nSymmetrized rates from j->i (W_ij):\n"
  //         W_str+='\t'+'\t'.join(self.alphabet)+'\n'
  //         for a,Wi in zip(self.alphabet, self.W):
  //             W_str+= '  '+a+'\t'+'\t'.join([str(np.round(max(0,p),4)) for p in Wi])+'\n'
  //
  //         if not multi_site:
  //             Q_str = "\nActual rates from j->i (Q_ij):\n"
  //             Q_str+='\t'+'\t'.join(self.alphabet)+'\n'
  //             for a,Qi in zip(self.alphabet, self.Q):
  //                 Q_str+= '  '+a+'\t'+'\t'.join([str(np.round(max(0,p),4)) for p in Qi])+'\n'
  //
  //         return eq_freq_str + W_str + Q_str
}

mod test {
  use crate::nuc_models::jc69::jc69;
  use eyre::Report;
  use ndarray::{array, Array2};

  fn evolves() -> Result<(), Report> {
    let gtr = jc69()?;

    let norm_prof: Array2<f32> = array![
      [0.25307157, 0.2191123, 0.01132857, 0.24025543, 0.27623213],
      [0.08326688, 0.11709925, 0.2998784, 0.2882319, 0.21152357],
      [0.11515529, 0.15260233, 0.06456645, 0.31272216, 0.35495376],
      [0.08100234, 0.32020589, 0.12664117, 0.23628513, 0.23586548],
      [0.11945823, 0.3206602, 0.15715599, 0.10981392, 0.29291166],
      [0.14397598, 0.23789498, 0.29021633, 0.21377755, 0.11413517],
      [0.11371981, 0.10308411, 0.30645272, 0.24627502, 0.23046833],
      [0.21878948, 0.41258937, 0.26041167, 0.05728378, 0.0509257,],
      [0.02167063, 0.11335226, 0.33757606, 0.1919524, 0.33544866],
      [0.02839864, 0.20064851, 0.06049814, 0.3828229, 0.32763181]
    ];

    gtr.evolve(&norm_prof, 0.1, false);
    // gtr.expQt(0.1)?;

    let eLambdaT = array![
      [0.8824969, 0., 0., 0., 0.],
      [0., 0.8824969, 0., 0., 0.],
      [0., 0., 0.8824969, 0., 0.],
      [0., 0., 0., 0.8824969, 0.],
      [0., 0., 0., 0., 1.]
    ];

    let dot1 = array![
      [
        -7.80313051e-01,
        -1.26800871e+00,
        6.82773919e-01,
        6.82773919e-01,
        6.82773919e-01
      ],
      [
        0.00000000e+00,
        0.00000000e+00,
        -8.03675458e-01,
        1.09784109e+00,
        -2.94165634e-01
      ],
      [
        0.00000000e+00,
        2.89769610e-17,
        -8.03675458e-01,
        -2.94165634e-01,
        1.09784109e+00
      ],
      [
        -1.11473293e+00,
        8.36049697e-01,
        9.28944108e-02,
        9.28944108e-02,
        9.28944108e-02
      ],
      [
        -1.00000000e+00,
        -1.00000000e+00,
        -1.00000000e+00,
        -1.00000000e+00,
        -1.00000000e+00
      ],
    ];

    let Qs = array![
      [0.90599752, 0.02350062, 0.02350062, 0.02350062, 0.02350062],
      [0.02350062, 0.90599752, 0.02350062, 0.02350062, 0.02350062],
      [0.02350062, 0.02350062, 0.90599752, 0.02350062, 0.02350062],
      [0.02350062, 0.02350062, 0.02350062, 0.90599752, 0.02350062],
      [0.02350062, 0.02350062, 0.02350062, 0.02350062, 0.90599752],
    ];

    Ok(())
  }
}
