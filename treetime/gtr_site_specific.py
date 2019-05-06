from __future__ import division, print_function, absolute_import
from collections import defaultdict
import numpy as np
from treetime import config as ttconf
from .seq_utils import alphabets, profile_maps, alphabet_synonyms
from .gtr import GTR

class GTR_site_specific(GTR):
    """
    Defines General-Time-Reversible model of character evolution that
    allows for different models at different sites in the alignment
    """
    def __init__(self, seq_len=1, approximate=True, **kwargs):
        """constructor for site specfic GTR models

        Parameters
        ----------
        seq_len : int, optional
            number of sites, determines dimensions of frequency vectors etc
        approximate : bool, optional
            use linear interpolation for exponentiated matrices to speed up calcuations
        **kwargs
            Description
        """
        self.seq_len=seq_len
        self.approximate = approximate
        super(GTR_site_specific, self).__init__(**kwargs)


    @property
    def Q(self):
        """function that return the product of the transition matrix
           and the equilibrium frequencies to obtain the rate matrix
           of the GTR model
        """
        tmp = np.einsum('ia,ij->ija', self.Pi, self.W)
        diag_vals = np.sum(tmp, axis=0)
        for x in range(tmp.shape[-1]):
            np.fill_diagonal(tmp[:,:,x], -diag_vals[:,x])
        return tmp


    def assign_rates(self, mu=1.0, pi=None, W=None):
        """
        Overwrite the GTR model given the provided data

        Parameters
        ----------

         mu : float
            Substitution rate

         W : nxn matrix
            Substitution matrix

         pi : n vector
            Equilibrium frequencies

        """
        n = len(self.alphabet)
        if np.isscalar(mu):
            self._mu = mu*np.ones(self.seq_len)
        else:
            self._mu = np.copy(mu)

        if pi is not None and pi.shape[0]==n:
            self.seq_len = pi.shape[-1]
            Pi = np.copy(pi)
        else:
            if pi is not None and len(pi)!=n:
                raise ArgumentError("GTR_site_specific: length of equilibrium frequency vector does not match alphabet length.")
            Pi = np.ones(shape=(n,self.seq_len))

        self._Pi = Pi/np.sum(Pi, axis=0)

        if W is None or W.shape!=(n,n):
            if (W is not None) and W.shape!=(n,n):
                raise ArgumentError("GTR_site_specific: Size of substitution matrix does not match alphabet length.")
            W = np.ones((n,n))
            np.fill_diagonal(W, 0.0)
            np.fill_diagonal(W, - W.sum(axis=0))
        else:
            W=0.5*(np.copy(W)+np.copy(W).T)

        np.fill_diagonal(W,0)
        avg_pi = self.Pi.mean(axis=-1)
        average_rate = W.dot(avg_pi).dot(avg_pi)
        self._W = W/average_rate
        self._mu *=average_rate

        self._eig()
        self._make_expQt_interpolator()


    @classmethod
    def random(cls, L=1, avg_mu=1.0, alphabet='nuc', pi_dirichlet_alpha=1,
               W_dirichlet_alpha=3.0, mu_gamma_alpha=3.0):
        """
        Creates a random GTR model

        Parameters
        ----------
        L : int, optional
            number of sites for which to generate a model
        avg_mu : float
           Substitution rate
        alphabet : str
           Alphabet name (should be standard: 'nuc', 'nuc_gap', 'aa', 'aa_gap')
        pi_dirichlet_alpha : float, optional
            parameter of dirichlet distribution
        W_dirichlet_alpha : float, optional
            parameter of dirichlet distribution
        mu_gamma_alpha : float, optional
            parameter of dirichlet distribution

        Returns
        -------
        GTR_site_specific
            model with randomly sampled frequencies
        """

        from scipy.stats import gamma
        alphabet=alphabets[alphabet]
        gtr = cls(alphabet=alphabet, seq_len=L)
        n = gtr.alphabet.shape[0]

        if pi_dirichlet_alpha:
            pi = 1.0*gamma.rvs(pi_dirichlet_alpha, size=(n,L))
        else:
            pi = np.ones((n,L))

        pi /= pi.sum(axis=0)
        if W_dirichlet_alpha:
            tmp = 1.0*gamma.rvs(W_dirichlet_alpha, size=(n,n))
        else:
            tmp = np.ones((n,n))
        tmp = np.tril(tmp,k=-1)
        W = tmp + tmp.T

        if mu_gamma_alpha:
            mu = gamma.rvs(mu_gamma_alpha, size=(L,))
        else:
            mu = np.ones(L)

        gtr.assign_rates(mu=mu, pi=pi, W=W)
        gtr.mu *= avg_mu/np.mean(gtr.mu)

        return gtr


    @classmethod
    def custom(cls, mu=1.0, pi=None, W=None, **kwargs):
        """
        Create a GTR model by specifying the matrix explicitly

        Parameters
        ----------

         mu : float
            Substitution rate

         W : nxn matrix
            Substitution matrix

         pi : n vector
            Equilibrium frequencies

         **kwargs:
            Key word arguments to be passed

        Keyword Args
        ------------

         alphabet : str
            Specify alphabet when applicable. If the alphabet specification is
            required, but no alphabet is specified, the nucleotide alphabet will be used as
            default.

        """
        gtr = cls(**kwargs)
        gtr.assign_rates(mu=mu, pi=pi, W=W)
        return gtr


    @classmethod
    def infer(cls, sub_ija, T_ia, root_state, pc=1.0,
              gap_limit=0.01, Nit=30, dp=1e-5, **kwargs):
        """
        Infer a GTR model by specifying the number of transitions and time spent in each
        character. The basic equation that is being solved is

        :math:`n_{ij} = pi_i W_{ij} T_j`

        where :math:`n_{ij}` are the transitions, :math:`pi_i` are the equilibrium
        state frequencies, :math:`W_{ij}` is the "substitution attempt matrix",
        while :math:`T_i` is the time on the tree spent in character state
        :math:`i`. To regularize the process, we add pseudocounts and also need
        to account for the fact that the root of the tree is in a particular
        state. the modified equation is

        :math:`n_{ij} + pc = pi_i W_{ij} (T_j+pc+root\_state)`

        Parameters
        ----------

         nija : nxn matrix
            The number of times a change in character state is observed
            between state j and i at position a

         Tia :n vector
            The time spent in each character state at position a

         root_state : np.array
            probability that site a is in state i.

         pc : float
            Pseudocounts, this determines the lower cutoff on the rate when
            no substitutions are observed

         **kwargs:
            Key word arguments to be passed

        Keyword Args
        ------------

         alphabet : str
            Specify alphabet when applicable. If the alphabet specification
            is required, but no alphabet is specified, the nucleotide alphabet will be used as default.

        """
        from scipy import linalg as LA
        gtr = cls(**kwargs)
        gtr.logger("GTR: model inference ",1)
        q = len(gtr.alphabet)
        L = sub_ija.shape[-1]

        n_iter = 0
        n_ija = np.copy(sub_ija)
        n_ija[range(q),range(q),:] = 0
        n_ij = n_ija.sum(axis=-1)

        m_ia = np.sum(n_ija,axis=1) + root_state + pc
        n_a = n_ija.sum(axis=1).sum(axis=0) + pc

        Lambda = np.sum(root_state,axis=0) + q*pc
        p_ia_old=np.zeros((q,L))
        p_ia = np.ones((q,L))/q
        mu_a = np.ones(L)

        W_ij = np.ones((q,q)) - np.eye(q)

        while (LA.norm(p_ia_old-p_ia)>dp) and n_iter<Nit:
            gtr.logger(' '.join(map(str, ['GTR inference iteration',n_iter,'change:',LA.norm(p_ia_old-p_ia)])), 3)
            n_iter += 1
            p_ia_old = np.copy(p_ia)
            S_ij = np.einsum('a,ia,ja',mu_a, p_ia, T_ia)
            W_ij = (n_ij + n_ij.T + pc)/(S_ij + S_ij.T + pc)
            
            avg_pi = p_ia.mean(axis=-1)
            average_rate = W_ij.dot(avg_pi).dot(avg_pi)
            W_ij = W_ij/average_rate
            mu_a *=average_rate
            
            p_ia = m_ia/(mu_a*np.dot(W_ij,T_ia)+Lambda)
            p_ia = p_ia/p_ia.sum(axis=0)
            
            mu_a = n_a/(pc+np.einsum('ia,ij,ja->a', p_ia, W_ij, T_ia))


        if n_iter >= Nit:
            gtr.logger('WARNING: maximum number of iterations has been reached in GTR inference',3, warn=True)
            if LA.norm(p_ia_old-p_ia) > dp:
                gtr.logger('the iterative scheme has not converged',3,warn=True)
        if gtr.gap_index is not None:
            for p in range(p_ia.shape[-1]):
                if p_ia[gtr.gap_index,p]<gap_limit:
                    gtr.logger('The model allows for gaps which are estimated to occur at a low fraction of %1.3e'%p_ia[gtr.gap_index,p]+
                           '\n\t\tthis can potentially result in artifacts.'+
                           '\n\t\tgap fraction will be set to %1.4f'%gap_limit,2,warn=True)
                p_ia[gtr.gap_index,p] = gap_limit
                p_ia[:,p] /= p_ia[:,p].sum()

        gtr.assign_rates(mu=mu_a, W=W_ij, pi=p_ia)
        return gtr


    def _eig(self):
        eigvals, vec, vec_inv = [], [], []
        for pi in range(self.seq_len):
            if len(self.W.shape)>2:
                W = np.copy(self.W[:,:,pi])
                np.fill_diagonal(W, 0)
            elif pi==0:
                np.fill_diagonal(self.W, 0)
                W=self.W

            ev, evec, evec_inv = self._eig_single_site(W,self.Pi[:,pi])
            eigvals.append(ev)
            vec.append(evec)
            vec_inv.append(evec_inv)

        self.eigenvals = np.array(eigvals).T
        self.v = np.swapaxes(vec,0,-1)
        self.v_inv = np.swapaxes(vec_inv, 0,-1)


    def _make_expQt_interpolator(self):
        """Function that evaluates the exponentiated substitution matrix at multiple
        time points and constructs a linear interpolation object
        """
        self.rate_scale = self.average_rate().mean()
        t_grid = (1.0/self.rate_scale)*np.concatenate((np.linspace(0,.1,11)[:-1],
                                                     np.linspace(.1,1,21)[:-1],
                                                     np.linspace(1,5,21)[:-1],
                                                     np.linspace(5,10,11)))
        stacked_expQT = np.stack([self._expQt(t) for t in t_grid], axis=0)

        from scipy.interpolate import interp1d
        self.expQt_interpolator = interp1d(t_grid, stacked_expQT, axis=0,
                                           assume_sorted=True, copy=False, kind='linear')


    def _expQt(self, t):
        """Raw numerical matrix exponentiation using the diagonalized matrix.
        This is the computational bottleneck in many simulations.

        Parameters
        ----------
        t : float
            time

        Returns
        -------
        np.array
            stack of matrices for each site
        """
        eLambdaT = np.exp(t*self.mu*self.eigenvals)
        return np.einsum('jia,ja,kja->ika', self.v, eLambdaT, self.v_inv)


    def expQt(self, t):
        if t*self.rate_scale<10 and self.approximate:
            return self.expQt_interpolator(t)
        else:
            return self._expQt(t)


    def prop_t_compressed(self, seq_pair, multiplicity, t, return_log=False):
        print("NOT IMPEMENTED")


    def propagate_profile(self, profile, t, return_log=False):
        """
        Compute the probability of the sequence state of the parent
        at time (t+t0, backwards), given the sequence state of the
        child (profile) at time t0.

        Parameters
        ----------

         profile : numpy.array
            Sequence profile. Shape = (L, a),
            where L - sequence length, a - alphabet size.

         t : double
            Time to propagate

         return_log: bool
            If True, return log-probability

        Returns
        -------
`
         res : np.array
            Profile of the sequence after time t in the past.
            Shape = (L, a), where L - sequence length, a - alphabet size.

        """
        Qt = self.expQt(t)
        res = np.einsum('ai,ija->aj', profile, Qt)

        return  np.log(np.maximum(ttconf.TINY_NUMBER,res)) if return_log else np.maximum(0,res)


    def evolve(self, profile, t, return_log=False):
        """
        Compute the probability of the sequence state of the child
        at time t later, given the parent profile.

        Parameters
        ----------

         profile : numpy.array
            Sequence profile. Shape = (L, a),
            where L - sequence length, a - alphabet size.

         t : double
            Time to propagate

         return_log: bool
            If True, return log-probability

        Returns
        -------

         res : np.array
            Profile of the sequence after time t in the future.
            Shape = (L, a), where L - sequence length, a - alphabet size.

        """
        Qt = self.expQt(t)
        res = np.einsum('ai,jia->aj', profile, Qt)

        return np.log(res) if return_log else res


    def prob_t(self, seq_p, seq_ch, t, pattern_multiplicity = None,
               return_log=False, ignore_gaps=True):
        """
        Compute the probability to observe seq_ch (child sequence) after time t starting from seq_p
        (parent sequence).

        Parameters
        ----------

         seq_p : character array
            Parent sequence

         seq_c : character array
            Child sequence

         t : double
            Time (branch len) separating the profiles.

         pattern_multiplicity : numpy array
            If sequences are reduced by combining identical alignment patterns,
            these multplicities need to be accounted for when counting the number
            of mutations across a branch. If None, all pattern are assumed to
            occur exactly once.

         return_log : bool
            It True, return log-probability.

        Returns
        -------
         prob : np.array
            Resulting probability

        """
        if t<0:
            logP = -ttconf.BIG_NUMBER
        else:
            tmp_eQT = self.expQt(t)
            bad_indices=(tmp_eQT==0)
            logQt = np.log(tmp_eQT + ttconf.TINY_NUMBER*(bad_indices))
            logQt[np.isnan(logQt) | np.isinf(logQt) | bad_indices] = -ttconf.BIG_NUMBER
            seq_indices_c = np.zeros(len(seq_ch), dtype=int)
            seq_indices_p = np.zeros(len(seq_p), dtype=int)
            for ai, a in self.alphabet:
                seq_indices_p[seq_p==a] = ai
                seq_indices_c[seq_ch==a] = ai

            if len(logQt.shape)==2:
                logP = np.sum(logQt[seq_indices_p, seq_indices_c]*pattern_multiplicity)
            else:
                logP = np.sum(logQt[seq_indices_p, seq_indices_c, np.arange(len(seq_c))]*pattern_multiplicity)

        return logP if return_log else np.exp(logP)


    def average_rate(self):
        if self.Pi.shape[1]>1:
            return np.einsum('a,ia,ij,ja->a',self.mu, self.Pi, self.W, self.Pi)
        else:
            return self.mu*np.einsum('ia,ij,ja->a',self.Pi, self.W, self.Pi)

