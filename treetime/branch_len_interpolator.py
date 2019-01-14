from __future__ import print_function, division, absolute_import
import numpy as np
from treetime import config as ttconf
from .distribution import Distribution

class BranchLenInterpolator (Distribution):
    """
    Tjis class defines the methods to manipulate the branch length probability
    distributions.

    """

    def __init__(self, node, gtr, one_mutation=None, min_width=ttconf.MIN_INTEGRATION_PEAK,
                 branch_length_mode = 'joint', pattern_multiplicity = None,
                 n_grid_points = ttconf.BRANCH_GRID_SIZE, ignore_gaps=True):

        self.node = node
        self.gtr = gtr
        if node.up is None:
            raise Exception("Cannot create branch length interpolator for the root node.")

        self._gamma = 1.0

        self._merger_cost = None
        if one_mutation is None:
            L = node.sequence.shape[0]
            one_mutation = 1.0/L

        # optimal branch length
        mutation_length = node.mutation_length
        if mutation_length < np.min((1e-5, 0.1*one_mutation)): # zero-length
            short_range = 10*one_mutation
            grid = np.concatenate([short_range*(np.linspace(0, 1.0 , n_grid_points/2)[:-1]),
                (short_range + (ttconf.MAX_BRANCH_LENGTH - short_range)*(np.linspace(0, 1.0 , n_grid_points/2+1)**2))])

        else: # branch length is not zero
            sigma = mutation_length #np.max([self.average_branch_len, mutation_length])
            # from zero to optimal branch length
            grid_left = mutation_length * (1 - np.linspace(1, 0.0, n_grid_points/3)**2.0)
            grid_zero = grid_left[1]*np.logspace(-20,0,6)[:5]
            grid_zero2 = grid_left[1]*np.linspace(0,1,10)[1:-1]
            # from optimal branch length to the right (--> 3*branch lengths),
            grid_right = mutation_length + (3*sigma*(np.linspace(0, 1, n_grid_points/3)**2))
            # far to the right (3*branch length ---> MAX_LEN), very sparse
            far_grid = grid_right.max() + ttconf.MAX_BRANCH_LENGTH*np.linspace(0, 1, n_grid_points/3)**2

            grid = np.concatenate((grid_zero,grid_zero2, grid_left,grid_right[1:],far_grid[1:]))
            grid.sort() # just for safety

        if branch_length_mode=='input':
            # APPROXIMATE HANDLING OF BRANCH LENGTH PROPAGATOR WHEN USING INPUT BRANCH LENGTH
            # branch length are estimated from as those maximizing the likelihood and the
            # sensitivity of the likelihood depends on the branch length (gets soft for long branches)
            # observed differences scale as p = p_0 (1-exp(-l/p_0)) where p_0 is the distance of random sequence
            # (3/4 for nucleotides, more like 0.9 for amino acids). The number of observable
            # substitutions fluctuates by dp = \sqrt{p(1-p)/L} which corresponds to fluctuation
            # in branch length of dp = dl exp(-l/p0). A Gaussian approximation for the branch length would
            # therefore have variance p(1-p)e^{2l/p0}/L. Substituting p results in
            # p_0(1-exp(-l/p0))(1-p_0(1-exp(-l/p0)))e^{2l/p0}/L which can be slightly rearranged to
            # p_0(exp(l/p0)-1)(exp(l/p0)-p_0(exp(l/p0)-1))/L

            p0 = 1.0-np.sum(self.gtr.Pi**2)
            # variance_scale = one_mutation*ttconf.OVER_DISPERSION
            if mutation_length<0.05:
                # for short branches, the number of mutations is poissonian. the prob of a branch to have l=mutation_length*L
                # mutations when its length is k, is therefor e^{-kL}(kL)^(Ll)/(Ll)!. Ignoring constants, the log is
                # -kL + lL\log(k)
                log_prob = np.array([ k - mutation_length*np.log(k+ttconf.MIN_BRANCH_LENGTH*one_mutation) for k in grid])/one_mutation
                log_prob -= log_prob.min()
            else:
                # make it a Gaussian
                #sigma_sq = (mutation_length+one_mutation)*variance_scale
                l = (mutation_length+one_mutation)
                nm_inv = np.exp(l/p0)
                sigma_sq = p0*(nm_inv-1)*(nm_inv - p0*(nm_inv-1))*one_mutation
                sigma = np.sqrt(sigma_sq+ttconf.MIN_BRANCH_LENGTH*one_mutation)
                log_prob = np.array(np.min([[ 0.5*(mutation_length-k)**2/sigma_sq for k in grid],
                                             100 + np.abs([(mutation_length-k)/sigma for k in grid])], axis=0))
        elif branch_length_mode=='marginal':
            if hasattr(node, 'profile_pair'):
                log_prob = np.array([-self.gtr.prob_t_profiles(node.profile_pair,
                                                        pattern_multiplicity,
                                                        k,
                                                        return_log=True)
                                    for k in grid])
            else:
                raise Exception("profile pairs need to be assigned to node")


        elif branch_length_mode=='joint':
            if not hasattr(node, 'compressed_sequence'):
                #FIXME: this assumes node.sequence is set, but this might not be the case if
                # ancestral reconstruction is run with final=False
                if hasattr(node, 'sequence'):
                    seq_pairs, multiplicity = self.gtr.compress_sequence_pair(node.up.sequence,
                                                                          node.sequence,
                                                                          ignore_gaps=ignore_gaps)
                    node.compressed_sequence = {'pair':seq_pairs, 'multiplicity':multiplicity}
                else:
                    raise Exception("uncompressed sequence needs to be assigned to nodes")

            log_prob = np.array([-self.gtr.prob_t_compressed(node.compressed_sequence['pair'],
                                                    node.compressed_sequence['multiplicity'],
                                                    k,
                                                    return_log=True)
                                for k in grid])
        else:
            raise Exception("unknown branch length mode! "+branch_length_mode)
        # tmp_dis = Distribution(grid, log_prob, is_log=True, kind='linear')
        # norm = tmp_dis.integrate(a=tmp_dis.xmin, b=tmp_dis.xmax, n=200)
        super(BranchLenInterpolator, self).__init__(grid, log_prob, is_log=True,
                                                    kind='linear', min_width=min_width)


    @property
    def gamma(self):
        return self._gamma

    @gamma.setter
    def gamma(self, value):
        self._gamma = max(ttconf.TINY_NUMBER, value)

    @property
    def merger_cost(self):
        return self._merger_cost

    @merger_cost.setter
    def merger_cost(self, cost_func):
        self._merger_cost = cost_func
        self._peak_idx = np.argmin(self.__call__(self.x))
        self._peak_pos = self.x[self._peak_idx]
        if self.kind=='linear': # can't mess like this with non-linear interpolation
            deltay = self.__call__(self.peak_pos) - self._peak_val
            self._peak_val += deltay
            self._func.y -= deltay

    @property
    def peak_pos(self):
        return super(BranchLenInterpolator,self).peak_pos/self.gamma

    @property
    def support(self):
        return self._support/self.gamma

    @property
    def fwhm(self):
        return super(BranchLenInterpolator,self).fwhm/self.gamma

    def __call__(self, x, tnode=None, multiplicity=None):
        res = super(BranchLenInterpolator, self).__call__(x*self.gamma)
        if self.merger_cost is not None:
            if tnode is None:
                tnode = self.node.time_before_present
            if multiplicity is None:
                multiplicity = len(self.node.up.clades)
            res += self.merger_cost(tnode, x, multiplicity=multiplicity)
        return res

    def __mul__(self, other):
        res = BranchLenInterpolator(super(BranchLenInterpolator, self).__mul__(other))
        return res


