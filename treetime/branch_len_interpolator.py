import numpy as np
import config as ttconf
from distribution import Distribution

class BranchLenInterpolator (Distribution):

    """
    """

    def __init__(self, node, gtr, one_mutation = None):
        """
        @brief      { constructor_description }

        @param      self  The object
        @param      node  The node

        """

        self.node = node
        self.gtr = gtr
        if node.up is None:
            raise Exception("Cannot create branch length interpolator for the root node.")

        self._gamma = 1.0

        self._merger_rate = ttconf.BRANCH_LEN_PENALTY
        if one_mutation is None:
            one_mutation = 1.0/node.sequence.shape[0]
        # optimal branch length
        mutation_length = node.mutation_length
        n_grid_points = ttconf.BRANCH_GRID_SIZE
        if mutation_length < np.min((1e-5, 0.1*one_mutation)): # zero-length
            grid = ttconf.MAX_BRANCH_LENGTH * (np.linspace(0, 1.0 , n_grid_points)**2)

        else: # branch length is not zero
            sigma = mutation_length #np.max([self.average_branch_len, mutation_length])
            # from zero to optimal branch length
            grid_left = mutation_length * (1 - np.linspace(1, 0.0, n_grid_points/3)**2.0)
            # from optimal branch length to the right (--> 3*branch lengths),
            grid_right = mutation_length + (3*sigma*(np.linspace(0, 1, n_grid_points/3)**2))
            # far to the right (3*branch length ---> MAX_LEN), very sparse
            far_grid = grid_right.max() + ttconf.MAX_BRANCH_LENGTH*np.linspace(0, 1, n_grid_points/3)**2

            grid = np.concatenate((grid_left,grid_right[1:],far_grid[1:]))
            grid.sort() # just for safety

        if not hasattr(node, 'compressed_sequence'):
            seq_pairs, multiplicity = self.gtr.compress_sequence_pair(node.up.sequence,
                                                                      node.sequence,
                                                                      ignore_gaps = self.ignore_gaps)
            node.compressed_sequence = {'pair':seq_pairs, 'multiplicity':multiplicity}

        log_prob = np.array([-self.gtr.prob_t_compressed(node.compressed_sequence['pair'],
                                                node.compressed_sequence['multiplicity'],
                                                k,
                                                return_log=True)
                    for k in grid])


        super(BranchLenInterpolator, self).__init__(grid, log_prob, is_log=True)



    @property
    def gamma(self):
       return self._gamma

    @gamma.setter
    def gamma(self, value):
        self._gamma = value

    @property
    def merger_rate(self):
       return self._merger_rate

    @merger_rate.setter
    def merger_rate(self, value):
        self._gamma = value

    def __call__(self, x):
        res = self.merger_rate*x
        res += super(BranchLenInterpolator, self).__call__(x/self.gamma)
        return res

    def __mul__(self, other):
        res = BranchLenInterpolator(super(BranchLenInterpolator, self).__mul__(other))


