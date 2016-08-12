from distribution import Distribution

class BranchLenInterpolator (Distribution):

    """
    """

    def __init__(self, node):
        """
        @brief      { constructor_description }

        @param      self  The object
        @param      node  The node

        """
        self.node = node
        if node.up is None:
            node.branch_neg_log_prob = None
            return None

        if not hasattr(node, 'gamma'):
            self.gamma = 1.0

        if not hasattr(node, 'merger_rate') or node.merger_rate is None:
            self.merger_rate = ttconf.BRANCH_LEN_PENALTY

        # optimal branch length
        mutation_length = node.mutation_length

        if mutation_length < np.min((1e-5, 0.1*self.one_mutation)): # zero-length
            grid = ttconf.MAX_BRANCH_LENGTH * (np.linspace(0, 1.0 , n)**2)

        else: # branch length is not zero

            sigma = mutation_length #np.max([self.average_branch_len, mutation_length])
            # from zero to optimal branch length
            grid_left = mutation_length * (1 - np.linspace(1, 0.0, n/3)**2.0)
            # from optimal branch length to the right (--> 3*branch lengths),
            grid_right = mutation_length + (3*sigma*(np.linspace(0, 1, n/3)**2))
            # far to the right (3*branch length ---> MAX_LEN), very sparse
            far_grid = grid_right.max() + ttconf.MAX_BRANCH_LENGTH*np.linspace(0, 1, n/3)**2

            grid = np.concatenate((grid_left,grid_right[1:],far_grid[1:])
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


        super(BranchLenInterpolator, self).__init__(grid, log_prob)

    @property
    def gamma(self):
        if hasattr(self.node, 'gamma'):
           return self.node.gamma
        else:
           return 1.0

    @gamma.setter
    def gamma(self, value):
        # rescale the interpolator
        #
        #

    def __call__(self, x):
        res = self.merger_rate*x
        res += super.__call__(x/self.gamma)

    def __mul__(self, other):

        res = BranchLenInterpolator(super(BranchLenInterpolator, self).__mul__(other))







