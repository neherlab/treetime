VERBOSE = 3

BIG_NUMBER = 1e10
TINY_NUMBER = 1e-12
MIN_LOG = -1e8 # minimal log value
MIN_BRANCH_LENGTH = 1e-3 # fraction of length 'one_mutation' that is used as lower cut-off for branch lengths in GTR

# distribution parameters
BRANCH_GRID_SIZE = 250
NODE_GRID_SIZE = 60
MIN_INTEGRATION_PEAK = 0.001

# clocktree parameters
BRANCH_LEN_PENALTY = 0
MAX_BRANCH_LENGTH = 1.5          # only relevent for time trees - upper boundary of interpolator objects
NINTEGRAL = 300
REL_TOL_PRUNE = 0.01
REL_TOL_REFINE = 0.05
NIQD = 3

# treetime
# autocorrelated molecular clock coefficients
MU_ALPHA = 1
MU_BETA = 1

