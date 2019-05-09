VERBOSE = 3

BIG_NUMBER = 1e10
TINY_NUMBER = 1e-12
SUPERTINY_NUMBER = 1e-24
MIN_LOG = -1e8 # minimal log value
MIN_BRANCH_LENGTH = 1e-3 # fraction of length 'one_mutation' that is used as lower cut-off for branch lengths in GTR
OVER_DISPERSION = 10

# distribution parameters
BRANCH_GRID_SIZE_ROUGH = 200
NODE_GRID_SIZE_ROUGH = 60
N_INTEGRAL_ROUGH = 60

BRANCH_GRID_SIZE = 250
NODE_GRID_SIZE = 100
N_INTEGRAL = 100

BRANCH_GRID_SIZE_FINE = 300
NODE_GRID_SIZE_FINE = 180
N_INTEGRAL_FINE = 150

BRANCH_GRID_SIZE_ULTRA = 400
NODE_GRID_SIZE_ULTRA = 400
N_INTEGRAL_ULTRA = 250

MIN_INTEGRATION_PEAK = 0.001

# clocktree parameters
BRANCH_LEN_PENALTY = 0
MAX_BRANCH_LENGTH = 4.0          # only relevant for branch length optimization and time trees - upper boundary of interpolator objects
NINTEGRAL = 300
REL_TOL_PRUNE = 0.01
REL_TOL_REFINE = 0.05
NIQD = 3

#
SUCCESS = "success"
ERROR = "error"

# treetime
# autocorrelated molecular clock coefficients
MU_ALPHA = 1
MU_BETA = 1

