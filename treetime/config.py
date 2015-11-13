import numpy as np

BRANCH_LEN_PENALTY = 0
MAX_BRANCH_LENGTH = 1
BIG_NUMBER = 1e10
TINY_NUMBER = 1e-50

MIN_T = -1e5
MAX_T =  1e5

WIDTH_DELTA = 1e-10 # width of the delta function
MIN_LOG = -1e8 # minimal log value

BRANCH_GRID_SIZE = 51
NODE_GRID_SIZE = 300
NODE_GRID_VAR = 0.5 # branch grid covers up to this ratio of the tree depth
