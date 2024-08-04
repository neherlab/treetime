import numpy as np

alphabet='ACGT'
default_profile = np.ones(len(alphabet), dtype=float)
#                        "A,C,G,T"
profile_map = {'A':np.array([1,0,0,0], dtype=float),
            'C':np.array([0,1,0,0], dtype=float),
            'G':np.array([0,0,1,0], dtype=float),
            'T':np.array([0,0,0,1], dtype=float),
            'R':np.array([1,0,1,0], dtype=float),
            'Y':np.array([0,1,0,1], dtype=float),
            'S':np.array([0,1,1,0], dtype=float),
            'W':np.array([1,0,0,1], dtype=float),
            'K':np.array([0,0,1,1], dtype=float),
            'M':np.array([1,1,0,0], dtype=float),
            'B':np.array([0,1,1,1], dtype=float),
            'D':np.array([1,0,1,1], dtype=float),
            'H':np.array([1,1,0,1], dtype=float),
            'V':np.array([1,1,1,0], dtype=float),
            'N':default_profile}
