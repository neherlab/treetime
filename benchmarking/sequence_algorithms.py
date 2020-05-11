from __future__ import print_function, division
import numpy as np
from Bio import Phylo


if __name__ == '__main__':
    from treetime.seq_utils import normalize_profile, prof2seq, seq2prof
    from treetime.gtr import GTR    

    gtr = GTR.standard('JC69')
    dummy_prof = np.random.random(size=(10000,5))

    # used a lot (300us)
    norm_prof = normalize_profile(dummy_prof)[0]

    # used less but still a lot (50us)
    gtr.evolve(norm_prof, 0.1)

    # used less but still a lot (50us)
    gtr.propagate_profile(norm_prof, 0.1)

    # used only in final, sample_from_prof=False speeds it up (600us or 300us)
    seq, p, seq_ii = prof2seq(norm_prof, gtr, sample_from_prof=True, normalize=False)

    # used only initially (slow, 5ms)
    tmp_prof = seq2prof(seq, gtr.profile_map)
