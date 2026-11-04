# Root dates move by years between v1 commits without a traced cause

## Problem

Between v1 at `b92574e3e` and v1 at `0e522d97e`, the root dates of five Nextstrain time trees moved by more than a year, on the same commands and inputs. The commits between them change convolutions at hard ends (exact tails, sections along long gaps), the FFT lengths (powers of two), the sparse reconstruction, and the clock filter. On the builds that these changes were measured on, they moved node dates by at most 0.45 years. A change whose intended effect is numerical, and that moves a root by a decade, shows that one of the two dates carries a large numerical error, or that the root posterior is so flat or multimodal that small changes move its summary. Either case limits the comparison of deep-node dates with v0 ([H-timetree-deep-node-dates-differ-from-v0-by-decades.md](H-timetree-deep-node-dates-differ-from-v0-by-decades.md)).

## Evidence

Root dates of the comparison of TreeTime v0 and v1 in the Nextstrain workflows (pipelandia). Round 1 ran an earlier v1 commit, round 2 `b92574e3e`, round 3 `0e522d97e`; v0 is the same frozen result in all rounds.

| Run | Build | v1, round 1 | v1, round 2 | v1, round 3 | v0 | Cause |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| norovirus, published | GII.17/genome | 1913.3 | 1851.8 | 1900.1 | 1908.3 | traced: a v1 binary of `3e9eb55be` gives 1851.84, a binary of `8364583af` (exact convolution samples at hard ends) gives 1900.09, on the same root split of 294 and 22 tips; the defect of empty time distributions moved the root in round 2 without stopping the run |
| rubella, published | E1 | | 1896.9 | 1902.9 | 1919.9 | partly traced: the clades of 33, 1,316, and 1,674 tips meet near the root; round 2 pairs the 33-tip clade with the 1,316-tip clade, as v0 does, and round 3 pairs it with the 1,674-tip clade. The cause of the changed pairing is not traced |
| yellow-fever, example | genome | 1506.4 | 1507.9 | 1519.3 | 1560.2 | not traced |
| dengue, published | all/genome | 1317.8 | 1289.2 | 1296.7 | 1315.4 | not traced |
| measles, published | genome/global | 1891.6 | 1875.5 | 1873.8 | 1787.4 | not traced (estimated rate) |

The median internal-node date difference against v0 stays at 0.13 to 0.14 years over the three rounds, so the movement concentrates at the deep nodes. The polytomy resolution can change the topology near the root after a round-off change of node times ([N-timetree-polytomy-resolution-flips-on-round-off-of-node-times.md](N-timetree-polytomy-resolution-flips-on-round-off-of-node-times.md)), and the slope of the straight lines at soft ends carries FFT round-off ([M-distribution-convolution-soft-end-line-slope-carries-fft-round-off.md](M-distribution-convolution-soft-end-line-slope-carries-fft-round-off.md)); neither is shown to cause these moves.

> [!IMPORTANT]
> **Investigation required.** The cause of the moves of yellow-fever genome (11.4 years), dengue all/genome (7.5 years), measles genome/global (1.7 years), and of the changed root pairing of rubella E1 is not located. Evidence to collect, with the method of the GII.17 attribution: v1 binaries of the single commits between `b92574e3e` and `0e522d97e`, each run on the `refine` step of the build with the command line of the replay; for the commit that moves the root, the root posterior of both binaries (one peak or several, width), the root split, and the resolved polytomies near the root.

## Reproduction

Each build runs from the data package of the pipelandia comparison with `just treetime1_external=path:<binary> check <pipeline> --profile production` in the pipelandia repository; the command line of each `refine` step is in its check log.
