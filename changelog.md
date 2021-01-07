# 0.8.1 -- bug fixe amino acid profile map.

# 0.8.0 -- drop python 2.7 support, bug fixes.

# 0.7.6 -- catch of distributions are too short for calculating confidence intervals.

# 0.7.5 -- fix desync of peak from grid of distributions after pruning

# 0.7.4 -- bug fix in reconstruct discrete trait routine

The `reconstruct_discrete_traits` wrapper function didn't handle missing data correctly (after the changed released in 0.7.2) which resulted in alphabets and weights of different lengths.


# 0.7.3 -- bug fix in average rate calculation

This release fixes a problem that surfaced when inferring GTR models from trees of very similar sequences but quite a few gaps. This resulted in mutation counts like so:

A: [[ 0.  1.  8.  3.  0.]
C:  [ 1.  0.  2.  7.  0.]
G:  [ 9.  0.  0.  2.  0.]
T:  [ 1. 23.  6.  0.  0.]
-:  [46. 22. 28. 38.  0.]]

As a result, the rate "to gap" is inferred quite high, while the equilibrium gap fraction is low. Since we cap the equilibrium gap fraction from below to avoid reconstruction problems when branches are very short, this resulted in an average rate that had substantial contribution from and assumed 1% equilibrum gap frequency where gaps mutate at 20times the rate as others. Since gaps are ignored in distance calculations anyway, it is more sensible to exclude these transitions from the calculation of the average rate. This is now happening in line 7 of treetime/gtr.py. The average rate is restricted to mutation substitutions from non-gap states to any state.


# 0.7.2 -- weights in discrete trait reconstruction
This release implements a more consistent handling of weights (fixed equilibrium frequencies) in discrete state reconstruction.
It also fixes a number of problems in who the arguments were processed.
TreeTime now allows
 * unobserved discrete states
 * uses expected time-in-tree instead of observed time-in-tree in GTR estimation when weights are fixed. The former resulted in very unstable rate estimates.


# 0.7.0 -- restructuring

## Major changes
This release largely includes changes under the hood, some of which also affect how treetime behaves. The biggest changes are
 * sequence data handling is now done by a separate class `SequenceData`. There is now a clear distinction between input data that is never changed and inferred sequences. This class also provides consolidated set of functions to convert sparse, compressed, and full sequence representations into each other.
 * sequences are now unicode when running from python3. This does not seem to come with a measurable performance hit compared to byte sequences as long as all characters are ASCII. Moving away from bytes to unicode proved much less hassle than converting sequences back and forth from unicode to bytes during IO.
 * Ancestral state reconstruction no longer reconstructs the state of terminal nodes by default and sequence accessors and output will return the input data by default. Reconstruction is optional.
 * The command-line mugration model inference now optimize the overall rate numerically and is hence no longer making a short-branch length assumption.
 * TreeTime raises now a number of custom errors rather than returning success or error codes. This should result in fewer "silent errors" that cause problems downstream.

## Minor new features
In addition, we implemented a number of other changes to the interface
 * `treetime`, `treetime clock` now accept the arguments `--name-column` and `-date-column` to explicitly specify the metadata columns to be used as name or date
 * `treetime mugration` accepts a `--name-column` argument.

## Bug fixes
 * scaling of skyline confidence intervals was wrong. It now reflects the inverse second derivative in log-space
 * catch problems after rerooting associated with missing attributes in the newly generated root node.
 * make conversion from calendar dates to numeric dates and vice versa compatible and remove approximate handling of leap-years.
 * avoid overwriting content of output directory with default names
 * don't export inferred dates of tips labeled as `bad_branch`.