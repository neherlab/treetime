# 0.9.3

 * Add extra error class for "unknown" (==unhandled) errors
 * Wrap `run` function and have it optionally raise unhandled exceptions as `TreeTimeUnknownError`.
   This is mainly done to improve interaction with `augur` that uses `TreeTime` internals as a library.
   (both by @anna-parker with input from @victorlin)

[PR #206](https://github.com/neherlab/treetime/pull/206)
[PR #208](https://github.com/neherlab/treetime/pull/208)

# 0.9.2
bug fix release:
 * CLI now works for windows (thanks @corneliusroemer for the fix)
 * fixes vcf parsing. haploid no-calls were not properly parsed and treated as reference (thanks @jodyphelan for the issue).
 * fix file names in CLI output. (thanks @gtonkinhill)

# 0.9.1
This release is mostly a bug-fix release and contains some additional safeguards against unwanted side-effects of greedy polytomy resolution.

 * resolve polytomies only when significant LH gain can be achieved
 * performance enhancement in pre-order iteration during marginal time tree estimate when hitting large polytomies.
 * allow users to set branch specific rates (only when used as a library)

# 0.9.0

This release contains several major changes to how TreeTime calculates time scaled phylogenies.
Most of this is work by @anna-parker!

 * implements convolutions needed for marginal time tree inference using FFT.
   Previously, these were calculated by explicit integration using optimized irregular grids.
   Using FFT requires regular (and hence much finer/larger) grids, but greatly reduces computational complexity from `n^2` to `n log(n)`, where `n` is the number of grid points.
   The FFT feature can be switched on an off with the `use_fft` attribute of the ClockTree class.
 * Using FFT in convolutions required moving the contributions of the coalescent models from th branches to the nodes.
   This should not change the results in any way, but cleans up the code.
 * The number concurrent of lineages determines the rate of coalescence.
   This can now optionally be calculated using the uncertainty of the timing of merger events, instead of the step functions used previously.
 * Adds a subcommand to read in ancestral reassortment graphs of two segments produced by [TreeKnit](https://github.com/PierreBarrat/TreeKnit.jl). This command takes two trees and a file with MCCs inferred by TreeKnit. See [these docs](https://treetime.readthedocs.io/en/latest/commands.html#arg) for command line usage.

# 0.8.6
 * optionally allow incomplete alignment [PR #178](https://github.com/neherlab/treetime/pull/178)
 * reduce memory footprint through better clean up and optimizing types. [PR #179](https://github.com/neherlab/treetime/pull/179)

# 0.8.5
 * bug fixes related to edge cases were sequences consist only of missing data
 * bug fix when the CLI command `treetime` is run without alignment
 * more robust behavior when parsing biopython alignments (id vs name of sequence records)
 * drop python 3.5 support

# 0.8.4 -- re-release of 0.8.3.1

# 0.8.3.1 -- bug fix related to Bio.Seq.Seq now bytearray

 * Biopython changed the representation of sequences from strings to bytearrays. This caused crashes of mugration inference with more than 62 states as states than exceeded the ascii range. This fix now bypasses Bio.Seq in the mugration analysis.

# 0.8.3 -- unpin biopython version

 * Biopython 1.77 and 1.78 had a bug in their nexus export. This is fixed in 1.79. We now explictly exclude the buggy versions but allow others.

# 0.8.2 -- bug fixes and small feature additions
This release fixes a few bugs and adds a few features

 * output statistics of different iterations of the treetime optimization loop (trace-log, thanks to @ktmeaton)
 * speed ups by @akislyuk
 * fix errors with dates in the distant future
 * better precision of tabular skyline output
 * adds clock-deviation to the root-to-tip output of the `clock` command


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