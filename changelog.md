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