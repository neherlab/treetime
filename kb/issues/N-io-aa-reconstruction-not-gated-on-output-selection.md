# AA reconstruction runs unconditionally regardless of output selection

When `ReconstructedAaFasta` is not in the resolved output set, the ancestral and timetree commands still run full amino acid reconstruction. The result is computed then discarded.

## Impact

Wasted compute on large alignments with many coding sequences. No correctness issue.

## Fix

Gate the AA reconstruction call on `resolved.contains(&OutputSelection::ReconstructedAaFasta)`. Skip entirely when the output is not selected.
