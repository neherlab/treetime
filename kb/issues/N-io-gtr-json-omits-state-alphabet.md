# The substitution model JSON of sequence models does not name its states

## Summary

The GTR output (`*.gtr.json`) of the sequence commands writes the equilibrium frequencies `pi` and the rate matrix `W` as bare arrays with `n_states`, but not the alphabet that orders them. The mugration model fills the optional `states` field; sequence models leave it out. A reader cannot tell which frequency belongs to which nucleotide or amino acid without knowing the internal state order of the command that wrote the file.

## Evidence

- `struct GtrOutput` serializes `states` only when it is set [`packages/treetime/src/gtr/get_gtr.rs#L30-L45`](../../packages/treetime/src/gtr/get_gtr.rs#L30-L45)
- `timetree`, `ancestral` and `optimize` on `data/zika/86` write `model_type`, `model_name`, `mu`, `pi`, `W` and `n_states` only

## Impact

- The app shows only the model name and the overall rate `mu` of an optimize run, because it cannot label the frequencies
- External tools that load the file must assume a state order
