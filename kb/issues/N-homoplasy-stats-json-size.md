# Homoplasy: the statistics JSON is large for big alignments

> [!IMPORTANT]
> **Decision required.** The statistics file of `treetime homoplasy` (`--output-homoplasy-stats`, `homoplasy.stats.json`) reaches hundreds of megabytes on large datasets, almost all of it the ranked list of changes involving ambiguous characters. Options and evidence are below.

## Symptom

`treetime homoplasy --output-all` with default settings writes a statistics file of:

- 0.4 MB for `data/flu/h3n2/500`
- 123 MB for `data/sc2/4500`
- 221 MB for `data/mpox/clade-ii/2000`

The ranked list of ambiguous changes (`ambiguous.all.ranked`, each entry with the name of every branch it occurs on) is 72 MB with 31,571 entries for `sc2/4500` and 97 MB with 189,389 entries for `mpox/clade-ii/2000`. The substitution lists are 0.1 to 1.3 MB, `taxa` is at most 0.65 MB, and `ambiguous.sites` has up to 188,636 rows (mpox).

Most of the ambiguous changes are the sequence ends that gap filling turns into `N` ([N-homoplasy-filled-overhangs-dominate-ambiguous-changes.md](N-homoplasy-filled-overhangs-dominate-ambiguous-changes.md)).

The web and desktop apps read the file for the results page (`fn homoplasy_statistics()` in `packages/app-commands/src/results/homoplasy.rs`). They skip `ambiguous.all.ranked` while deserializing, which avoids allocating its entries but not reading and scanning its bytes. Without whitespace, the list is 72.8 MB of the 77.3 MB of the `sc2/4500` file and 100.3 MB of the 107.9 MB of the `mpox/clade-ii/2000` file. With the stream reader of `json_read_file()` in the dev profile, the `sc2/4500` file takes 9.8 s to read as written, and 0.42 s without whitespace and without the list ([M-io-json-read-from-reader-slow.md](M-io-json-read-from-reader-slow.md)). In the dev profile, the mpox file takes 0.37 s to read and 2.55 s to deserialize without the ranked ambiguous list (2.62 s with it); the sc2 file takes 0.27 s, 0.59 s, and 0.86 s.

## Options

- O1. Keep the file as it is. It is complete, and the apps skip the large list
- O2. Drop the branch names from the entries of `ambiguous.all.ranked`, keeping the mutation and its multiplicity
- O3. List only ambiguous changes that occur on two or more branches, as the substitution report focuses on recurrent changes
