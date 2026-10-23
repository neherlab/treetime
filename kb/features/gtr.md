# GTR Substitution Models

## Implemented Models

- [x] Named nucleotide models: JC69, K80, F81, HKY85, T92, TN93
- [x] Amino acid model: JTT92

## Model Inference

- [x] Model inference from data (sparse and dense paths)
- [x] Eigendecomposition (symmetric, caching)
- [x] expQt (matrix exponential for transition probabilities)
- [x] Equilibrium frequencies (pi vector, normalized)
- [x] Exchangeability matrix (W, symmetric, normalized by avg_transition)
- [x] Rate scaling (mu parameter)
- [x] Model document (`--output-model`): one JSON document `ModelOutput` with `alignments` (each file with its number, path, number of sites and indel rate; absent without alignments) and `models`, one entry per model, tagged `kind`. A `sequence` entry holds the partition number, the partition name, the alignment number, the number of sites, in timetree the clock rate of the partition (clock rate of the run times μ), the state characters in the order of π and W, and the model fields W, π, μ and number of states; a `trait` entry (mugration) holds the attribute, the ordered state labels and the model fields. A run without any model reports the output as unavailable
- [x] One model per partition: without entries `--model` and `--model-params` apply to every file; with `--partition` entries each entry carries its own `model` and `model_params`, and `--model` or `--model-params` next to entries is an argument error
- [x] Model check against the alphabet: `k80`, `hky85`, `t92` and `tn93` need `nuc`, `jtt92` needs `aa-no-stop`, `pi` needs one value per state ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md))
- [x] One relative rate μ per partition in optimize and timetree with several partitions (`--partition-rates estimate|equal`; rate step with W and π fixed, normalization to a site-weighted mean of 1: [kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md))

## Not Implemented

- [ ] Custom GTR from file (`--custom-gtr` / `GTR.from_file()` in v0)
- [ ] Random GTR generation (`GTR.random()` in v0, used for testing)
- [ ] GTR save to file (`save_to_npz()` in v0)
- [ ] Site-specific models (multi-dimensional pi/W in v0)
- [ ] Branch length optimization at GTR level (v0 `optimal_t()` - moved to partition layer in v1)
- [ ] Sequence probability at GTR level (v0 `prob_t()` - moved to partition layer in v1)
- [ ] State pair compression (v0 `state_pair()` - moved to partition layer in v1)

## References

- Jukes & Cantor (1969). "Evolution of Protein Molecules." In Munro (ed.), Mammalian Protein Metabolism, vol. 3, pp. 21-132. Academic Press
- Kimura (1980). "A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences." J Mol Evol, 16(2):111-120. https://doi.org/10.1007/BF01731581
- Felsenstein (1981). "Evolutionary trees from DNA sequences: a maximum likelihood approach." J Mol Evol, 17(6):368-376. https://doi.org/10.1007/BF01734359
- Hasegawa, Kishino & Yano (1985). "Dating of the human-ape splitting by a molecular clock of mitochondrial DNA." J Mol Evol, 22(2):160-174. https://doi.org/10.1007/BF02101694
- Tamura (1992). "Estimation of the number of nucleotide substitutions when there are strong transition-transversion and G+C-content biases." Mol Biol Evol, 9(4):678-687. https://doi.org/10.1093/oxfordjournals.molbev.a040752
- Jones, Taylor & Thornton (1992). "The rapid generation of mutation data matrices from protein sequences." CABIOS, 8(3):275-282. https://doi.org/10.1093/bioinformatics/8.3.275
- Tamura & Nei (1993). "Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees." Mol Biol Evol, 10(3):512-526. https://doi.org/10.1093/oxfordjournals.molbev.a040023
