# Configuration file format for multi-partition analysis

## Motivation

The `optimize.md` design document notes: "For more complex inputs (multiple alignments and models along with discrete characters we might need a config file format)." Currently each command accepts a single alignment via `--aln`. Multi-segment genomes, mixed data types, and per-partition model specification all require a configuration system that the CLI flag approach cannot scale to.

## Use cases

### Multi-segment genomes

Influenza has 8 genome segments, each with independent evolution but shared tree topology. Currently requires concatenation, losing per-segment model information. The v1 partition system already supports multiple independent partitions on the same tree -- the gap is input specification.

### Mixed data types

Phylogeographic analysis combines nucleotide sequences (4-state alphabet, GTR model) with discrete location traits (variable-size alphabet, symmetric or asymmetric rate matrix). The `mugration` command handles location separately from `ancestral`/`timetree`. A unified config would enable joint inference.

### Per-partition model specification

Different genome regions evolve under different substitution models. Protein-coding regions may use codon models or amino acid models, while non-coding regions use nucleotide models. Each partition needs its own model assignment.

## Prior art

### NEXUS SETS block

The NEXUS file format (Maddison, Swofford, and Maddison 1997, doi:10.1093/sysbio/46.4.590) defines `CHARSET` and `CHARPARTITION` commands for specifying data subsets:

```
#nexus
begin sets;
  charset gene1 = 1-500;
  charset gene2 = 501-1200;
  charpartition mine = GTR+G:gene1, HKY+I:gene2;
end;
```

MrBayes (Ronquist and Huelsenbeck 2003, doi:10.1093/bioinformatics/btg180) popularized NEXUS charsets for mixed-model Bayesian partitioned analysis. PartitionFinder (Lanfear et al. 2012, doi:10.1093/molbev/mss020) automates partition scheme selection.

### IQ-TREE partition modes

IQ-TREE 2 (Minh et al. 2020, doi:10.1093/molbev/msaa015) supports three partition modes: edge-linked proportional (`-spp`), edge-linked equal (`-q`), and edge-unlinked (`-sp`). Terrace-aware search (Chernomor, von Haeseler, and Minh 2016, doi:10.1093/sysbio/syw037) handles partitioned data with missing taxa.

### RAxML-NG partition format

```
GTR+G, gene1 = 1-500
HKY+I, gene2 = 501-1200
```

Column-range-based partitioning of a concatenated alignment. Kozlov et al. 2019, doi:10.1093/bioinformatics/btz305.

### BEAST XML

Full model hierarchy in XML: sequence data, substitution models, clock models, tree priors, operators. Complete but verbose. Drummond and Rambaut 2007, doi:10.1186/1471-2148-7-214.

## Proposed design

YAML or JSON configuration file specifying:

- Partition list: each with alignment path, alphabet, substitution model
- Shared parameters: tree path, clock model, coalescent model
- Output specification

This leverages the existing v1 partition architecture (`PartitionMarginalDense`, `PartitionMarginalSparse`) which already supports multiple independent partitions on the same tree.

## Related known issues

- [N-optimize-multi-alignment-input](../port-known-issues/N-optimize-multi-alignment-input.md) -- optimize accepts only a single alignment
- [N-io-multi-segment-genome-input](../port-known-issues/N-io-multi-segment-genome-input.md) -- multi-segment genome input not wired

## References

1. <a id="ref-1"></a> Maddison, David R., David L. Swofford, and Wayne P. Maddison. 1997. "NEXUS: An Extensible File Format for Systematic Information." _Systematic Biology_ 46(4):590-621. https://doi.org/10.1093/sysbio/46.4.590
2. <a id="ref-2"></a> Ronquist, Fredrik, and John P. Huelsenbeck. 2003. "MrBayes 3: Bayesian Phylogenetic Inference under Mixed Models." _Bioinformatics_ 19(12):1572-1574. https://doi.org/10.1093/bioinformatics/btg180
3. <a id="ref-3"></a> Lanfear, Robert, Brett Calcott, Simon Y. W. Ho, and Stephane Guindon. 2012. "PartitionFinder: Combined Selection of Partitioning Schemes and Substitution Models for Phylogenetic Analyses." _Molecular Biology and Evolution_ 29(6):1695-1701. https://doi.org/10.1093/molbev/mss020
4. <a id="ref-4"></a> Chernomor, Olga, Arndt von Haeseler, and Bui Quang Minh. 2016. "Terrace Aware Data Structure for Phylogenomic Inference from Supermatrices." _Systematic Biology_ 65(6):997-1008. https://doi.org/10.1093/sysbio/syw037
5. <a id="ref-5"></a> Minh, Bui Quang, Heiko A. Schmidt, Olga Chernomor, Dominik Schrempf, Michael D. Woodhams, Arndt von Haeseler, and Robert Lanfear. 2020. "IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era." _Molecular Biology and Evolution_ 37(5):1530-1534. https://doi.org/10.1093/molbev/msaa015
6. <a id="ref-6"></a> Kozlov, Alexey M., Diego Darriba, Tomas Flouri, Benoit Morel, and Alexandros Stamatakis. 2019. "RAxML-NG: A Fast, Scalable and User-Friendly Tool for Maximum Likelihood Phylogenetic Inference." _Bioinformatics_ 35(21):4453-4455. https://doi.org/10.1093/bioinformatics/btz305
7. <a id="ref-7"></a> Drummond, Alexei J., and Andrew Rambaut. 2007. "BEAST: Bayesian Evolutionary Analysis by Sampling Trees." _BMC Evolutionary Biology_ 7:214. https://doi.org/10.1186/1471-2148-7-214
