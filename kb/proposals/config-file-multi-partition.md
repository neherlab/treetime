# Configuration file format for multi-partition analysis

## Motivation

Complex inputs combine several alignments, their substitution models, and discrete characters, and need a config file format. Partitions of sequences are implemented: the `partition` setting of ancestral, optimize and timetree (see Implemented part). Mixed data types and configurations that span several commands still need a configuration system beyond the per-command settings.

## Use cases

### Multi-segment genomes

Influenza has 8 genome segments, each with independent evolution but shared tree topology. One `--alignment` file per segment gives one partition per segment with its own model and relative rate; `partition` entries choose a different model per segment. Segments with different trees remain open ([kb/issues/N-io-multi-segment-genome-input.md](../issues/N-io-multi-segment-genome-input.md)).

### Mixed data types

Phylogeographic analysis combines nucleotide sequences (4-state alphabet, GTR model) with discrete location traits (variable-size alphabet, symmetric or asymmetric rate matrix). The `mugration` command handles location separately from `ancestral`/`timetree`. A unified config would enable joint inference.

### Per-partition model specification

Different genome regions evolve under different substitution models. Protein-coding regions may use codon models or amino acid models, while non-coding regions use nucleotide models. Each partition needs its own model assignment. The `partition` setting assigns a nucleotide model to each set of sites; codon models and amino-acid partitions next to nucleotide partitions remain open.

## Prior art

### NEXUS SETS block

The NEXUS file format (Maddison, Swofford, and Maddison 1997, https://doi.org/10.1093/sysbio/46.4.590) defines `CHARSET` and `CHARPARTITION` commands for specifying data subsets:

```
#nexus
begin sets;
  charset gene1 = 1-500;
  charset gene2 = 501-1200;
  charpartition mine = GTR+G:gene1, HKY+I:gene2;
end;
```

MrBayes (Ronquist and Huelsenbeck 2003, https://doi.org/10.1093/bioinformatics/btg180) popularized NEXUS charsets for mixed-model Bayesian partitioned analysis. PartitionFinder (Lanfear et al. 2012, https://doi.org/10.1093/molbev/mss020) automates partition scheme selection.

### IQ-TREE partition modes

IQ-TREE 2 (Minh et al. 2020, https://doi.org/10.1093/molbev/msaa015) supports three partition modes: edge-linked proportional (`-spp`), edge-linked equal (`-q`), and edge-unlinked (`-sp`). Terrace-aware search (Chernomor, von Haeseler, and Minh 2016, https://doi.org/10.1093/sysbio/syw037) handles partitioned data with missing taxa.

### RAxML-NG partition format

```
GTR+G, gene1 = 1-500
HKY+I, gene2 = 501-1200
```

Column-range-based partitioning of a concatenated alignment. Kozlov et al. 2019, https://doi.org/10.1093/bioinformatics/btz305.

### BEAST XML

Full model hierarchy in XML: sequence data, substitution models, clock models, tree priors, operators. Complete but verbose. Drummond and Rambaut 2007, https://doi.org/10.1186/1471-2148-7-214.

## Proposed design

YAML or JSON configuration file specifying:

- Partition list: each with alignment path, alphabet, substitution model
- Shared parameters: tree path, clock model, coalescent model
- Output specification

This leverages the existing v1 partition architecture (`enum PartitionSet` of dense or sparse partitions) which already supports multiple partitions on the same tree.

## Implemented part

The `partition` setting of ancestral, optimize and timetree (flag `--partition`, one entry per flag, or a list under the config key `partition`) describes each partition as a set of sites of one `--alignment` file with its own model ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)):

```yaml
alignment:
  - genome.fasta.xz
partition:
  - {name: cp1, sites: ['101-1600\3']}
  - {name: cp2, sites: ['102-1600\3']}
  - {name: cp3, sites: ['103-1600\3'], model: hky85, model_params: ['kappa=4']}
  - {name: noncoding}
```

- **Sites**: site ranges `a`, `a-b` and `a-b\k`, numbered from 1, in the syntax of RAxML-NG partition files and NEXUS charsets. An entry without `sites` takes the sites that the other entries of its alignment leave, and every site belongs to one entry
- **Models**: each entry carries its own `model` and `model_params`, so `--model` and `--model-params` cannot be given together with entries
- **Shared branch lengths**: the partitions share the branch lengths, with one relative rate per partition (proportional scaling, as IQ-TREE `-p`); indel events and rates stay per file
- **Outputs**: sequence outputs stay one file per alignment, in the coordinates of the file

The apps edit the entries in the settings form, and one Rust resolver checks them for the CLI and the apps.

## Implementation pointers

Open design questions:

- YAML vs JSON vs TOML? The project already uses `serde` for all three. YAML is human-friendliest for nested config. JSON is simplest to validate. TOML is used for Cargo but less natural for lists-of-partitions.
- Should the config file replace CLI flags or complement them? IQ-TREE uses partition file alongside CLI flags. BEAST uses XML for everything. A hybrid (config file for partitions, CLI for global options) is the lightest integration.
- Edge-linked vs edge-unlinked partitions: should partitions share branch lengths (proportional scaling) or have independent branch lengths? IQ-TREE supports both via `-q` (linked) and `-sp` (unlinked). All partitions currently share the branch lengths, with one relative rate per partition (proportional scaling, as IQ-TREE `-p`).


## Related

### Known issues

- [N-io-sequence-name-matching-unreliable](../issues/N-io-sequence-name-matching-unreliable.md) -- config file could include explicit name mappings
- [N-io-large-dataset-memory-constraint](../issues/N-io-large-dataset-memory-constraint.md) -- config could specify lazy loading options
- [N-io-multi-segment-genome-input](../issues/N-io-multi-segment-genome-input.md) -- segments with different trees cannot share a run

### Proposals

- [unified-input-format-support](unified-input-format-support.md) -- alternative input path via Auspice/MAT/PhyloXML

## References

1. <a id="ref-1"></a> Maddison, David R., David L. Swofford, and Wayne P. Maddison. 1997. "NEXUS: An Extensible File Format for Systematic Information." _Systematic Biology_ 46(4):590-621. https://doi.org/10.1093/sysbio/46.4.590
2. <a id="ref-2"></a> Ronquist, Fredrik, and John P. Huelsenbeck. 2003. "MrBayes 3: Bayesian Phylogenetic Inference under Mixed Models." _Bioinformatics_ 19(12):1572-1574. https://doi.org/10.1093/bioinformatics/btg180
3. <a id="ref-3"></a> Lanfear, Robert, Brett Calcott, Simon Y. W. Ho, and Stephane Guindon. 2012. "PartitionFinder: Combined Selection of Partitioning Schemes and Substitution Models for Phylogenetic Analyses." _Molecular Biology and Evolution_ 29(6):1695-1701. https://doi.org/10.1093/molbev/mss020
4. <a id="ref-4"></a> Chernomor, Olga, Arndt von Haeseler, and Bui Quang Minh. 2016. "Terrace Aware Data Structure for Phylogenomic Inference from Supermatrices." _Systematic Biology_ 65(6):997-1008. https://doi.org/10.1093/sysbio/syw037
5. <a id="ref-5"></a> Minh, Bui Quang, Heiko A. Schmidt, Olga Chernomor, Dominik Schrempf, Michael D. Woodhams, Arndt von Haeseler, and Robert Lanfear. 2020. "IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era." _Molecular Biology and Evolution_ 37(5):1530-1534. https://doi.org/10.1093/molbev/msaa015
6. <a id="ref-6"></a> Kozlov, Alexey M., Diego Darriba, Tomas Flouri, Benoit Morel, and Alexandros Stamatakis. 2019. "RAxML-NG: A Fast, Scalable and User-Friendly Tool for Maximum Likelihood Phylogenetic Inference." _Bioinformatics_ 35(21):4453-4455. https://doi.org/10.1093/bioinformatics/btz305
7. <a id="ref-7"></a> Drummond, Alexei J., and Andrew Rambaut. 2007. "BEAST: Bayesian Evolutionary Analysis by Sampling Trees." _BMC Evolutionary Biology_ 7:214. https://doi.org/10.1186/1471-2148-7-214
