## Ancestral sequence reconstruction using TreeTime

At the core of TreeTime is a class that models how sequences change along the tree.
This class allows to reconstruct likely sequences of internal nodes of the tree.
On the command-line, ancestral reconstruction can be done via the command
```bash
treetime ancestral --aln data/h3n2_na/h3n2_na_20.fasta --tree data/h3n2_na/h3n2_na_20.nwk --outdir ancestral_results
```
This command will generate the output
```
Inferred GTR model:
Substitution rate (mu): 1.0

Equilibrium frequencies (pi_i):
  A: 0.2983
  C: 0.1986
  G: 0.2353
  T: 0.2579
  -: 0.01

Symmetrized rates from j->i (W_ij):
  A C G T -
  A 0 0.8273  2.8038  0.4525  1.031
  C 0.8273  0 0.5688  2.8435  1.0561
  G 2.8038  0.5688  0 0.6088  1.0462
  T 0.4525  2.8435  0.6088  0 1.0418
  - 1.031 1.0561  1.0462  1.0418  0

Actual rates from j->i (Q_ij):
  A C G T -
  A 0 0.2468  0.8363  0.135 0.3075
  C 0.1643  0 0.1129  0.5646  0.2097
  G 0.6597  0.1338  0 0.1432  0.2462
  T 0.1167  0.7332  0.157 0 0.2686
  - 0.0103  0.0106  0.0105  0.0104  0

--- alignment including ancestral nodes saved as
   ancestral_results/ancestral_sequences.fasta

--- tree saved in nexus format as
   ancestral_results/annotated_tree.nexus
```
TreeTime has inferred a GTR model and used it to reconstruct the most likely ancestral sequences.
The reconstructed sequences will be written to a file ending in `_ancestral.fasta` and a tree with mutations mapped to branches will be saved in nexus format in a file ending on `_mutation.nexus`.
Mutations are added as comments to the nexus file like `[&mutations="G27A,A58G,A745G,G787A,C1155T,G1247A,G1272A"]`.

### Amino acid sequences
Ancestral reconstruction of amino acid sequences works analogously to nucleotide sequences.
However, the user has to either explicitly choose an amino acid substitution model (JTT92)
```bash
treetime ancestral --tree data/h3n2_na/h3n2_na_20.nwk  --aln data/h3n2_na/h3n2_na_20_aa.fasta --gtr JTT92
```
or specify that this is a protein sequence alignment using the flag `--aa`:
```bash
treetime ancestral --tree data/h3n2_na/h3n2_na_20.nwk  --aln data/h3n2_na/h3n2_na_20_aa.fasta --aa
```


### VCF files as input
In addition to standard fasta files, TreeTime can ingest sequence data in form of vcf files which is common for bacterial data sets where short reads are mapped against a reference and only variable sites are reported.
In this case, an additional argument specifying the mapping reference is required.
```bash
treetime ancestral --aln data/tb/lee_2015.vcf.gz --vcf-reference data/tb/tb_ref.fasta --tree data/tb/lee_2015.nwk
```
The ancestral reconstruction is saved as a vcf files with the name `ancestral_sequences.vcf`.

### Additional options

 * `--marginal`: By default, TreeTime will determine the *jointly* most likely sequences. In some cases, the *marginally* most likely sequences (that is after summing over all possible sequences at all other nodes) are more appropriate. Typically marginal and joint reconstruction will give similar results.
 * `--gtr` and `--gtr-params`: By default, TreeTime will infer a GTR model. A few default ones can be specified via command-line arguments, for example `--gtr K80 --gtr-params kappa=0.2 pis=0.3,0.25,0.25,0.2` would specify the model K80 with a tv/ts ratio of 0.2 and a set of nucleotide frequencies.
 * `--report-ambiguous`: by default, sequence changes involving characters indicating ambiguous states (like 'N') are not included in the output. To include them, add this flag.

Documentation of the complete set of options is available by typing
```bash
treetime ancestral -h
```
which yields
```
usage: TreeTime: Maximum Likelihood Phylodynamics ancestral
       [-h] --aln ALN [--vcf-reference VCF_REFERENCE] [--tree TREE]
       [--gtr GTR] [--gtr-params GTR_PARAMS [GTR_PARAMS ...]] [--marginal]
       [--keep-overhangs] [--zero-based] [--report-ambiguous]
       [--verbose VERBOSE] [--outdir OUTDIR]

Reconstructs ancestral sequences and maps mutations to the tree. The output
consists of a file ending with _ancestral.fasta with ancestral sequences and a
tree ending with _mutation.nexus with mutations added as comments like
_A45G_..., number in SNPs used 1-based index by default. The inferred GTR
model is written to stdout.

optional arguments:
  -h, --help            show this help message and exit
  --aln ALN             alignment file (fasta)
  --vcf-reference VCF_REFERENCE
                        only for vcf input: fasta file of the sequence the VCF
                        was mapped to.
  --tree TREE           Name of file containing the tree in newick, nexus, or
                        phylip format. If none is provided, treetime will
                        attempt to build a tree from the alignment using
                        fasttree, iqtree, or raxml (assuming they are
                        installed)
  --gtr GTR             GTR model to use. '--gtr infer' will infer a model
                        from the data. Alternatively, specify the model type.
                        If the specified model requires additional options,
                        use '--gtr-params' to specify those.
  --gtr-params GTR_PARAMS [GTR_PARAMS ...]
                        GTR parameters for the model specified by the --gtr
                        argument. The parameters should be feed as 'key=value'
                        list of parameters. Example: '--gtr K80 --gtr-params
                        kappa=0.2 pis=0.25,0.25,0.25,0.25'. See the exact
                        definitions of the parameters in the GTR creation
                        methods in treetime/nuc_models.py or
                        treetime/aa_models.py
  --aa                  use aminoacid alphabet
  --marginal            marginal reconstruction of ancestral sequences
  --keep-overhangs      do not fill terminal gaps
  --zero-based          zero based mutation indexing
  --report-ambiguous    include transitions involving ambiguous states
  --verbose VERBOSE     verbosity of output 0-6
  --outdir OUTDIR       directory to write the output to

```
