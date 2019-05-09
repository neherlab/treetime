## Inference of transition between discrete characters and 'mugration' models

Metadata often consists of discrete characters such as 'country' or 'host' and transitions between these states happen along branches of the phylogenetic tree.
If these transitions are plausibly modeled as a time reversible process with comparable sampling probabilities of the different states, these transitions can be treated as if they were mutations between different sequence states.
Due to its analogy to mutation, this type of inference of migration dynamics between geographic regions is often called "mugration".
TreeTime can be re-purposed to do mugration inference by coding discrete states as single letters in "sequences" of length one.

Such 'mugration' model inference is implemented as a special subcommand that is called as
```bash
treetime mugration --tree data/zika/zika.nwk --states data/zika/zika.metadata.csv --attribute country
```
The last parameter `--attribute country` specifies that the column 'country' in the metadata table `zika.metadata.csv` is to be used as discrete character.

This command will produce an annotated nexus tree with state of the attribute added as comment to each node (for example `[&country="brazil"]`).
In addition, an inferred GTR model is between the different states is written to file.

### Marginal distributions of ancestral states
If the additional flag `--confidence` is added to the command, TreeTime will output a file with the inferred probability of finding an internal node in a particular discrete state.
A typical output would be
```
#name     american_samoa  brazil  china   colombia ...
NODE_00   0.0001          0.0003  0.0002  0        ...
NODE_05   0.594           0       0.406   0        ...
```
Note, however, that these probabilities depend strongly on the model that TreeTime inferred to estimate the ancestral states.
Biased sampling of different states (e.g. a human case might be sampled with higher probability than a bird case) violate the model assumptions and will produce unreliable inferences.

### Specifying equilibrium frequencies
To partly address the problems associated with biased sampling, the user can specify the equilibrium frequencies using the flag `--weights`.
This parameter expects a csv or tsv file specifying the relative weights for each discrete state (they will be normalized to 1.0).
These weights correspond to equilibrium frequencies in a time-reversible model.
This has sometimes slightly bizarre implications and should be used with caution.

### Complete set of options
```
\:> treetime mugration -h
usage: TreeTime: Maximum Likelihood Phylodynamics mugration
       [-h] --tree TREE [--attribute ATTRIBUTE] --states STATES
       [--weights WEIGHTS] [--confidence] [--pc PC]
       [--missing-data MISSING_DATA] [--verbose VERBOSE] [--outdir OUTDIR]

Reconstructs discrete ancestral states, for example geographic location, host,
or similar. In addition to ancestral states, a GTR model of state transitions
is inferred.

optional arguments:
  -h, --help            show this help message and exit
  --tree TREE           Name of file containing the tree in newick, nexus, or
                        phylip format. If none is provided, treetime will
                        attempt to build a tree from the alignment using
                        fasttree, iqtree, or raxml (assuming they are
                        installed)
  --attribute ATTRIBUTE
                        attribute to reconstruct, e.g. country
  --states STATES       csv or tsv file with discrete characters.
                        #name,country,continent taxon1,micronesia,oceania ...
  --weights WEIGHTS     csv or tsv file with probabilities of that a randomly
                        sampled sequence at equilibrium has a particular
                        state. E.g. population of different continents or
                        countries. E.g.: #country,weight micronesia,0.1 ...
  --confidence          output confidence of mugration inference
  --pc PC               pseudo-counts higher numbers will results in 'flatter'
                        models
  --missing-data MISSING_DATA
                        string indicating missing data
  --verbose VERBOSE     verbosity of output 0-6
  --outdir OUTDIR       directory to write the output to
```