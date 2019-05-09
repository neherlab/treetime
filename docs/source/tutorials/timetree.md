## Estimation of time scaled phylogenies

The principal functionality of TreeTime is estimating time trees from an initial tree topology, a set of date constraints (e.g. tip dates), and an alignment (optional)
```bash
treetime --tree data/h3n2_na/h3n2_na_20.nwk --dates data/h3n2_na/h3n2_na_20.metadata.csv --aln data/h3n2_na/h3n2_na_20.fasta --outdir h3n2_timetree
```
This command will estimate an GTR model, a molecular clock model, and a time-stamped phylogeny.
Most results are saved to file, but rate and GTR model are printed to the console:
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
	A	C	G	T	-
  A	0	0.8273	2.8038	0.4525	1.031
  C	0.8273	0	0.5688	2.8435	1.0561
  G	2.8038	0.5688	0	0.6088	1.0462
  T	0.4525	2.8435	0.6088	0	1.0418
  -	1.031	1.0561	1.0462	1.0418	0

Actual rates from j->i (Q_ij):
	A	C	G	T	-
  A	0	0.2468	0.8363	0.135	0.3075
  C	0.1643	0	0.1129	0.5646	0.2097
  G	0.6597	0.1338	0	0.1432	0.2462
  T	0.1167	0.7332	0.157	0	0.2686
  -	0.0103	0.0106	0.0105	0.0104	0

Root-Tip-Regression:
 --rate:	2.613e-03
 --chi^2:	22.16

 --r^2:	0.98

--- saved tree as
	 h3n2_timetree/timetree.pdf

--- alignment including ancestral nodes saved as
	 h3n2_timetree/ancestral_sequences.fasta

--- saved divergence times in
	 h3n2_timetree/dates.tsv

--- tree saved in nexus format as
	 h3n2_timetree/timetree.nexus

--- root-to-tip plot saved to
	h3n2_timetree/root_to_tip_regression.pdf
```
The script saved an alignment with reconstructed ancestral sequences, an annotated tree in nexus format in which branch length correspond to years and mutations and node dates are added as comments to each node.
In addition, the root-to-tip vs time regression and the tree are drawn and saved to file.

![rtt](figures/timetree.png)

### Fixed evolutionary rate
If the temporal signal in the data is weak and the clock rate can't be estimated confidently from the data, it is advisable to specify the rate explicitly.
This can be done using the argument
```
--clock-rate <rate>
```

### Specify or estimate coalescent models
TreeTime can be run either without a tree prior or with a Kingman coalescent tree prior.
The later is parameterized by a time scale 'Tc' which can vary in time.
This time scale is often called 'effective population size' Ne, but the appropriate Tc has very little to do with census population sizes.
To activate the Kingman Coalescent model in TreeTime, you need to add the flag
```
 --coalescent <arg>
```
where the argument is either a floating point number giving the time scale of coalescence in units of divergence, 'const' to have TreeTime estimate a constant merger rate, or 'skyline'.
In the latter case, TreeTime will estimate a piece-wise linear merger rate trajectory and save this in files ending on 'skyline.tsv' and 'skyline.pdf'

The following command will run TreeTime on the ebola example data set and estimate a time tree along with a skyline (this will take a few minutes).
```bash
treetime --tree data/ebola/ebola.nwk --dates data/ebola/ebola.metadata.csv --aln data/ebola/ebola.fasta --outdir ebola  --coalescent skyline
```
![ebola_skyline](figures/ebola_skyline.png)


### Confidence intervals
In its default setting, `treetime` does not estimate confidence intervals of divergence times.
Such estimates require calculation of the marginal probability distributions of the dates of the internal nodes that take about 2-3 times as long as calculating only the jointly maximally likely dates.
To switch on confidence estimation, pass the flag `--confidence`.
TreeTime will run another round of marginal timetree reconstruction and determine the region that contains 90% of the marginal probability distribution of the node dates.
These intervals are drawn into the tree graph and written to the dates file.

### VCF files as input
In addition to standard fasta files, TreeTime can ingest sequence data in form of vcf files which is common for bacterial data sets where short reads are mapped against a reference and only variable sites are reported.
The following example with a set of MtB sequences uses a fixed evolutionary rate of 1e-7 per site and year.
```bash
treetime --aln data/tb/lee_2015.vcf.gz --vcf-reference data/tb/tb_ref.fasta --tree data/tb/lee_2015.nwk --clock-rate 1e-7 --dates data/tb/lee_2015.metadata.tsv
```
For many bacterial data set were the temporal signal in the data is weak, it is advisable to fix the rate of the molecular clock explicitly.
Divergence times, however, will depend on this choice.

### Additional options
Documentation of the complete set of options is available by typing
```bash
treetime -h
```
which yields
```
usage: TreeTime: Maximum Likelihood Phylodynamics

positional arguments:
  {homoplasy,ancestral,mugration,clock,version}

optional arguments:
  -h, --help            show this help message and exit
  --tree TREE           Name of file containing the tree in newick, nexus, or
                        phylip format. If none is provided, treetime will
                        attempt to build a tree from the alignment using
                        fasttree, iqtree, or raxml (assuming they are
                        installed)
  --sequence-length SEQUENCE_LENGTH
                        length of the sequence, used to calculate expected
                        variation in branch length. Not required if alignment
                        is provided.
  --aln ALN             alignment file (fasta)
  --vcf-reference VCF_REFERENCE
                        only for vcf input: fasta file of the sequence the VCF
                        was mapped to.
  --dates DATES         csv file with dates for nodes with 'node_name, date'
                        where date is float (as in 2012.15)
  --clock-filter CLOCK_FILTER
                        ignore tips that don't follow a loose clock, 'clock-
                        filter=number of interquartile ranges from regression'
  --reroot REROOT       reroot the tree. Valid choices are 'ML', 'ML-rough',
                        'least-squares', 'min_dev', 'midpoint' or a node name
                        to be used as outgroup. Use --keep-root to keep the
                        current root.
  --keep-root           don't reroot the tree. Otherwise, reroot to minimize
                        the the residual of the regression of root-to-tip
                        distance and sampling time
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
  --clock-rate CLOCK_RATE
                        if specified, the rate of the molecular clock won't be
                        optimized.
  --clock-std-dev CLOCK_STD_DEV
                        standard deviation of the provided clock rate estimate
  --branch-length-mode {auto,input,joint,marginal}
                        If set to 'input', the provided branch length will be
                        used without modification. Note that branch lengths
                        optimized by treetime are only accurate at short
                        evolutionary distances.
  --confidence          estimate confidence intervals of divergence times.
  --keep-polytomies     Don't resolve polytomies using temporal information.
  --relax [RELAX [RELAX ...]]
                        use an autocorrelated molecular clock. Prior strength
                        and coupling of parent and offspring rates can be
                        specified e.g. as --relax 1.0 0.5
  --max-iter MAX_ITER   maximal number of iterations the inference cycle is
                        run. Note that for polytomy resolution and coalescence
                        models max_iter should be at least 2
  --coalescent COALESCENT
                        coalescent time scale -- sensible values are on the
                        order of the average hamming distance of
                        contemporaneous sequences. In addition, 'opt'
                        'skyline' are valid options and estimate a constant
                        coalescent rate or a piecewise linear coalescent rate
                        history
  --plot-tree PLOT_TREE
                        filename to save the plot to. Suffix will determine
                        format (choices pdf, png, svg, default=pdf)
  --plot-rtt PLOT_RTT   filename to save the plot to. Suffix will determine
                        format (choices pdf, png, svg, default=pdf)
  --tip-labels          add tip labels (default for small trees with <30
                        leaves)
  --no-tip-labels       don't show tip labels (default for small trees with
                        >=30 leaves)
  --keep-overhangs      do not fill terminal gaps
  --zero-based          zero based mutation indexing
  --report-ambiguous    include transitions involving ambiguous states
  --verbose VERBOSE     verbosity of output 0-6
  --outdir OUTDIR       directory to write the output to
```