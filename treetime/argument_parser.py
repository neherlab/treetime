#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
import sys, argparse, os
from treetime.wrappers import ancestral_reconstruction, mugration, scan_homoplasies, timetree, estimate_clock_model
import treetime

py2 = sys.version_info.major==2

def set_default_subparser(self, name, args=None, positional_args=0):
    """default subparser selection. Call after setup, just before parse_args()
    name: is the name of the subparser to call by default
    args: if set is the argument list handed to parse_args()
    https://stackoverflow.com/questions/6365601/default-sub-command-or-handling-no-sub-command-with-argparse
    """
    subparser_found = False
    if len(sys.argv)==1:
        sys.argv.append('-h')
    else:
        for x in self._subparsers._actions:
            if not isinstance(x, argparse._SubParsersAction):
                continue
            for sp_name in x._name_parser_map.keys():
                if sp_name in sys.argv[1:]:
                    subparser_found = True
        if not subparser_found:
            # insert default subcommand in first position
            if args is None:
                sys.argv.insert(1, name)
            else:
                args.insert(1, name)


if py2:
    argparse.ArgumentParser.set_default_subparser = set_default_subparser


treetime_description = \
    "TreeTime: Maximum Likelihood Phylodynamics\n\n"
subcommand_description = \
    "In addition, TreeTime implements several sub-commands:\n\n"\
    "\t ancestral\tinfer ancestral sequences maximizing the joint or marginal likelihood.\n"\
    "\t homoplasy\tanalyze patterns of recurrent mutations aka homoplasies.\n"\
    "\t clock\t\testimate molecular clock parameters and reroot the tree.\n"\
    "\t mugration\tmap discrete character such as host or country to the tree.\n\n"\
    "(note that 'tt' is a default subcommand in python2 that doesn't need to be specified).\n"\
    "To print a description and argument list of the individual sub-commands, type:\n\n"\
    "\t treetime <subcommand> -h\n\n"

ref_msg = \
    "If you use results from treetime in a publication, please cite:"\
    "\n\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"\
    "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n"

timetree_description=\
    "TreeTime infers a time scaled phylogeny given a tree topology, an alignment, "\
    "and tip dates. Reconstructs ancestral sequences and infers a molecular clock tree. "\
    "TreeTime will reroot the tree and resolve polytomies by default. "\
    "In addition, treetime will infer ancestral sequences and a GTR substitution model. "\
    "Inferred mutations are included as comments in the output tree.\n\n"

gtr_description = "GTR model to use. '--gtr infer' will infer a model "\
    "from the data. Alternatively, specify the model type. If the specified model "\
    "requires additional options, use '--gtr-params' to specify those."

gtr_params_description =  "GTR parameters for the model specified by "\
    "the --gtr argument. The parameters should be feed as 'key=value' "\
    "list of parameters. Example: '--gtr K80 --gtr-params kappa=0.2 "\
    "pis=0.25,0.25,0.25,0.25'. See the exact definitions of the "\
    "parameters in the GTR creation methods in treetime/nuc_models.py "\
    "or treetime/aa_models.py"

reroot_description = "Reroot the tree using root-to-tip regression. Valid choices are "\
    "'min_dev', 'least-squares', and 'oldest'. 'least-squares' adjusts the root to "\
    "minimize residuals of the root-to-tip vs sampling time regression, " \
    "'min_dev' minimizes variance of root-to-tip distances. "\
    "'least-squares' can be combined with --covariation to account for shared ancestry. "\
    "Alternatively, you can specify a node name or a list of node names "\
    "to be used as outgroup or use 'oldest' to reroot to the oldest node. "\
    "By default, TreeTime will reroot using 'least-squares'. "\
    "Use --keep-root to keep the current root."

tree_description = "Name of file containing the tree in "\
    "newick, nexus, or phylip format. If none is provided, "\
    "treetime will attempt to build a tree from the alignment "\
    "using fasttree, iqtree, or raxml (assuming they are installed)"

aln_description = "alignment file (fasta)"

dates_description = "csv file with dates for nodes with 'node_name, date' where date is float (as in 2012.15)"

coalescent_description = \
    "coalescent time scale -- sensible values are on the order of the average "\
    "hamming distance of contemporaneous sequences. In addition, 'opt' "\
    "'skyline' are valid options and estimate a constant coalescent rate "\
    "or a piecewise linear coalescent rate history"

ancestral_description = \
    "Reconstructs ancestral sequences and maps mutations to the tree. "\
    "The output consists of a file 'ancestral.fasta' with ancestral sequences "\
    "and a tree 'annotated_tree.nexus' with mutations added as comments "\
    "like A45G,G136T,..., number in SNPs used 1-based index by default. "\
    "The inferred GTR model is written to stdout."

homoplasy_description = \
    "Reconstructs ancestral sequences and maps mutations to the tree. "\
    "The tree is then scanned for homoplasies. An excess number of homoplasies "\
    "might suggest contamination, recombination, culture adaptation or similar."

mugration_description = \
    "Reconstructs discrete ancestral states, for example "\
    "geographic location, host, or similar. In addition to ancestral states, "\
    "a GTR model of state transitions is inferred."

def add_seq_len_aln_group(parser):
    parser.add_argument('--sequence-length', type=int, help="length of the sequence, "
                              "used to calculate expected variation in branch length. "
                              "Not required if alignment is provided.")
    add_aln_group(parser, required=False)
    # seq_group_ex.add_argument('--aln',  type=str, help=aln_description)

def add_aln_group(parser, required=True):
    parser.add_argument('--aln', required=required, type=str, help=aln_description)
    parser.add_argument('--vcf-reference', type=str, help='only for vcf input: fasta file of the sequence the VCF was mapped to.')


def add_reroot_group(parser):
    parser.add_argument('--clock-filter', type=float, default=3,
                              help="ignore tips that don't follow a loose clock, "
                                   "'clock-filter=number of interquartile ranges from regression'. "
                                   "Default=3.0, set to 0 to switch off.")
    reroot_group = parser.add_mutually_exclusive_group()
    reroot_group.add_argument('--reroot', nargs='+', default='best', help=reroot_description)
    reroot_group.add_argument('--keep-root', required = False, action="store_true", default=False,
            help ="don't reroot the tree. Otherwise, reroot to minimize the "
                  "the residual of the regression of root-to-tip distance and sampling time")
    parser.add_argument('--tip-slack', type=float, default=3,
                              help="excess variance associated with terminal nodes accounting for "
                                   " overdisperion of the molecular clock")
    parser.add_argument('--covariation', action='store_true', help="Account for covariation when estimating rates "
                        "or rerooting using root-to-tip regression, default False.")

def add_gtr_arguments(parser):
    parser.add_argument('--gtr', default='infer', help=gtr_description)
    parser.add_argument('--gtr-params', nargs='+', help=gtr_params_description)
    parser.add_argument('--aa', action='store_true', help="use aminoacid alphabet")


def add_anc_arguments(parser):
    parser.add_argument('--keep-overhangs', default = False, action='store_true', help='do not fill terminal gaps')
    parser.add_argument('--zero-based', default = False, action='store_true', help='zero based mutation indexing')
    parser.add_argument('--report-ambiguous', default=False, action="store_true", help='include transitions involving ambiguous states')


def add_common_args(parser):
    parser.add_argument('--verbose', default=1, type=int,  help='verbosity of output 0-6')
    parser.add_argument('--outdir', type=str,  help='directory to write the output to')


def make_parser():
    parser = argparse.ArgumentParser(description = "",
                                     usage=treetime_description)

    subparsers = parser.add_subparsers()

    if py2:
        t_parser = subparsers.add_parser('tt', description=timetree_description)
    else:
        t_parser = parser
    t_parser.add_argument('--tree', type=str, help=tree_description)
    add_seq_len_aln_group(t_parser)
    t_parser.add_argument('--dates', type=str, help=dates_description)
    add_reroot_group(t_parser)
    add_gtr_arguments(t_parser)
    t_parser.add_argument('--clock-rate', type=float, help="if specified, the rate of the molecular clock won't be optimized.")
    t_parser.add_argument('--clock-std-dev', type=float, help="standard deviation of the provided clock rate estimate")
    t_parser.add_argument('--branch-length-mode', default='auto', type=str, choices=['auto', 'input', 'joint', 'marginal'],
                        help="If set to 'input', the provided branch length will be used without modification. "
                             "Note that branch lengths optimized by treetime are only accurate at short evolutionary distances.")
    t_parser.add_argument('--confidence', action='store_true', help="estimate confidence intervals of divergence times.")
    t_parser.add_argument('--keep-polytomies', default=False, action='store_true',
                        help="Don't resolve polytomies using temporal information.")
    t_parser.add_argument('--relax',nargs=2, type=float,
                        help='use an autocorrelated molecular clock. Strength of the gaussian priors on'
                             ' branch specific rate deviation and the coupling of parent and offspring'
                             ' rates can be specified e.g. as --relax 1.0 0.5. Values around 1.0 correspond'
                             ' to weak priors, larger values constrain rate deviations more strongly.'
                             ' Coupling 0 (--relax 1.0 0) corresponds to an un-correlated clock.')
    t_parser.add_argument('--max-iter', default=2, type=int,
                        help='maximal number of iterations the inference cycle is run. Note that for polytomy resolution and coalescence models max_iter should be at least 2')
    t_parser.add_argument('--coalescent', default="0.0", type=str,
                          help=coalescent_description)
    t_parser.add_argument('--plot-tree', default="timetree.pdf",
                            help = "filename to save the plot to. Suffix will determine format"
                                   " (choices pdf, png, svg, default=pdf)")
    t_parser.add_argument('--plot-rtt', default="root_to_tip_regression.pdf",
                            help = "filename to save the plot to. Suffix will determine format"
                                   " (choices pdf, png, svg, default=pdf)")
    t_parser.add_argument('--tip-labels', action='store_true',
                            help = "add tip labels (default for small trees with <30 leaves)")
    t_parser.add_argument('--no-tip-labels', action='store_true',
                            help = "don't show tip labels (default for small trees with >=30 leaves)")
    add_anc_arguments(t_parser)
    add_common_args(t_parser)

    def toplevel(params):
        if (params.aln or params.tree) and params.dates:
            timetree(params)
        else:
            print(treetime_description+timetree_description+subcommand_description+
                  "'--dates' and '--aln' or '--tree' are REQUIRED inputs, type 'treetime -h' for a full list of arguments.\n")

    t_parser.set_defaults(func=toplevel)


    ## HOMOPLASY SCANNER
    h_parser = subparsers.add_parser('homoplasy', description=homoplasy_description)
    add_aln_group(h_parser)
    h_parser.add_argument('--tree', type = str,  help=tree_description)
    h_parser.add_argument('--const', type = int, default=0, help ="number of constant sites not included in alignment")
    h_parser.add_argument('--rescale', type = float, default=1.0, help ="rescale branch lengths")
    h_parser.add_argument('--detailed', required = False, action="store_true",  help ="generate a more detailed report")
    add_gtr_arguments(h_parser)
    h_parser.add_argument('--zero-based', default = False, action='store_true', help='zero based mutation indexing')
    h_parser.add_argument('-n', default = 10, type=int, help='number of mutations/nodes that are printed to screen')
    h_parser.add_argument('--drms', type=str, help='TSV file containing DRM info. columns headers: GENOMIC_POSITION, ALT_BASE, DRUG, GENE, SUBSTITUTION')
    add_common_args(h_parser)
    h_parser.set_defaults(func=scan_homoplasies)

    ## ANCESTRAL RECONSTRUCTION
    a_parser = subparsers.add_parser('ancestral', description=ancestral_description)
    add_aln_group(a_parser)
    a_parser.add_argument('--tree', type = str,  help =tree_description)
    add_gtr_arguments(a_parser)
    a_parser.add_argument('--marginal', default = False, action="store_true", help ="marginal reconstruction of ancestral sequences")
    add_anc_arguments(a_parser)
    add_common_args(a_parser)
    a_parser.set_defaults(func=ancestral_reconstruction)

    ## MUGRATION
    m_parser = subparsers.add_parser('mugration', description=mugration_description)
    m_parser.add_argument('--tree', required = True, type=str, help=tree_description)
    m_parser.add_argument('--attribute', type=str, help ="attribute to reconstruct, e.g. country")
    m_parser.add_argument('--states', required = True, type=str, help ="csv or tsv file with discrete characters."
                                    "\n#name,country,continent\ntaxon1,micronesia,oceania\n...")
    m_parser.add_argument('--weights', type=str, help="csv or tsv file with probabilities of that a randomly sampled "
                        "sequence at equilibrium has a particular state. E.g. population of different continents or countries. E.g.:"
                        "\n#country,weight\nmicronesia,0.1\n...")
    m_parser.add_argument('--confidence', action="store_true", help="output confidence of mugration inference")
    m_parser.add_argument('--pc', type=float, default=1.0, help ="pseudo-counts higher numbers will results in 'flatter' models")
    m_parser.add_argument('--missing-data', type=str, default='?', help ="string indicating missing data")
    m_parser.add_argument('--sampling-bias-correction', type=float,
                        help='a rough estimate of how many more events would have been observed'
                             ' if sequences represented an even sample. This should be'
                             ' roughly the (1-sum_i p_i^2)/(1-sum_i t_i^2), where p_i'
                             ' are the equilibrium frequencies and t_i are apparent ones.'
                             '(or rather the time spent in a particular state on the tree)')
    add_common_args(m_parser)
    m_parser.set_defaults(func=mugration)


    ## CLOCKSIGNAL
    c_parser = subparsers.add_parser('clock',
            description="Calculates the root-to-tip regression and quantifies the 'clock-i-ness' of the tree. "
                        "It will reroot the tree to maximize the clock-like "
                        "signal and recalculate branch length unless run with --keep_root.")
    c_parser.add_argument('--tree', required=True, type=str,  help=tree_description)
    c_parser.add_argument('--dates', required=True, type=str, help=dates_description)
    add_seq_len_aln_group(c_parser)

    add_reroot_group(c_parser)
    c_parser.add_argument('--allow-negative-rate', required = False, action="store_true", default=False,
                          help="By default, rates are forced to be positive. For trees with little temporal "
                               "signal it is advisable to remove this restriction to achieve essentially mid-point rooting.")
    c_parser.add_argument('--plot-rtt', default="root_to_tip_regression.pdf",
                            help = "filename to save the plot to. Suffix will determine format"
                                   " (choices pdf, png, svg, default=pdf)")
    add_common_args(c_parser)
    c_parser.set_defaults(func=estimate_clock_model)

    # make a version subcommand
    v_parser = subparsers.add_parser('version', description='print version')
    v_parser.set_defaults(func=lambda x: print(treetime.version))

    ## call the relevant function and return
    if py2:
        parser.set_default_subparser('tt')

    return parser
