from __future__ import print_function, division, absolute_import
import os, shutil, sys
import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo, AlignIO
from Bio import __version__ as bioversion
from treetime import TreeAnc, GTR, TreeTime
from treetime import config as ttconf
from treetime import utils
from treetime.vcf_utils import read_vcf, write_vcf
from treetime.seq_utils import alphabets

def assure_tree(params, tmp_dir='treetime_tmp'):
    """
    Function that attempts to load a tree and build it from the alignment
    if no tree is provided.
    """
    if params.tree is None:
        params.tree = os.path.basename(params.aln)+'.nwk'
        print("No tree given: inferring tree")
        utils.tree_inference(params.aln, params.tree, tmp_dir = tmp_dir)

    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)

    try:
        tt = TreeAnc(params.tree)
    except:
        print("Tree loading/building failed.")
        return 1
    return 0

def create_gtr(params):
    """
    parse the arguments referring to the GTR model and return a GTR structure
    """
    model = params.gtr
    gtr_params = params.gtr_params
    if model == 'infer':
        gtr = GTR.standard('jc', alphabet='aa' if params.aa else 'nuc')
    else:
        try:
            kwargs = {}
            if gtr_params is not None:
                for param in gtr_params:
                    keyval = param.split('=')
                    if len(keyval)!=2: continue
                    if keyval[0] in ['pis', 'pi', 'Pi', 'Pis']:
                        keyval[0] = 'pi'
                        keyval[1] = list(map(float, keyval[1].split(',')))
                    elif keyval[0] not in ['alphabet']:
                        keyval[1] = float(keyval[1])
                    kwargs[keyval[0]] = keyval[1]
            else:
                print ("GTR params are not specified. Creating GTR model with default parameters")

            gtr = GTR.standard(model, **kwargs)
            infer_gtr = False
        except:
            print ("Could not create GTR model from input arguments. Using default (Jukes-Cantor 1969)")
            gtr = GTR.standard('jc', alphabet='aa' if params.aa else 'nuc')
            infer_gtr = False
    return gtr

def get_outdir(params, suffix='_treetime'):
    if params.outdir:
        if os.path.exists(params.outdir):
            if os.path.isdir(params.outdir):
                return params.outdir.rstrip('/') + '/'
            else:
                print("designated output location %s is not a directory"%params.outdir, file=stderr)
        else:
            os.makedirs(params.outdir)
            return params.outdir.rstrip('/') + '/'

    from datetime import datetime
    outdir = datetime.now().date().isoformat()+suffix.rstrip('/') + '/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    return outdir

def get_basename(params, outdir):
    # if params.aln:
    #     basename = outdir + '.'.join(params.aln.split('/')[-1].split('.')[:-1])
    # elif params.tree:
    #     basename = outdir + '.'.join(params.tree.split('/')[-1].split('.')[:-1])
    # else:
    basename = outdir
    return basename

def read_in_DRMs(drm_file, offset):
    import pandas as pd

    DRMs = {}
    drmPositions = []

    df = pd.read_csv(drm_file, sep='\t')
    for mi, m in df.iterrows():
        pos = m.GENOMIC_POSITION-1+offset #put in correct numbering
        drmPositions.append(pos)

        if pos in DRMs:
            DRMs[pos]['alt_base'][m.ALT_BASE] = m.SUBSTITUTION
        else:
            DRMs[pos] = {}
            DRMs[pos]['drug'] = m.DRUG
            DRMs[pos]['alt_base'] = {}
            DRMs[pos]['alt_base'][m.ALT_BASE] = m.SUBSTITUTION
            DRMs[pos]['gene'] = m.GENE

    drmPositions = np.array(drmPositions)
    drmPositions = np.unique(drmPositions)
    drmPositions = np.sort(drmPositions)

    DRM_info = {'DRMs': DRMs,
            'drmPositions': drmPositions}

    return DRM_info


def read_if_vcf(params):
    """
    Checks if input is VCF and reads in appropriately if it is
    """
    ref = None
    aln = params.aln
    fixed_pi = None
    if hasattr(params, 'aln') and params.aln is not None:
        if any([params.aln.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
            if not params.vcf_reference:
                print("ERROR: a reference Fasta is required with VCF-format alignments")
                return -1
            compress_seq = read_vcf(params.aln, params.vcf_reference)
            sequences = compress_seq['sequences']
            ref = compress_seq['reference']
            aln = sequences

            if not hasattr(params, 'gtr') or params.gtr=="infer": #if not specified, set it:
                alpha = alphabets['aa'] if params.aa else alphabets['nuc']
                fixed_pi = [ref.count(base)/len(ref) for base in alpha]
                if fixed_pi[-1] == 0:
                    fixed_pi[-1] = 0.05
                    fixed_pi = [v-0.01 for v in fixed_pi]

    return aln, ref, fixed_pi


def plot_rtt(tt, fname):
    tt.plot_root_to_tip()

    from matplotlib import pyplot as plt
    plt.savefig(fname)
    print("--- root-to-tip plot saved to  \n\t"+fname)


def export_sequences_and_tree(tt, basename, is_vcf=False, zero_based=False,
                              report_ambiguous=False, timetree=False, confidence=False):
    seq_info = is_vcf or tt.aln
    if is_vcf:
        tt.recover_var_ambigs()
        outaln_name = basename + 'ancestral_sequences.vcf'
        write_vcf(tt.get_tree_dict(keep_var_ambigs=True), outaln_name)
    elif tt.aln:
        outaln_name = basename + 'ancestral_sequences.fasta'
        AlignIO.write(tt.get_reconstructed_alignment(), outaln_name, 'fasta')
    if seq_info:
        print("\n--- alignment including ancestral nodes saved as  \n\t %s\n"%outaln_name)

    # decorate tree with inferred mutations
    terminal_count = 0
    offset = 0 if zero_based else 1
    if timetree:
        dates_fname = basename + 'dates.tsv'
        fh_dates = open(dates_fname, 'w')
        if confidence:
            fh_dates.write('#Lower and upper bound delineate the 90% max posterior region\n')
            fh_dates.write('#node\tdate\tnumeric date\tlower bound\tupper bound\n')
        else:
            fh_dates.write('#node\tdate\tnumeric date\n')
    for n in tt.tree.find_clades():
        if timetree:
            if confidence:
                conf = tt.get_max_posterior_region(n, fraction=0.9)
                fh_dates.write('%s\t%s\t%f\t%f\t%f\n'%(n.name, n.date, n.numdate,conf[0], conf[1]))
            else:
                fh_dates.write('%s\t%s\t%f\n'%(n.name, n.date, n.numdate))

        n.confidence=None
        # due to a bug in older versions of biopython that truncated filenames in nexus export
        # we truncate them by hand and make them unique.
        if n.is_terminal() and len(n.name)>40 and bioversion<"1.69":
            n.name = n.name[:35]+'_%03d'%terminal_count
            terminal_count+=1
        n.comment=''
        if seq_info and len(n.mutations):
            if report_ambiguous:
                n.comment= '&mutations="' + ','.join([a+str(pos + offset)+d for (a,pos, d) in n.mutations])+'"'
            else:
                n.comment= '&mutations="' + ','.join([a+str(pos + offset)+d for (a,pos, d) in n.mutations
                                                      if tt.gtr.ambiguous not in [a,d]])+'"'
        if timetree:
            n.comment+=(',' if n.comment else '&') + 'date=%1.2f'%n.numdate

    # write tree to file
    fmt_bl = "%1.6f" if tt.seq_len<1e6 else "%1.8e"
    if timetree:
        outtree_name = basename + 'timetree.nexus'
        print("--- saved divergence times in \n\t %s\n"%dates_fname)
        Phylo.write(tt.tree, outtree_name, 'nexus')
    else:
        outtree_name = basename + 'annotated_tree.nexus'
        Phylo.write(tt.tree, outtree_name, 'nexus', format_branch_length=fmt_bl)
    print("--- tree saved in nexus format as  \n\t %s\n"%outtree_name)

    if timetree:
        for n in tt.tree.find_clades():
            n.branch_length = n.mutation_length
        outtree_name = basename + 'divergence_tree.nexus'
        Phylo.write(tt.tree, outtree_name, 'nexus', format_branch_length=fmt_bl)
        print("--- divergence tree saved in nexus format as  \n\t %s\n"%outtree_name)


def print_save_plot_skyline(tt, n_std=2.0, screen=True, save='', plot=''):
    if plot:
        import matplotlib.pyplot as plt

    skyline, conf = tt.merger_model.skyline_inferred(gen=50, confidence=n_std)
    if save: fh = open(save, 'w')
    header1 = "Skyline assuming 50 gen/year and approximate confidence bounds (+/- %f standard deviations of the LH)\n"%n_std
    header2 = "date \tN_e \tlower \tupper"
    if screen: print('\t'+header1+'\t'+header2)
    if save: fh.write("#"+ header1+'#'+header2+'\n')
    for (x,y, y1, y2) in zip(skyline.x, skyline.y, conf[0], conf[1]):
        if screen: print("\t%1.1f\t%1.1f\t%1.1f\t%1.1f"%(x,y, y1, y2))
        if save: fh.write("%1.1f\t%1.1f\t%1.1f\t%1.1f\n"%(x,y, y1, y2))

    if save:
        print("\n --- written skyline to %s\n"%save)
        fh.close()

    if plot:
        plt.figure()
        plt.fill_between(skyline.x, conf[0], conf[1], color=(0.8, 0.8, 0.8))
        plt.plot(skyline.x, skyline.y, label='maximum likelihood skyline')
        plt.yscale('log')
        plt.legend()
        plt.ticklabel_format(axis='x',useOffset=False)
        plt.savefig(plot)


def scan_homoplasies(params):
    """
    the function implementing treetime homoplasies
    """
    if assure_tree(params, tmp_dir='homoplasy_tmp'):
        return 1

    gtr = create_gtr(params)

    ###########################################################################
    ### READ IN VCF
    ###########################################################################
    #sets ref and fixed_pi to None if not VCF
    aln, ref, fixed_pi = read_if_vcf(params)
    is_vcf = True if ref is not None else False

    ###########################################################################
    ### ANCESTRAL RECONSTRUCTION
    ###########################################################################
    treeanc = TreeAnc(params.tree, aln=aln, ref=ref, gtr=gtr, verbose=1,
                      fill_overhangs=True)
    if treeanc.aln is None: # if alignment didn't load, exit
        return 1

    if is_vcf:
        L = len(ref) + params.const
    else:
        L = treeanc.aln.get_alignment_length() + params.const

    N_seq = len(treeanc.aln)
    N_tree = treeanc.tree.count_terminals()
    if params.rescale!=1.0:
        for n in treeanc.tree.find_clades():
            n.branch_length *= params.rescale
            n.mutation_length = n.branch_length

    print("read alignment from file %s with %d sequences of length %d"%(params.aln,N_seq,L))
    print("read tree from file %s with %d leaves"%(params.tree,N_tree))
    print("\ninferring ancestral sequences...")

    ndiff = treeanc.infer_ancestral_sequences('ml', infer_gtr=params.gtr=='infer',
                                      marginal=False, fixed_pi=fixed_pi)
    print("...done.")
    if ndiff==ttconf.ERROR: # if reconstruction failed, exit
        print("Something went wrong during ancestral reconstruction, please check your input files.", file=sys.stderr)
        return 1
    else:
        print("...done.")

    if is_vcf:
        treeanc.recover_var_ambigs()

    ###########################################################################
    ### analysis of reconstruction
    ###########################################################################
    from collections import defaultdict
    from scipy.stats import poisson
    offset = 0 if params.zero_based else 1

    if params.drms:
        DRM_info = read_in_DRMs(params.drms, offset)
        drms = DRM_info['DRMs']

    # construct dictionaries gathering mutations and positions
    mutations = defaultdict(list)
    positions = defaultdict(list)
    terminal_mutations = defaultdict(list)
    for n in treeanc.tree.find_clades():
        if n.up is None:
            continue

        if len(n.mutations):
            for (a,pos, d) in n.mutations:
                if '-' not in [a,d] and 'N' not in [a,d]:
                    mutations[(a,pos+offset,d)].append(n)
                    positions[pos+offset].append(n)
            if n.is_terminal():
                for (a,pos, d) in n.mutations:
                    if '-' not in [a,d] and 'N' not in [a,d]:
                        terminal_mutations[(a,pos+offset,d)].append(n)

    # gather homoplasic mutations by strain
    mutation_by_strain = defaultdict(list)
    for n in treeanc.tree.get_terminals():
        for a,pos,d in n.mutations:
            if pos+offset in positions and len(positions[pos+offset])>1:
                if '-' not in [a,d] and 'N' not in [a,d]:
                    mutation_by_strain[n.name].append([(a,pos+offset,d), len(positions[pos])])


    # total_branch_length is the expected number of substitutions
    # corrected_branch_length is the expected number of observable substitutions
    # (probability of an odd number of substitutions at a particular site)
    total_branch_length = treeanc.tree.total_branch_length()
    corrected_branch_length = np.sum([np.exp(-x.branch_length)*np.sinh(x.branch_length)
                                      for x in treeanc.tree.find_clades()])
    corrected_terminal_branch_length = np.sum([np.exp(-x.branch_length)*np.sinh(x.branch_length)
                                      for x in treeanc.tree.get_terminals()])
    expected_mutations = L*corrected_branch_length
    expected_terminal_mutations = L*corrected_terminal_branch_length

    # make histograms and sum mutations in different categories
    multiplicities = np.bincount([len(x) for x in mutations.values()])
    total_mutations = np.sum([len(x) for x in mutations.values()])

    multiplicities_terminal = np.bincount([len(x) for x in terminal_mutations.values()])
    terminal_mutation_count = np.sum([len(x) for x in terminal_mutations.values()])

    multiplicities_positions = np.bincount([len(x) for x in positions.values()])
    multiplicities_positions[0] = L - np.sum(multiplicities_positions)

    ###########################################################################
    ### Output the distribution of times particular mutations are observed
    ###########################################################################
    print("\nThe TOTAL tree length is %1.3e and %d mutations were observed."
          %(total_branch_length,total_mutations))
    print("Of these %d mutations,"%total_mutations
            +"".join(['\n\t - %d occur %d times'%(n,mi)
                      for mi,n in enumerate(multiplicities) if n]))
    # additional optional output this for terminal mutations only
    if params.detailed:
        print("\nThe TERMINAL branch length is %1.3e and %d mutations were observed."
              %(corrected_terminal_branch_length,terminal_mutation_count))
        print("Of these %d mutations,"%terminal_mutation_count
                +"".join(['\n\t - %d occur %d times'%(n,mi)
                          for mi,n in enumerate(multiplicities_terminal) if n]))


    ###########################################################################
    ### Output the distribution of times mutations at particular positions are observed
    ###########################################################################
    print("\nOf the %d positions in the genome,"%L
            +"".join(['\n\t - %d were hit %d times (expected %1.2f)'%(n,mi,L*poisson.pmf(mi,1.0*total_mutations/L))
                      for mi,n in enumerate(multiplicities_positions) if n]))


    # compare that distribution to a Poisson distribution with the same mean
    p = poisson.pmf(np.arange(10*multiplicities_positions.max()),1.0*total_mutations/L)
    print("\nlog-likelihood difference to Poisson distribution with same mean: %1.3e"%(
            - L*np.sum(p*np.log(p+1e-100))
            + np.sum(multiplicities_positions*np.log(p[:len(multiplicities_positions)]+1e-100))))


    ###########################################################################
    ### Output the mutations that are observed most often
    ###########################################################################
    if params.drms:
        print("\n\nThe ten most homoplasic mutations are:\n\tmut\tmultiplicity\tDRM details (gene drug AAmut)")
        mutations_sorted = sorted(mutations.items(), key=lambda x:len(x[1])-0.1*x[0][1]/L, reverse=True)
        for mut, val in mutations_sorted[:params.n]:
            if len(val)>1:
                print("\t%s%d%s\t%d\t%s"%(mut[0], mut[1], mut[2], len(val),
                    " ".join([drms[mut[1]]['gene'], drms[mut[1]]['drug'], drms[mut[1]]['alt_base'][mut[2]]]) if mut[1] in drms else ""))
            else:
                break
    else:
        print("\n\nThe ten most homoplasic mutations are:\n\tmut\tmultiplicity")
        mutations_sorted = sorted(mutations.items(), key=lambda x:len(x[1])-0.1*x[0][1]/L, reverse=True)
        for mut, val in mutations_sorted[:params.n]:
            if len(val)>1:
                print("\t%s%d%s\t%d"%(mut[0], mut[1], mut[2], len(val)))
            else:
                break

    # optional output specifically for mutations on terminal branches
    if params.detailed:
        if params.drms:
            print("\n\nThe ten most homoplasic mutation on terminal branches are:\n\tmut\tmultiplicity\tDRM details (gene drug AAmut)")
            terminal_mutations_sorted = sorted(terminal_mutations.items(), key=lambda x:len(x[1])-0.1*x[0][1]/L, reverse=True)
            for mut, val in terminal_mutations_sorted[:params.n]:
                if len(val)>1:
                    print("\t%s%d%s\t%d\t%s"%(mut[0], mut[1], mut[2], len(val),
                        " ".join([drms[mut[1]]['gene'], drms[mut[1]]['drug'], drms[mut[1]]['alt_base'][mut[2]]]) if mut[1] in drms else ""))
                else:
                    break
        else:
            print("\n\nThe ten most homoplasic mutation on terminal branches are:\n\tmut\tmultiplicity")
            terminal_mutations_sorted = sorted(terminal_mutations.items(), key=lambda x:len(x[1])-0.1*x[0][1]/L, reverse=True)
            for mut, val in terminal_mutations_sorted[:params.n]:
                if len(val)>1:
                    print("\t%s%d%s\t%d"%(mut[0], mut[1], mut[2], len(val)))
                else:
                    break

    ###########################################################################
    ### Output strains that have many homoplasic mutations
    ###########################################################################
    # TODO: add statistical criterion
    if params.detailed:
        if params.drms:
            print("\n\nTaxons that carry positions that mutated elsewhere in the tree:\n\ttaxon name\t#of homoplasic mutations\t# DRM")
            mutation_by_strain_sorted = sorted(mutation_by_strain.items(), key=lambda x:len(x[1]), reverse=True)
            for name, val in mutation_by_strain_sorted[:params.n]:
                if len(val):
                    print("\t%s\t%d\t%d"%(name, len(val),
                        len([mut for mut,l in val if mut[1] in drms])))
        else:
            print("\n\nTaxons that carry positions that mutated elsewhere in the tree:\n\ttaxon name\t#of homoplasic mutations")
            mutation_by_strain_sorted = sorted(mutation_by_strain.items(), key=lambda x:len(x[1]), reverse=True)
            for name, val in mutation_by_strain_sorted[:params.n]:
                if len(val):
                    print("\t%s\t%d"%(name, len(val)))


    return 0


def timetree(params):
    """
    implementeing treetime tree
    """
    if params.relax is None:
        relaxed_clock_params = None
    elif params.relax==[]:
        relaxed_clock_params=True
    elif len(params.relax)==2:
        relaxed_clock_params={'slack':params.relax[0], 'coupling':params.relax[1]}


    dates = utils.parse_dates(params.dates)
    if len(dates)==0:
        print("No valid dates -- exiting.")
        return 1

    if assure_tree(params, tmp_dir='timetree_tmp'):
        print("No tree -- exiting.")
        return 1

    outdir = get_outdir(params, '_treetime')

    gtr = create_gtr(params)
    infer_gtr = params.gtr=='infer'

    ###########################################################################
    ### READ IN VCF
    ###########################################################################
    #sets ref and fixed_pi to None if not VCF
    aln, ref, fixed_pi = read_if_vcf(params)
    is_vcf = True if ref is not None else False
    branch_length_mode = params.branch_length_mode
    #variable-site-only trees can have big branch lengths, the auto setting won't work.
    if is_vcf or (params.aln and params.sequence_length):
        if branch_length_mode == 'auto':
            branch_length_mode = 'joint'



    ###########################################################################
    ### SET-UP and RUN
    ###########################################################################
    if params.aln is None and params.sequence_length is None:
        print("one of arguments '--aln' and '--sequence-length' is required.", file=sys.stderr)
        return 1
    myTree = TreeTime(dates=dates, tree=params.tree, ref=ref,
                      aln=aln, gtr=gtr, seq_len=params.sequence_length,
                      verbose=params.verbose)
    myTree.tip_slack=params.tip_slack
    if not myTree.one_mutation:
        print("TreeTime setup failed, exiting")
        return 1

    # coalescent model options
    try:
        coalescent = float(params.coalescent)
        if coalescent<10*myTree.one_mutation:
            coalescent = None
    except:
        if params.coalescent in ['opt', 'const', 'skyline']:
            coalescent = params.coalescent
        else:
            print("unknown coalescent model specification, has to be either "
                  "a float, 'opt', 'const' or 'skyline' -- exiting")
            return 1

    # determine whether confidence intervals are to be computed and how the
    # uncertainty in the rate estimate should be treated
    calc_confidence = params.confidence
    if params.clock_std_dev:
        vary_rate = params.clock_std_dev if calc_confidence else False
    elif params.confidence and params.covariation:
        vary_rate = True
    elif params.confidence:
        print("\nOutside of covariation aware mode TreeTime cannot estimate confidence intervals "
                "without specified standard deviation of the clock rate. \nPlease specify '--clock-std-dev' "
                "or rerun with '--covariation'. \nWill proceed without confidence estimation")
        vary_rate = False
        calc_confidence = False
    else:
        vary_rate = False

    # RUN
    root = None if params.keep_root else params.reroot
    success = myTree.run(root=root, relaxed_clock=relaxed_clock_params,
               resolve_polytomies=(not params.keep_polytomies),
               Tc=coalescent, max_iter=params.max_iter,
               fixed_clock_rate=params.clock_rate,
               n_iqd=params.clock_filter,
               time_marginal="assign" if calc_confidence else False,
               vary_rate = vary_rate,
               branch_length_mode = branch_length_mode,
               fixed_pi=fixed_pi,
               use_covariation = params.covariation)
    if success==ttconf.ERROR: # if TreeTime.run failed, exit
        print("\nTreeTime run FAILED: please check above for errors and/or rerun with --verbose 4.\n")
        return 1

    ###########################################################################
    ### OUTPUT and saving of results
    ###########################################################################
    if infer_gtr:
        fname = outdir+'sequence_evolution_model.txt'
        with open(fname, 'w') as ofile:
            ofile.write(str(myTree.gtr)+'\n')
        print('\nInferred sequence evolution model (saved as %s):'%fname)
        print(myTree.gtr)

    fname = outdir+'molecular_clock.txt'
    with open(fname, 'w') as ofile:
        ofile.write(str(myTree.date2dist)+'\n')
    print('\nInferred sequence evolution model (saved as %s):'%fname)
    print(myTree.date2dist)

    basename = get_basename(params, outdir)
    if coalescent in ['skyline', 'opt', 'const']:
        print("Inferred coalescent model")
        if coalescent=='skyline':
            print_save_plot_skyline(myTree, plot=basename+'skyline.pdf', save=basename+'skyline.tsv', screen=True)
        else:
            Tc = myTree.merger_model.Tc.y[0]
            print(" --T_c: \t %1.2e \toptimized inverse merger rate in units of substitutions"%Tc)
            print(" --T_c: \t %1.2e \toptimized inverse merger rate in years"%(Tc/myTree.date2dist.clock_rate))
            print(" --N_e: \t %1.2e \tcorresponding 'effective population size' assuming 50 gen/year\n"%(Tc/myTree.date2dist.clock_rate*50))

    # plot
    import matplotlib.pyplot as plt
    from .treetime import plot_vs_years
    leaf_count = myTree.tree.count_terminals()
    label_func = lambda x: (x.name if x.is_terminal() and ((leaf_count<30
                                        and (not params.no_tip_labels))
                                      or params.tip_labels) else '')

    plot_vs_years(myTree, show_confidence=False, label_func=label_func,
                  confidence=0.9 if params.confidence else None)
    tree_fname = (outdir + params.plot_tree)
    plt.savefig(tree_fname)
    print("--- saved tree as \n\t %s\n"%tree_fname)

    plot_rtt(myTree, outdir + params.plot_rtt)
    if params.relax:
        fname = outdir+'substitution_rates.tsv'
        print("--- wrote branch specific rates to\n\t %s\n"%fname)
        with open(fname, 'w') as fh:
            fh.write("#node\tclock_length\tmutation_length\trate\tfold_change\n")
            for n in myTree.tree.find_clades(order="preorder"):
                if n==myTree.tree.root:
                    continue
                g = n.branch_length_interpolator.gamma
                fh.write("%s\t%1.3e\t%1.3e\t%1.3e\t%1.2f\n"%(n.name, n.clock_length, n.mutation_length, myTree.date2dist.clock_rate*g, g))

    export_sequences_and_tree(myTree, basename, is_vcf, params.zero_based,
                              timetree=True, confidence=calc_confidence)

    return 0


def ancestral_reconstruction(params):
    """
    implementing treetime ancestral
    """

    # set up
    if assure_tree(params, tmp_dir='ancestral_tmp'):
        return 1

    outdir = get_outdir(params, '_ancestral')
    basename = get_basename(params, outdir)

    gtr = create_gtr(params)

    ###########################################################################
    ### READ IN VCF
    ###########################################################################
    #sets ref and fixed_pi to None if not VCF
    aln, ref, fixed_pi = read_if_vcf(params)
    is_vcf = True if ref is not None else False

    treeanc = TreeAnc(params.tree, aln=aln, ref=ref, gtr=gtr, verbose=1,
                      fill_overhangs=not params.keep_overhangs)
    ndiff =treeanc.infer_ancestral_sequences('ml', infer_gtr=params.gtr=='infer',
                                             marginal=params.marginal, fixed_pi=fixed_pi)
    if ndiff==ttconf.ERROR: # if reconstruction failed, exit
        return 1

    ###########################################################################
    ### OUTPUT and saving of results
    ###########################################################################
    if params.gtr=='infer':
        fname = outdir+'/sequence_evolution_model.txt'
        with open(fname, 'w') as ofile:
            ofile.write(str(treeanc.gtr)+'\n')
        print('\nInferred sequence evolution model (saved as %s):'%fname)
        print(treeanc.gtr)

    export_sequences_and_tree(treeanc, basename, is_vcf, params.zero_based,
                              report_ambiguous=params.report_ambiguous)

    return 0

def reconstruct_discrete_traits(tree, traits, missing_data='?', pc=1.0, sampling_bias_correction=None,
                                weights=None, verbose=0):
    unique_states = sorted(set(traits.values()))
    nc = len(unique_states)
    if nc>180:
        print("mugration: can't have more than 180 states!", file=sys.stderr)
        return 1
    elif nc<2:
        print("mugration: only one or zero states found -- this doesn't make any sense", file=sys.stderr)
        return 1

    ###########################################################################
    ### make a single character alphabet that maps to discrete states
    ###########################################################################
    alphabet = [chr(65+i) for i,state in enumerate(unique_states)]
    missing_char = chr(65+nc)
    letter_to_state = {a:unique_states[i] for i,a in enumerate(alphabet)}
    letter_to_state[missing_char]=missing_data
    reverse_alphabet = {v:k for k,v in letter_to_state.items()}

    ###########################################################################
    ### construct gtr model
    ###########################################################################
    if type(weights)==str:
        tmp_weights = pd.read_csv(weights, sep='\t' if weights[-3:]=='tsv' else ',',
                             skipinitialspace=True)
        weights = {row[0]:row[1] for ri,row in tmp_weights.iterrows()}
        mean_weight = np.mean(list(weights.values()))
        weights = np.array([weights[c] if c in weights else mean_weight for c in unique_states], dtype=float)
        weights/=weights.sum()
    else:
        weights = None

    # set up dummy matrix
    W = np.ones((nc,nc), dtype=float)

    mugration_GTR = GTR.custom(pi = weights, W=W, alphabet = np.array(alphabet))
    mugration_GTR.profile_map[missing_char] = np.ones(nc)
    mugration_GTR.ambiguous=missing_char

    ###########################################################################
    ### set up treeanc
    ###########################################################################
    treeanc = TreeAnc(tree, gtr=mugration_GTR, verbose=verbose,
                      convert_upper=False, one_mutation=0.001)
    treeanc.use_mutation_length = False
    pseudo_seqs = [SeqRecord(id=n.name,name=n.name,
                   seq=Seq(reverse_alphabet[traits[n.name]]
                           if n.name in traits else missing_char))
                   for n in treeanc.tree.get_terminals()]
    treeanc.aln = MultipleSeqAlignment(pseudo_seqs)

    ndiff = treeanc.infer_ancestral_sequences(method='ml', infer_gtr=True,
            store_compressed=False, pc=pc, marginal=True, normalized_rate=False,
            fixed_pi=weights)

    if ndiff==ttconf.ERROR: # if reconstruction failed, exit
        return 1

    if sampling_bias_correction:
        treeanc.gtr.mu *= sampling_bias_correction
        treeanc.infer_ancestral_sequences(infer_gtr=False, store_compressed=False,
                                     marginal=True, normalized_rate=False)
    return treeanc, letter_to_state, reverse_alphabet


def mugration(params):
    """
    implementing treetime mugration
    """

    ###########################################################################
    ### Parse states
    ###########################################################################
    if os.path.isfile(params.states):
        states = pd.read_csv(params.states, sep='\t' if params.states[-3:]=='tsv' else ',',
                             skipinitialspace=True)
    else:
        print("file with states does not exist")
        return 1

    outdir = get_outdir(params, '_mugration')

    taxon_name = 'name' if 'name' in states.columns else states.columns[0]
    if params.attribute:
        if params.attribute in states.columns:
            attr = params.attribute
        else:
            print("The specified attribute was not found in the metadata file "+params.states, file=sys.stderr)
            print("Available columns are: "+", ".join(states.columns), file=sys.stderr)
            return 1
    else:
        attr = states.columns[1]
        print("Attribute for mugration inference was not specified. Using "+attr, file=sys.stderr)

    leaf_to_attr = {x[taxon_name]:x[attr] for xi, x in states.iterrows()
                    if x[attr]!=params.missing_data}

    mug, letter_to_state, reverse_alphabet = reconstruct_discrete_traits(params.tree, leaf_to_attr, missing_data=params.missing_data,
            pc=params.pc, sampling_bias_correction=params.sampling_bias_correction, verbose=params.verbose)

    unique_states = sorted(letter_to_state.values())
    ###########################################################################
    ### output
    ###########################################################################
    print("\nCompleted mugration model inference of attribute '%s' for"%attr,params.tree)

    basename = get_basename(params, outdir)
    gtr_name = basename + 'GTR.txt'
    with open(gtr_name, 'w') as ofile:
        ofile.write('Character to attribute mapping:\n')
        for state in unique_states:
            ofile.write('  %s: %s\n'%(reverse_alphabet[state], state))
        ofile.write('\n\n'+str(mug.gtr)+'\n')
        print("\nSaved inferred mugration model as:", gtr_name)

    terminal_count = 0
    for n in mug.tree.find_clades():
        if n.up is None:
            continue
        n.confidence=None
        # due to a bug in older versions of biopython that truncated filenames in nexus export
        # we truncate them by hand and make them unique.
        if n.is_terminal() and len(n.name)>40 and bioversion<"1.69":
            n.name = n.name[:35]+'_%03d'%terminal_count
            terminal_count+=1
        n.comment= '&%s="'%attr + letter_to_state[n.sequence[0]] +'"'

    if params.confidence:
        conf_name = basename+'confidence.csv'
        with open(conf_name, 'w') as ofile:
            ofile.write('#name, '+', '.join(unique_states)+'\n')
            for n in mug.tree.find_clades():
                ofile.write(n.name + ', '+', '.join([str(x) for x in n.marginal_profile[0]])+'\n')
        print("Saved table with ancestral state confidences as:", conf_name)

    # write tree to file
    outtree_name = basename+'annotated_tree.nexus'
    Phylo.write(mug.tree, outtree_name, 'nexus')
    print("Saved annotated tree as:", outtree_name)

    return 0


def estimate_clock_model(params):
    """
    implementing treetime clock
    """

    if assure_tree(params, tmp_dir='clock_model_tmp'):
        return 1
    dates = utils.parse_dates(params.dates)
    if len(dates)==0:
        return 1

    outdir = get_outdir(params, '_clock')

    ###########################################################################
    ### READ IN VCF
    ###########################################################################
    #sets ref and fixed_pi to None if not VCF
    aln, ref, fixed_pi = read_if_vcf(params)
    is_vcf = True if ref is not None else False

    ###########################################################################
    ### ESTIMATE ROOT (if requested) AND DETERMINE TEMPORAL SIGNAL
    ###########################################################################
    if params.aln is None and params.sequence_length is None:
        print("one of arguments '--aln' and '--sequence-length' is required.", file=sys.stderr)
        return 1

    basename = get_basename(params, outdir)
    myTree = TreeTime(dates=dates, tree=params.tree, aln=aln, gtr='JC69',
                      verbose=params.verbose, seq_len=params.sequence_length,
                      ref=ref)
    myTree.tip_slack=params.tip_slack
    if myTree.tree is None:
        print("ERROR: tree loading failed. exiting...")
        return 1

    if params.clock_filter:
        n_bad = [n.name for n in myTree.tree.get_terminals() if n.bad_branch]
        myTree.clock_filter(n_iqd=params.clock_filter, reroot=params.reroot or 'least-squares')
        n_bad_after = [n.name for n in myTree.tree.get_terminals() if n.bad_branch]
        if len(n_bad_after)>len(n_bad):
            print("The following leaves don't follow a loose clock and "
                  "will be ignored in rate estimation:\n\t"
                  +"\n\t".join(set(n_bad_after).difference(n_bad)))

    if not params.keep_root:
        # reroot to optimal root, this assigns clock_model to myTree
        if params.covariation: # this requires branch length estimates
            myTree.run(root="least-squares", max_iter=0,
                       use_covariation=params.covariation)

        res = myTree.reroot(params.reroot,
                      force_positive=not params.allow_negative_rate)
        myTree.get_clock_model(covariation=params.covariation)

        if res==ttconf.ERROR:
            print("ERROR: unknown root or rooting mechanism!\n"
                  "\tvalid choices are 'least-squares', 'ML', and 'ML-rough'")
            return 1
    else:
        myTree.get_clock_model(covariation=params.covariation)

    d2d = utils.DateConversion.from_regression(myTree.clock_model)
    print('\n',d2d)
    print('The R^2 value indicates the fraction of variation in'
          '\nroot-to-tip distance explained by the sampling times.'
          '\nHigher values corresponds more clock-like behavior (max 1.0).')

    print('\nThe rate is the slope of the best fit of the date to'
          '\nthe root-to-tip distance and provides an estimate of'
          '\nthe substitution rate. The rate needs to be positive!'
          '\nNegative rates suggest an inappropriate root.\n')

    print('\nThe estimated rate and tree correspond to a root date:')
    if params.covariation:
        reg = myTree.clock_model
        dp = np.array([reg['intercept']/reg['slope']**2,-1./reg['slope']])
        droot = np.sqrt(reg['cov'][:2,:2].dot(dp).dot(dp))
        print('\n--- root-date:\t %3.2f +/- %1.2f (one std-dev)\n\n'%(-d2d.intercept/d2d.clock_rate, droot))
    else:
        print('\n--- root-date:\t %3.2f\n\n'%(-d2d.intercept/d2d.clock_rate))

    if not params.keep_root:
        # write rerooted tree to file
        outtree_name = basename+'rerooted.newick'
        Phylo.write(myTree.tree, outtree_name, 'newick')
        print("--- re-rooted tree written to \n\t%s\n"%outtree_name)

    table_fname = basename+'rtt.csv'
    with open(table_fname, 'w') as ofile:
        ofile.write("#name, date, root-to-tip distance\n")
        ofile.write("#Dates of nodes that didn't have a specified date are inferred from the root-to-tip regression.\n")
        for n in myTree.tree.get_terminals():
            if hasattr(n, "raw_date_constraint") and (n.raw_date_constraint is not None):
                if np.isscalar(n.raw_date_constraint):
                    tmp_str = str(n.raw_date_constraint)
                elif len(n.raw_date_constraint):
                    tmp_str = str(n.raw_date_constraint[0])+'-'+str(n.raw_date_constraint[1])
                else:
                    tmp_str = ''
                ofile.write("%s, %s, %f\n"%(n.name, tmp_str, n.dist2root))
            else:
                ofile.write("%s, %f, %f\n"%(n.name, d2d.numdate_from_dist2root(n.dist2root), n.dist2root))
        for n in myTree.tree.get_nonterminals(order='preorder'):
            ofile.write("%s, %f, %f\n"%(n.name, d2d.numdate_from_dist2root(n.dist2root), n.dist2root))
        print("--- wrote dates and root-to-tip distances to \n\t%s\n"%table_fname)


    ###########################################################################
    ### PLOT AND SAVE RESULT
    ###########################################################################
    plot_rtt(myTree, outdir+params.plot_rtt)
    return 0
