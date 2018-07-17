from __future__ import print_function, division, absolute_import
import os, shutil
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
        gtr = GTR.standard('jc')
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
            gtr = GTR.standard('jc')
            infer_gtr = False
    return gtr

def parse_dates(params):
    """
    parse dates from the arguments and return a dictionary mapping
    taxon names to numerical dates.
    """
    dates = {}
    if not os.path.isfile(params.dates):
        print("\n ERROR: file %s does not exist, exiting..."%params.dates)
        return dates
    with open(params.dates) as date_file:
        failed_dates = 0
        for line in date_file:
            try:
                name, date = line.strip().split(',')[:2]
                dates[name] = float(date)
            except:
                failed_dates+=1

        if len(dates)<failed_dates:
            print("\n\nDATE PARSING FAILED, ABORTING...")

    return dates

def read_if_vcf(params):
    """
    Checks if input is VCF and reads in appropriately if it is
    """
    ref = None
    aln = params.aln
    fixed_pi = None
    if any([params.aln.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
        if not params.vcf_reference:
            print("ERROR: a reference Fasta is required with VCF-format alignments")
            return -1
        compress_seq = read_vcf(params.aln, params.vcf_reference)
        sequences = compress_seq['sequences']
        ref = compress_seq['reference']
        aln = sequences

        if not hasattr(params, 'gtr') or params.gtr=="infer": #if not specified, set it:
            fixed_pi = [ref.count(base)/len(ref) for base in ['A','C','G','T','-']]
            if fixed_pi[-1] == 0:
                fixed_pi[-1] = 0.05
                fixed_pi = [v-0.01 for v in fixed_pi]

    return aln, ref, fixed_pi

def scan_homoplasies(params):
    """
    the function implementing treetime homoplasies
    """
    if assure_tree(params, tmp_dir='homoplasy_tmp'):
        return 1

    gtr = create_gtr(params)

    ###########################################################################
    ### ANCESTRAL RECONSTRUCTION
    ###########################################################################
    treeanc = TreeAnc(params.tree, aln=params.aln, gtr=gtr, verbose=1,
                      fill_overhangs=True)
    if treeanc.aln is None: # if alignment didn't load, exit
        return 1

    # FIXME: resetting one_mutation is no longer possible
    L = treeanc.aln.get_alignment_length() + params.const
    treeanc.one_mutation = 1.0/L
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
                                      marginal=False)
    print("...done.")
    if ndiff==ttconf.ERROR: # if reconstruction failed, exit
        return 1
    else:
        print("...done.")

    ###########################################################################
    ### analysis of reconstruction
    ###########################################################################
    from collections import defaultdict
    from scipy.stats import poisson
    offset = 0 if params.zero_based else 1

    # construct dictionaries gathering mutations and positions
    mutations = defaultdict(list)
    positions = defaultdict(list)
    terminal_mutations = defaultdict(list)
    for n in treeanc.tree.find_clades():
        if n.up is None:
            continue

        if len(n.mutations):
            for (a,pos, d) in n.mutations:
                if '-' not in [a,d]:
                    mutations[(a,pos+offset,d)].append(n)
                    positions[pos+offset].append(n)
            if n.is_terminal():
                for (a,pos, d) in n.mutations:
                    if '-' not in [a,d]:
                        terminal_mutations[(a,pos+offset,d)].append(n)

    # gather homoplasic mutations by strain
    mutation_by_strain = defaultdict(list)
    for n in treeanc.tree.get_terminals():
        for a,pos,d in n.mutations:
            if pos in positions and len(positions[pos])>1:
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
    print("\nThe TOTAL tree length is %1.3e, expecting %1.1f mutations vs an observed %d"
          %(total_branch_length,expected_mutations,total_mutations))
    print("Of these %d mutations,"%total_mutations
            +"".join(['\n\t - %d occur %d times'%(n,mi)
                      for mi,n in enumerate(multiplicities) if n]))
    # additional optional output this for terminal mutations only
    if params.detailed:
        print("\nThe TERMINAL branch length is %1.3e, expecting %1.1f mutations vs an observed %d"
              %(corrected_terminal_branch_length,expected_terminal_mutations,terminal_mutation_count))
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
    print("\n\nThe ten most homoplasic mutations are:\n\tmut\tmultiplicity")
    mutations_sorted = sorted(mutations.items(), key=lambda x:len(x[1])-0.1*x[0][1]/L, reverse=True)
    for mut, val in mutations_sorted[:params.n]:
        if len(val)>1:
            print("\t%s%d%s\t%d"%(mut[0], mut[1], mut[2], len(val)))
        else:
            break

    # optional output specifically for mutations on terminal branches
    if params.detailed:
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
    if params.relax==[]:
        params.relax=True

    dates = parse_dates(params)
    if len(dates)==0:
        return 1

    if assure_tree(params, tmp_dir='timetree_tmp'):
        return 1

    gtr = create_gtr(params)
    infer_gtr = params.gtr=='infer'

    ###########################################################################
    # PARSING OPTIONS
    ###########################################################################
    try:
        Tc = float(params.Tc)
        if Tc<1e-5:
            Tc = None
    except:
        if params.Tc in ['opt', 'skyline']:
            Tc = params.Tc
        else:
            Tc = None

    ###########################################################################
    ### READ IN VCF
    ###########################################################################
    #sets ref and fixed_pi to None if not VCF
    aln, ref, fixed_pi = read_if_vcf(params)
    is_vcf = True if ref is not None else False

    ###########################################################################
    ### SET-UP and RUN
    ###########################################################################
    myTree = TreeTime(dates=dates, tree=params.tree, ref=ref,
                       aln=aln, gtr=gtr, verbose=params.verbose)
    root = None if params.keep_root else params.reroot
    success = myTree.run(root=root, relaxed_clock=params.relax,
               resolve_polytomies=(not params.keep_polytomies),
               Tc=Tc, max_iter=params.max_iter,
               fixed_clock_rate=params.clock_rate,
               branch_length_mode = params.branch_length_mode,
               fixed_pi=fixed_pi)
    if success==ttconf.ERROR: # if TreeTime.run failed, exit
        return 1

    ###########################################################################
    ### OUTPUT and saving of results
    ###########################################################################
    if infer_gtr:
        print('\nInferred GTR model:')
        print(myTree.gtr)

    print(myTree.date2dist)

    if Tc=='skyline':
        skyline = myTree.merger_model.skyline_inferred(gen=50)
        print("inferred skyline assuming 50 generations per year:")
        for (x,y) in zip(skyline.x, skyline.y):
            print("%1.3f\t%1.3f"%(x,y))


    base_name = '.'.join(params.aln.split('/')[-1].split('.')[:-1])
    # plot
    if params.plot:
        import matplotlib.pyplot as plt
        from .treetime import plot_vs_years
        plt.ion()
        leaf_count = myTree.tree.count_terminals()
        label_func = lambda x: x.name[:20] if (leaf_count<30 & x.is_terminal()) else ''
        # branch_label_func = lambda x: (','.join([a+str(pos)+d for a,pos, d in x.mutations[:10]])
        #                                +('...' if  len(x.mutations)>10 else '')) if leaf_count<30 else ''
        plot_vs_years(myTree, show_confidence=False, label_func = label_func) #, branch_labels=branch_label_func)
        plt.savefig(base_name+'_tree.pdf')
        print("--- saved tree as pdf in \n\t %s\n"%(base_name+'_tree.pdf'))
    else:
        # convert branch length to years (this is implicit in the above plot)
        myTree.branch_length_to_years()

    # decorate tree with inferred mutations
    if is_vcf:
        outaln_name = base_name+'_ancestral.vcf'
        write_vcf(myTree.get_tree_dict(keep_var_ambigs=True), outaln_name)
    else:
        outaln_name = base_name+'_ancestral.fasta'
        AlignIO.write(myTree.get_reconstructed_alignment(), outaln_name, 'fasta')
    print("--- alignment including ancestral nodes saved as  \n\t %s\n"%outaln_name)

    terminal_count = 0
    for n in myTree.tree.find_clades():
        if n.up is None:
            continue
        n.confidence=None
        # due to a bug in older versions of biopython that truncated filenames in nexus export
        # we truncate them by hand and make them unique.
        if n.is_terminal() and len(n.name)>40 and bioversion<"1.69":
            n.name = n.name[:35]+'_%03d'%terminal_count
            terminal_count+=1
        if len(n.mutations):
            n.comment= '&mutations="' + '_'.join([a+str(pos)+d for (a,pos, d) in n.mutations])+'"'

    # write tree to file
    outtree_name = '.'.join(params.tree.split('/')[-1].split('.')[:-1])+'_timetree.nexus'
    Phylo.write(myTree.tree, outtree_name, 'nexus')

    print("--- tree saved in nexus format as  \n\t %s\n"%outtree_name)
    return 0


def ancestral_reconstruction(params):
    """
    implementing treetime ancestral
    """

    # set up
    if assure_tree(params, tmp_dir='ancestral_tmp'):
        return 1

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
                                             marginal=params.marginal)
    if ndiff==ttconf.ERROR: # if reconstruction failed, exit
        return 1

    ###########################################################################
    ### OUTPUT and saving of results
    ###########################################################################
    if params.gtr=="infer":
        print('\nInferred GTR model:')
        print(treeanc.gtr)

    if is_vcf:
        outaln_name = '.'.join(params.aln.split('/')[-1].split('.')[:-1])+'_ancestral.vcf'
        write_vcf(treeanc.get_tree_dict(keep_var_ambigs=True), outaln_name)
    else:
        outaln_name = '.'.join(params.aln.split('/')[-1].split('.')[:-1])+'_ancestral.fasta'
        AlignIO.write(treeanc.get_reconstructed_alignment(), outaln_name, 'fasta')
    print("--- alignment including ancestral nodes saved as  \n\t %s\n"%outaln_name)

    # decorate tree with inferred mutations
    terminal_count = 0
    offset = 0 if params.zero_based else 1
    for n in treeanc.tree.find_clades():
        if n.up is None:
            continue
        n.confidence=None
        # due to a bug in older versions of biopython that truncated filenames in nexus export
        # we truncate them by hand and make them unique.
        if n.is_terminal() and len(n.name)>40 and bioversion<"1.69":
            n.name = n.name[:35]+'_%03d'%terminal_count
            terminal_count+=1
        if len(n.mutations):
            n.comment= '&mutations="' + '_'.join([a+str(pos + offset)+d for (a,pos, d) in n.mutations])+'"'

    # write tree to file
    outtree_name = '.'.join(params.tree.split('/')[-1].split('.')[:-1])+'_mutation.nexus'
    Phylo.write(treeanc.tree, outtree_name, 'nexus')
    print("--- tree saved in nexus format as  \n\t %s\n"%outtree_name)

    return 0

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

    taxon_name = 'name' if 'name' in states.columns else states.columns[0]
    if params.attribute and params.attribute in states.columns:
        attr = params.attribute
    else:
        attr = states.columns[1]

    leaf_to_attr = {x[taxon_name]:x[attr] for xi, x in states.iterrows()
                    if x[attr]!=params.missing_data}
    unique_states = sorted(set(leaf_to_attr.values()))
    nc = len(unique_states)
    if nc>180:
        print("mugration: can't have more than 180 states!")
        exit(1)
    elif nc<2:
        print("mugration: only one or zero states found -- this doesn't make any sense")
        exit(1)

    ###########################################################################
    ### make a single character alphabet that maps to discrete states
    ###########################################################################
    alphabet = [chr(65+i) for i,state in enumerate(unique_states)]
    missing_char = chr(65+nc)
    letter_to_state = {a:unique_states[i] for i,a in enumerate(alphabet)}
    letter_to_state[missing_char]=params.missing_data
    reverse_alphabet = {v:k for k,v in letter_to_state.items()}

    ###########################################################################
    ### construct gtr model
    ###########################################################################
    if params.weights:
        params.infer_gtr = True
        tmp_weights = pd.read_csv(params.weights, sep='\t' if params.states[-3:]=='tsv' else ',',
                             skipinitialspace=True)
        weights = {row[0]:row[1] for ri,row in tmp_weights.iterrows()}
        mean_weight = np.mean(list(weights.values()))
        weights = np.array([weights[c] if c in weights else mean_weight for c in unique_states], dtype=float)
        weights/=weights.sum()
    else:
        weights = np.ones(nc, dtype=float)/nc

    # set up dummy matrix
    W = np.ones((nc,nc), dtype=float)

    mugration_GTR = GTR.custom(pi = weights, W=W, alphabet = np.array(alphabet))
    mugration_GTR.profile_map[missing_char] = np.ones(nc)
    mugration_GTR.ambiguous=missing_char

    ###########################################################################
    ### set up treeanc
    ###########################################################################
    treeanc = TreeAnc(params.tree, gtr=mugration_GTR, verbose=params.verbose,
                      convert_upper=False, one_mutation=0.001)
    pseudo_seqs = [SeqRecord(id=n.name,name=n.name,
                   seq=Seq(reverse_alphabet[leaf_to_attr[n.name]]
                           if n.name in leaf_to_attr else missing_char))
                   for n in treeanc.tree.get_terminals()]
    treeanc.aln = MultipleSeqAlignment(pseudo_seqs)

    ndiff = treeanc.infer_ancestral_sequences(method='ml', infer_gtr=True,
            store_compressed=False, pc=params.pc, marginal=True, normalized_rate=False,
            fixed_pi=weights if params.weights else None)
    if ndiff==ttconf.ERROR: # if reconstruction failed, exit
        return 1


    ###########################################################################
    ### output
    ###########################################################################
    print("\nCompleted mugration model inference of attribute '%s' for"%attr,params.tree)

    bname = './'+os.path.basename(params.tree)
    gtr_name = bname + '.GTR.txt'
    with open(gtr_name, 'w') as ofile:
        ofile.write('Character to attribute mapping:\n')
        for state in unique_states:
            ofile.write('  %s: %s\n'%(reverse_alphabet[state], state))
        ofile.write('\n\n'+str(treeanc.gtr)+'\n')
        print("\nSaved inferred mugration model as:", gtr_name)

    terminal_count = 0
    for n in treeanc.tree.find_clades():
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
        conf_name = bname+'.confidence.csv'
        with open(conf_name, 'w') as ofile:
            ofile.write('#name, '+', '.join(unique_states)+'\n')
            for n in treeanc.tree.find_clades():
                ofile.write(n.name + ', '+', '.join([str(x) for x in n.marginal_profile[0]])+'\n')
        print("Saved table with ancestral state confidences as:", conf_name)

    # write tree to file
    outtree_name = bname+'.mugration.nexus'
    Phylo.write(treeanc.tree, outtree_name, 'nexus')
    print("Saved annotated tree as:",outtree_name)

    return 0

def estimate_clock_model(params):
    """
    implementing treetime clock
    """

    if assure_tree(params, tmp_dir='clock_model_tmp'):
        return 1
    dates = parse_dates(params)
    if len(dates)==0:
        return 1

    ###########################################################################
    ### READ IN VCF
    ###########################################################################
    #sets ref and fixed_pi to None if not VCF
    aln, ref, fixed_pi = read_if_vcf(params)
    is_vcf = True if ref is not None else False

    ###########################################################################
    ### ESTIMATE ROOT (if requested) AND DETERMINE TEMPORAL SIGNAL
    ###########################################################################
    base_name = '.'.join(params.tree.split('/')[-1].split('.')[:-1])
    myTree = TreeTime(dates=dates, tree=params.tree, aln=aln, gtr='JC69',
                      verbose=params.verbose, seq_len=params.sequence_length,
                      ref=ref)
    if myTree.tree is None:
        print("ERROR: tree loading failed. exiting...")
        return 1

    if not params.keep_root:
        # reroot to optimal root, this assigns clock_model to myTree
        myTree.reroot(params.reroot, force_positive=not params.allow_negative_rate)
    else:
        Treg = myTree.setup_TreeRegression(covariation=True)
        myTree.clock_model = Treg.regression()

    d2d = utils.DateConversion.from_regression(myTree.clock_model)
    print('\n',d2d)
    print('The R^2 value indicates the fraction of variation in'
          '\nroot-to-tip distance explained by the sampling times.'
          '\nHigher values corresponds more clock-like behavior (max 1.0).')

    print('\nThe rate is the slope of the best fit of the date to'
          '\nthe root-to-tip distance and provides an estimate of'
          '\nthe substitution rate. The rate needs to be positive!'
          '\nNegative rates suggest an inappropriate root.\n\n')

    print('\nThe estimated rate and tree correspond to a root date:\n')
    print('\n--root-date:\t %3.2f\n\n'%(-d2d.intercept/d2d.clock_rate))

    if not params.keep_root:
        # write rerooted tree to file
        outtree_name = base_name+'_rerooted.newick'
        Phylo.write(myTree.tree, outtree_name, 'newick')
        print("--- re-rooted tree written to \n\t %s\n"%outtree_name)

    table_fname = base_name+'_rtt.csv'
    with open(table_fname, 'w') as ofile:
        ofile.write("#name, date, root-to-tip distance\n")
        ofile.write("#Dates of nodes that didn't have a specified date are inferred from the root-to-tip regression.\n")
        for n in myTree.tree.get_terminals():
            if hasattr(n, "numdate_given"):
                ofile.write("%s, %f, %f\n"%(n.name, n.numdate_given, n.dist2root))
            else:
                ofile.write("%s, %f, %f\n"%(n.name, d2d.numdate_from_dist2root(n.dist2root), n.dist2root))
        for n in myTree.tree.get_nonterminals(order='preorder'):
            ofile.write("%s, %f, %f\n"%(n.name, d2d.numdate_from_dist2root(n.dist2root), n.dist2root))
        print("--- wrote dates and root-to-tip distances to \n\t %s\n"%table_fname)


    ###########################################################################
    ### PLOT AND SAVE RESULT
    ###########################################################################
    if params.plot:
        import matplotlib.pyplot as plt
        myTree.plot_root_to_tip()
        if params.output:
            fname = params.output
        else:
            fname = base_name+'_root_to_tip_regression.pdf'
        plt.savefig(fname)
        print("--- root-to-tip plot saved to  \n\t"+fname)

    return 0
