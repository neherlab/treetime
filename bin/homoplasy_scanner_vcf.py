#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
from treetime import TreeAnc, GTR
from Bio import Phylo, AlignIO
from Bio import __version__ as bioversion
import os,shutil
from treetime.vcf_utils import read_vcf

def read_in_DRMs(drm_file):
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


if __name__=="__main__":
    ###########################################################################
    ### parameter parsing
    ###########################################################################
    import argparse
    parser = argparse.ArgumentParser(
            description='Reconstructs ancestral sequences and maps mutations to the tree.'
                        ' The tree is then scanned for homoplasies. An excess number of homoplasies'
                        ' might suggest contamination, recombination, culture adaptation or similar. ')
    parser.add_argument('--aln', required = True, type = str,  help ="fasta/vcf file with input nucleotide sequences")
    parser.add_argument('--ref', required=False, type=str, help ="for VCF files, the reference fasta for the VCF")
    parser.add_argument('--tree', type = str,  help ="newick file with tree (optional if tree builders installed)")
    parser.add_argument('--detailed', required = False, action="store_true",  help ="generate a more detailed report")
    parser.add_argument('--gtr', required=False, type = str, default='infer', help="GTR model to use. "
        " Type 'infer' to infer the model from the data. Or, specify the model type. "
        " If the specified model requires additional options, use '--gtr_args' to specify those")

    parser.add_argument('--gtr_params', type=str, nargs='+', help="GTR parameters for the model "
        "specified by the --gtr argument. The parameters should be feed as 'key=value' list of parameters. "
        "Example: '--gtr K80 --gtr_params kappa=0.2 pis=0.25,0.25,0.25,0.25'. See the exact definitions of "
        " the parameters in the GTR creation methods in treetime/nuc_models.py. Only nucleotide models supported at present")

    parser.add_argument('--prot', default = False, action="store_true", help ="protein alignment")
    parser.add_argument('--zero_based', default = False, action='store_true', help='zero based SNP indexing')
    parser.add_argument('-n', default = 10, type=int, help='number of mutations/nodes that are printed to screen')
    parser.add_argument('--verbose', default = 1, type=int, help='verbosity of output 0-6')
    parser.add_argument('--drm', type=str, help="file with drug resistance mutation information")
    parser.add_argument('--vcf_ann', type=str, help="VCF file with 'ANN' field in INFO column - can be same as input file")
    params = parser.parse_args()


    ###########################################################################
    ### CHECK FOR TREE, build if not in place
    ### Doubt this will work for VCF files - unchecked.
    ###########################################################################
    if params.tree is None:
        from treetime.utils import tree_inference
        params.tree = os.path.basename(params.aln)+'.nwk'
        print("No tree given: inferring tree")
        tmp_dir = 'homoplasy_scanner_tmp_files'
        tree_inference(params.aln, params.tree, tmp_dir = tmp_dir)
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)


    ###########################################################################
    ### GTR SET-UP
    ###########################################################################
    model = params.gtr
    gtr_params = params.gtr_params
    if model == 'infer':
        gtr = GTR.standard('jc')
        infer_gtr = True
    else:
        infer_gtr = False
        try:
            kwargs = {}
            if gtr_params is not None:
                for param in gtr_params:
                    keyval = param.split('=')
                    if len(keyval)!=2: continue
                    if keyval[0] in ['pis', 'pi', 'Pi', 'Pis']:
                        keyval[1] = map(float, keyval[1].split(','))
                    elif keyval[0] not in ['alphabet']:
                        keyval[1] = float(keyval[1])
                    kwargs[keyval[0]] = keyval[1]
            else:
                print ("GTR params are not specified. Creating GTR model with default parameters")

            gtr = GTR.standard(model, **kwargs)
        except:
            print ("Could not create GTR model from input arguments. Using default (Jukes-Cantor 1969)")
            gtr = GTR.standard('jc')


    ###########################################################################
    ### ANCESTRAL RECONSTRUCTION
    ###########################################################################
    aln = params.aln
    ref = params.ref
    if params.ref:
        compress_seq = read_vcf(params.aln, params.ref)
        aln = compress_seq['sequences']
        ref = compress_seq['reference']

    treeanc = TreeAnc(params.tree, aln=aln, ref=ref, gtr=gtr, #verbose=1,
                      fill_overhangs=True)

    L = len(ref)
    N_seq = len(treeanc.aln)
    N_tree = treeanc.tree.count_terminals()

    print("read alignment from file %s with %d sequences of length %d"%(params.aln,N_seq,L))
    print("read tree from file %s with %d leaves"%(params.tree,N_tree))
    print("\ninferring ancestral sequences...")

    pi = None
    if params.ref: #if VCF, we need to fix pi
        #otherwise because of sequence length, mutation rate TO gap is overestimated
        fixed_pi = [ref.count(base)/len(ref) for base in ['A','C','G','T','-']]
        if fixed_pi[-1] == 0: #if no gaps in ref, set ~4% gaps.. adjust other bases accordingly
            fixed_pi[-1] = 0.05
            fixed_pi = [v-0.01 for v in fixed_pi]
    treeanc.infer_ancestral_sequences('ml', infer_gtr=infer_gtr, fixed_pi=fixed_pi, marginal=False)
    if params.ref:
        treeanc.recover_var_ambigs() #put Ns back on tips!
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

        #now doesn't track N's either.
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
                if 'N' not in [a,d]:
                    mutation_by_strain[n.name].append([(a,pos+offset,d), len(positions[pos+offset])])


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
    #Read in DRM information if supplied
    if params.drm:
        DRM_info = read_in_DRMs(params.drm)
        drms = DRM_info['DRMs']


    #read in vcf with annotations if supplied
    import vcf
    vcf_reader = vcf.Reader(filename=params.vcf_ann, compressed=True)
    byPos = {}
    for record in vcf_reader:
        byPos[record.POS] = record#.INFO['ANN']


    ###########################################################################
    ### Output the mutations that are observed most often
    ###########################################################################
    # print("\n\nThe ten most homoplasic mutations are:\n\tmut\tmultiplicity\tDRM details")
    # mutations_sorted = sorted(mutations.items(), key=lambda x:len(x[1])-0.1*x[0][1]/L, reverse=True)
    # for mut, val in mutations_sorted[:params.n]:
        # if len(val)>1:
            # print("\t%s%d%s\t%d\t%s"%(mut[0], mut[1], mut[2], len(val),
                # " ".join([drms[mut[1]]['gene'], drms[mut[1]]['drug'], drms[mut[1]]['alt_base'][mut[2]]]) if mut[1] in drms else ""))
        # else:
            # break

    #10 most homoplasic mutations, but with DRM and ANN information!
    mutations_sorted = sorted(mutations.items(), key=lambda x:len(x[1])-0.1*x[0][1]/L, reverse=True)
    tmp = []
    i=0
    annots = {}
    for mut, val in mutations_sorted[:params.n]:
        if len(val)>1:
            deet = [mut[0]+str(mut[1])+mut[2], byPos[mut[1]].REF, len(val)]
            if mut[1] in drms:
                deet.extend((drms[mut[1]]['gene'], drms[mut[1]]['drug'], drms[mut[1]]['alt_base'][mut[2]]))
                deet.append("")
            else:
                deet.extend(("","",""))
                if 'ANN' in byPos[mut[1]].INFO:
                    i+=1
                    deet.append(str(i))
                    annots[i] = byPos[mut[1]].INFO['ANN'][0]
                else:
                    deet.append("")
            tmp.append(deet)
        else:
            break

    print("\n\nThe {} most homoplasic mutations are:".format(params.n))
    firstTen = np.array(tmp)
    colN = ["Mutation", "Ref", "Multiplicity", "DRM-gene", "DRM-drug", "DRM-AA", "Annotation"]
    import pandas as pd
    print(pd.DataFrame(data=firstTen, columns=colN))
    print("")
    for k,v in annots.iteritems():
        print("%d:\t%s\n"%(k,v))



    # optional output specifically for mutations on terminal branches
    # if params.detailed:
        # print("\n\nThe ten most homoplasic mutations on terminal branches are:\n\tmut\tmultiplicity\tDRM details")
        # terminal_mutations_sorted = sorted(terminal_mutations.items(), key=lambda x:len(x[1])-0.1*x[0][1]/L, reverse=True)
        # for mut, val in terminal_mutations_sorted[:params.n]:
            # if len(val)>1:
                # print("\t%s%d%s\t%d\t%s"%(mut[0], mut[1], mut[2], len(val),
                    # " ".join([drms[mut[1]]['gene'], drms[mut[1]]['drug'], drms[mut[1]]['alt_base'][mut[2]]]) if mut[1] in drms else ""))
            # else:
                # break

    if params.detailed:
        terminal_mutations_sorted = sorted(terminal_mutations.items(), key=lambda x:len(x[1])-0.1*x[0][1]/L, reverse=True)
        tmp = []
        i=0
        annots = {}
        for mut, val in terminal_mutations_sorted[:params.n]:
            if len(val)>1:
                deet = [mut[0]+str(mut[1])+mut[2], byPos[mut[1]].REF, len(val)]
                if deet[0] in firstTen:
                    deet.append("Y")
                else:
                    deet.append("N")
                if mut[1] in drms:
                    deet.extend((drms[mut[1]]['gene'], drms[mut[1]]['drug'], drms[mut[1]]['alt_base'][mut[2]]))
                    deet.append("")
                else:
                    deet.extend(("","",""))
                    if 'ANN' in byPos[mut[1]].INFO:
                        i+=1
                        deet.append(str(i))
                        annots[i] = byPos[mut[1]].INFO['ANN'][0]
                    else:
                        deet.append("")
                tmp.append(deet)
            else:
                break

        print("\n\nThe {} most homoplasic mutations on terminal branches are:".format(params.n))
        secondTen = np.array(tmp)
        colN = ["Mutation", "Ref", "Multiplicity", "Seen-Above?", "DRM-gene", "DRM-drug", "DRM-AA", "Annotation"]
        import pandas as pd
        print(pd.DataFrame(data=secondTen, columns=colN))
        print("")
        for k,v in annots.iteritems():
            print("%d:\t%s\n"%(k,v))


    ###########################################################################
    ### Output strains that have many homoplasic mutations
    ###########################################################################
    # TODO: add statistical criterion
    # if params.detailed:
        # print("\n\nTaxons that carry positions that mutated elsewhere in the tree:\n\ttaxon_name\t#of homoplasic mutations \t# DRM")
        # mutation_by_strain_sorted = sorted(mutation_by_strain.items(), key=lambda x:len(x[1]), reverse=True)
        # for name, val in mutation_by_strain_sorted[:params.n]:
            # if len(val):
                # print("\t%s\t%d\t%d"%(name, len(val),
                        # len([mut for mut,l in val if mut[1] in drms])))

    if params.detailed:
        mutation_by_strain_sorted = sorted(mutation_by_strain.items(), key=lambda x:len(x[1]), reverse=True)
        tmp = []
        i=0
        for name, val in mutation_by_strain_sorted[:params.n]:
            if len(val):
                deet = [name, len(val), len([mut for mut,l in val if mut[1] in drms])]
                muts = [mut[0]+str(mut[1])+mut[2] for mut,count in val ]
                deet.append(",".join(muts))
                tmp.append(deet)

        print("\n\nTaxons that carry positions that mutated elsewhere in the tree:".format(params.n))
        secondTen = np.array(tmp)
        colN = ["TaxonName", "#_of_homoplasic_mutations", "#_of_those_DRM", "Mutations"]
        import pandas as pd
        print(pd.DataFrame(data=secondTen, columns=colN))

