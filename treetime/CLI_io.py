import os, sys
from Bio import AlignIO, Phylo
from .vcf_utils import read_vcf, write_vcf
from .seq_utils import alphabets
from Bio import __version__ as bioversion
from . import version as treetime_version
import numpy as np

def get_outdir(params, suffix='_treetime'):
    if params.outdir:
        if os.path.exists(params.outdir):
            if os.path.isdir(params.outdir):
                return params.outdir.rstrip('/') + '/'
            else:
                print("designated output location %s is not a directory"%params.outdir, file=sys.stderr)
        else:
            os.makedirs(params.outdir)
            return params.outdir.rstrip('/') + '/'

    from datetime import datetime
    outdir_stem = datetime.now().date().isoformat()
    outdir = outdir_stem + suffix.rstrip('/')+'/'
    count = 1
    while os.path.exists(outdir):
        outdir = outdir_stem + '-%04d'%count + suffix.rstrip('/')+'/'
        count += 1

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
                              report_ambiguous=False, timetree=False, confidence=False,
                              reconstruct_tip_states=False, tree_suffix=''):
    seq_info = is_vcf or tt.aln
    if is_vcf:
        outaln_name = basename + f'ancestral_sequences{tree_suffix}.vcf'
        write_vcf(tt.get_reconstructed_alignment(reconstruct_tip_states=reconstruct_tip_states), outaln_name)
    elif tt.aln:
        outaln_name = basename + f'ancestral_sequences{tree_suffix}.fasta'
        AlignIO.write(tt.get_reconstructed_alignment(reconstruct_tip_states=reconstruct_tip_states), outaln_name, 'fasta')
    if seq_info:
        print("\n--- alignment including ancestral nodes saved as  \n\t %s\n"%outaln_name)

    # decorate tree with inferred mutations
    terminal_count = 0
    offset = 0 if zero_based else 1
    if timetree:
        dates_fname = basename + f'dates{tree_suffix}.tsv'
        fh_dates = open(dates_fname, 'w', encoding='utf-8')
        if confidence:
            fh_dates.write('#Lower and upper bound delineate the 90% max posterior region\n')
            fh_dates.write('#node\tdate\tnumeric date\tlower bound\tupper bound\n')
        else:
            fh_dates.write('#node\tdate\tnumeric date\n')

    mutations_out = open(basename + "branch_mutations.txt", "w")
    mutations_out.write("node\tstate1\tpos\tstate2\n")
    for n in tt.tree.find_clades():
        if timetree:
            if confidence:
                if n.bad_branch:
                    fh_dates.write('%s\t--\t--\t--\t--\n'%(n.name))
                else:
                    conf = tt.get_max_posterior_region(n, fraction=0.9)
                    fh_dates.write('%s\t%s\t%f\t%f\t%f\n'%(n.name, n.date, n.numdate,conf[0], conf[1]))
            else:
                if n.bad_branch:
                    fh_dates.write('%s\t--\t--\n'%(n.name))
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
            if n.mask is None:
                if report_ambiguous:
                    n.comment= '&mutations="' + ','.join([a+str(pos + offset)+d for (a,pos, d) in n.mutations])+'"'
                else:
                    n.comment= '&mutations="' + ','.join([a+str(pos + offset)+d for (a,pos, d) in n.mutations
                                                        if tt.gtr.ambiguous not in [a,d]])+'"'
            else:
                if report_ambiguous:
                    n.comment= '&mutations="' + ','.join([a+str(pos + offset)+d for (a,pos, d) in n.mutations if n.mask[pos]>0])+f'",mcc="{n.mcc}"'
                else:
                    n.comment= '&mutations="' + ','.join([a+str(pos + offset)+d for (a,pos, d) in n.mutations
                                                        if tt.gtr.ambiguous not in [a,d] and n.mask[pos]>0])+f'",mcc="{n.mcc}"'

                for (a, pos, d) in n.mutations:
                    if tt.gtr.ambiguous not in [a,d] or report_ambiguous:
                        mutations_out.write("%s\t%s\t%s\t%s\n" %(n.name, a, pos + 1, d))
        if timetree:
            n.comment+=(',' if n.comment else '&') + 'date=%1.2f'%n.numdate
    mutations_out.close()

    # write tree to file
    fmt_bl = "%1.6f" if tt.data.full_length<1e6 else "%1.8e"
    if timetree:
        outtree_name = basename + f'timetree{tree_suffix}.nexus'
        print("--- saved divergence times in \n\t %s\n"%dates_fname)
        Phylo.write(tt.tree, outtree_name, 'nexus')
    else:
        outtree_name = basename + f'annotated_tree{tree_suffix}.nexus'
        Phylo.write(tt.tree, outtree_name, 'nexus', format_branch_length=fmt_bl)
    print("--- tree saved in nexus format as  \n\t %s\n"%outtree_name)

    auspice = create_auspice_json(tt, timetree=timetree, confidence=confidence)
    outtree_name_json = basename + f'auspice_tree{tree_suffix}.json'
    with open(outtree_name_json, 'w') as fh:
        import json
        json.dump(auspice, fh, indent=0)
        print("--- tree saved in auspice json format as  \n\t %s\n"%outtree_name_json)

    if timetree:
        for n in tt.tree.find_clades():
            n.branch_length = n.mutation_length
        outtree_name = basename + f'divergence_tree{tree_suffix}.nexus'
        Phylo.write(tt.tree, outtree_name, 'nexus', format_branch_length=fmt_bl)
        print("--- divergence tree saved in nexus format as  \n\t %s\n"%outtree_name)


def print_save_plot_skyline(tt, n_std=2.0, screen=True, save='', plot='', gen=50):
    if plot:
        import matplotlib.pyplot as plt

    skyline, conf = tt.merger_model.skyline_inferred(gen=gen, confidence=n_std)
    if save: fh = open(save, 'w', encoding='utf-8')
    header1 = "Skyline assuming "+ str(gen)+" gen/year and approximate confidence bounds (+/- %f standard deviations of the LH)\n"%n_std
    header2 = "date \tN_e \tlower \tupper"
    if screen: print('\t'+header1+'\t'+header2)
    if save: fh.write("#"+ header1+'#'+header2+'\n')
    for (x,y, y1, y2) in zip(skyline.x, skyline.y, conf[0], conf[1]):
        if screen: print("\t%1.3f\t%1.3e\t%1.3e\t%1.3e"%(x,y, y1, y2))
        if save: fh.write("%1.3f\t%1.3e\t%1.3e\t%1.3e\n"%(x,y, y1, y2))

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



def create_auspice_json(tt, timetree=False, confidence=False):
    # mock up meta data for auspice json
    from datetime import datetime
    meta = {
        "title": f"Auspice visualization of TreeTime (v{treetime_version}) analysis",
        "build_url": "https://github.com/neherlab/treetime",
        "last_updated": datetime.now().strftime("%Y-%m-%d"),
        "treetime_version": treetime_version,
        "genome_annotations": {
            "nuc":{"start":1, "end":int(tt.data.full_length), "type":"source", "strand":"+:"}
        },
        "panels":["tree", "entropy"],
        "colorings": [
            {
                "title": "Date",
                "type": "continuous",
                "key": "num_date",
            },
            {
                "title": "Genotype",
                "type": "categorical",
                "key": "gt",
            },
            {
                "title": "Excluded",
                "type": "categorical",
                "key": "bad_branch"
            },
            {
                "title": "Branch Support",
                "type": "continuous",
                "key": "confidence"
            }
        ],
        "display_defaults": {"color_by":"bad_branch"},
        "filters": ["bad_branch"]
    }

    def node_to_json(n, pdiv=0.0):
        j = {"name":n.name, "node_attrs":{}, "branch_attrs":{}}
        if n.clades:
            j["children"] = []

        if timetree:
            j["node_attrs"]["num_date"] = {"value":float(n.numdate)}
            if confidence:
                conf = tt.get_max_posterior_region(n, fraction=0.9)
                j["node_attrs"]["num_date"]["confidence"] = (float(conf[0]), float(conf[1]))
        j["node_attrs"]["div"] = float(pdiv + n.mutation_length)
        j["node_attrs"]["bad_branch"] = {"value": "Yes" if n.bad_branch else "No"}

        j["branch_attrs"]["mutations"] = {"nuc": [f"{a}{pos+1}{d}" for a,pos,d in n.mutations if d in "ACGT-"]}
        # generate bootstrap confidence substitute via the negative exponential of the number of mutations
        # this is the bootstrap confidence for iid mutations (only ACGT mutations)
        j["node_attrs"]["confidence"] = {"value":round(1-np.exp(-len([pos for a,pos,d in n.mutations if d in "ACGT"])),3)
                                          if not n.is_terminal() else 1.0}
        return j

    # create the tree data structure from the Biopython tree
    tree = node_to_json(tt.tree.root, 0.0)
    # dictionary to look up nodes by name
    node_lookup = {tt.tree.root.name: tree}
    for n in tt.tree.get_nonterminals():
        n_json = node_lookup[n.name]
        for c in n.clades:
            # generate node jsons for all children and attach them the to parent
            n_json["children"].append(node_to_json(c, n_json["node_attrs"]["div"]))
            node_lookup[c.name] = n_json["children"][-1]

    return {"meta":meta, "tree":tree}
