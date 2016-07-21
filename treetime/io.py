import Bio
from Bio import Phylo, Align, AlignIO
import numpy as np
import json, copy, datetime
from treetime import TreeTime
import utils
import seq_utils
import os
import StringIO


def treetime_from_newick(gtr, infile):
    """
    Create TreeTime object and load phylogenetic tree from newick file
    Args:
     - infile(str): path to the newick file.
    Returns:
     - tanc(TreeTime): tree time object with phylogenetic tree set and required
     parameters assigned to the nodes.
    """
    tanc = TreeTime(gtr)
    tanc.tree = Phylo.read(infile, 'newick')
    tanc.set_additional_tree_params()
    return tanc

def treetime_to_newick(tt, outf):
    Phylo.write(tt.tree, outf, 'newick')

def _layout(tree):
    """Add clade, xvalue, yvalue, mutation and trunk attributes to all nodes in tree"""
    tree.root.branch_length = 0.01
    clade = 0
    yvalue = 0
    for node in tree.find_clades(order="preorder"):
        # set mutations
        if node.up is not None:
            node.muts = ', '.join([node.up.sequence[p] + str(p) + node.sequence[p]
                for p in np.where(node.up.sequence != node.sequence)[0]])

        # set sequences
        node.strseq = "".join(node.sequence)

        # set clade No
        node.clade = clade
        clade += 1
        if node.up is not None: #try:
            # Set xValue, tValue, yValue
            node.xvalue = node.up.xvalue+node.opt_branch_length
            node.tvalue = node.numdate - tree.root.numdate
        else:
            node.xvalue = 0.0
            node.tvalue = 0.0
        if node.is_terminal():
            node.yvalue = yvalue
            yvalue += 1
        # check numdate
        if not hasattr(node, 'numdate'):
            node.numdate = 0.0
    for node in tree.get_nonterminals(order="postorder"):
        node.yvalue = np.mean([x.yvalue for x in node.clades])

def treetime_to_json(tt, outf):

    def _node_to_json(node):

        tree_json = {}
        str_attr = ['clade','strain', 'date', 'muts', 'strseq']
        num_attr = ['xvalue', 'yvalue', 'tvalue', 'numdate']

        if hasattr(node, 'name'):
            tree_json['strain'] = node.name
            tree_json['name'] = node.name

        for prop in str_attr:
            if hasattr(node, prop):
                tree_json[prop] = node.__getattribute__(prop)
        for prop in num_attr:
            if hasattr(node, prop):
                try:
                    tree_json[prop] = round(node.__getattribute__(prop),5)
                except:
                    print "cannot round:", node.__getattribute__(prop), "assigned as is"
                    tree_json[prop] = node.__getattribute__(prop)

        if node.clades: # node is internal
            tree_json["children"] = []
            for ch in node.clades:
                tree_json["children"].append(_node_to_json(ch))
        else:
            # node is terminal, set both terminal and internal metadata
            tree_json["terminal_metadata"] = [{'name': k.name, 'value': k.attr(node)} for k in tt._terminal_metadata_names]

        return tree_json

    _layout(tt.tree)
    tree_json = _node_to_json(tt.tree.root)
    with open (outf,'w') as of:
        json.dump(tree_json, of, indent=False)

def tips_data_to_json(tt, outf):

    if not hasattr(tt.tree.get_terminals()[0], 'xvalue'):
        _layout(tt.tree);

    arr = [
    {
        'name': k.name,
        'strain':k.name,
        'numdate_given': k.numdate_given if hasattr(k, 'numdate_given') else 0.0,
        'numdate': k.numdate if hasattr(k, 'numdate') else 0.0,
        'xValue': k.xvalue if hasattr(k, 'xvalue') else 0.0,

    } for k in tt.tree.get_terminals()]

    with open (outf,'w') as of:
        json.dump(arr, of, indent=True)

def root_pos_lh_to_human_readable(tt, cutoff=1e-4):

    mtp = mtp = tt.tree.root.msg_to_parent
    mtp_min = mtp.y.min()

    mtpy = np.array([np.exp(-k+mtp_min) for k in mtp.y])
    mtpx = mtp.x

    # cut and center
    maxy_idx = mtpy.argmax()
    val_right = (mtpy[maxy_idx:] > cutoff)
    if (val_right.sum() == 0):
        right_dist = 0
    else:
        # left, actually (time is in the opposite direction)
        right_dist = - mtpx[maxy_idx] + mtpx[maxy_idx + val_right.argmin()]

    val_left = mtpy[:maxy_idx] > cutoff
    if (val_left.sum() == 0):
        left_dist = 0.0
    else:
        left_dist =  mtpx[maxy_idx] - mtpx[maxy_idx - (maxy_idx - val_left.argmax())]


    dist = np.max((left_dist, right_dist))
    center = mtpx[maxy_idx]

    # final x-y scatter
    #import ipdb; ipdb.set_trace()
    raw_x = np.unique(np.concatenate(([center-dist], [center], [center+dist], mtpx[(mtpx < dist + center) & (mtpx > center-dist)])))

    x = utils.numeric_date() -  np.array(map(tt.date2dist.get_date, raw_x))
    y = np.exp(-(mtp(raw_x) - mtp_min))
    return x, y

def root_lh_to_json(tt, outf):

    x,y = root_pos_lh_to_human_readable(tt)
    arr = [{"x":f, "y":b} for f, b in zip(x, y)]

    with open (outf,'w') as of:
        json.dump(arr, of, indent=True)

    print (', '.join([str(k) for k in x]))
    print (', '.join([str(k) for k in y]))

    #import ipdb; ipdb.set_trace()
def root_lh_to_csv(tt, outf):
    """Save node position likelihood distribution to CSV file"""
    x,y = root_pos_lh_to_human_readable(tt)
    arr = np.zeros((x.shape[0], 2))
    arr[:, 0] = x[:]
    arr[:, 1] = y[:]
    np.savetxt(outf, arr, delimiter=",",header="#Numdate,Likelihood_normalized")



def save_all_nodes_metadata(tt, outfile):

    import pandas
    metadata = tt._terminal_metadata_names
    d = [[k.attr(n) for k in metadata] for n in tt.tree.find_clades()]
    df = pandas.DataFrame(d, index=[k.name for k in tt.tree.find_clades()], columns=[k.name for k in metadata])
    df.sort_index(inplace=True)
    df.to_csv(outfile)

def save_timetree_results(tree, outfile_prefix):
    """
    First, it scans the tree and assigns the namesto every node with no name
    then, it saves the information as the csv table
    """
    import pandas
    df = pandas.DataFrame(columns=["Given_date", "Initial_root_dist", "Inferred_date"])
    aln = Align.MultipleSeqAlignment([])

    i = 0

    # save everything
    df.to_csv(outfile_prefix + ".meta.csv")
    #  TODO save variance to the metadata
    Phylo.write(tree.tree, outfile_prefix + ".tree.nwk", "newick")
    AlignIO.write(aln, outfile_prefix + ".aln.fasta", "fasta")

    # save root distibution
    mtp = tree.tree.root.msg_to_parent
    threshold = mtp.y.min() + 1000
    idxs = [mtp.y < threshold]
    mtpy = mtp.y[idxs]
    mtpx = utils.numeric_date() -  np.array(map(tree.date2dist.get_date, mtp.x[idxs]))
    mtpy[0] = threshold
    mtpy[-1] = threshold

    np.savetxt(outfile_prefix + ".root_dist.csv",
        np.hstack((mtpx[:, None], mtpy[:, None])),
        header="Root date,-log(LH)", delimiter=',')

    # zip results to one file
    import zipfile
    outzip = outfile_prefix + ".zip"
    zipf = zipfile.ZipFile(outzip, 'w')
    zipf.write(outfile_prefix + ".meta.csv")
    zipf.write(outfile_prefix + ".aln.fasta")
    zipf.write(outfile_prefix + ".tree.nwk")
    zipf.write(outfile_prefix + ".root_dist.csv")

def set_seqs_to_leaves(tree, aln):
    """
    Set sequences from the alignment to the leaves of the tree of the TreeAnc class.
    Args:
     - tree (TreeAnc): instance of the treeAnc class with the tree loaded. The
     names of the tree leaves must match exactly with those of the alignment
     sequences.
     - aln(Bio.MultipleSequenceAlignment): alignment ogbject
    Returns:
     - failed_leaves(int): number of leaves which could not be assigned with
     sequences.
    Note:
     - If there are more than 100 leaves failed to get sequences, the function
     breaks, returning 100.
    """
    failed_leaves= 0
    dic_aln = {k.name: seq_utils.prepare_seq(k.seq) for k in aln} #
    for l in tree.tree.get_terminals():
        if l.name in dic_aln:
            l.state_seq = dic_aln[l.name]
            l.sequence=l.state_seq
        else:
            print ("Cannot find sequence for leaf: %s" % l.name)
            failed_leaves += 1
            if failed_leaves == 100:
                print ("Error: cannot set sequences to the terminal nodes.\n"
                    "Are you sure the alignment belongs to the tree?")
                break
    tree.L = aln.get_alignment_length()
    tree.one_mutation = 1.0/tree.L
    tree.tree.root.branch_length = tree.one_mutation
    return failed_leaves

def read_dates_file(self, inf, **kwargs):
        """
        Read dates from the file into python dictionary. The input file should
        be in csv format 'node name, date'. The date will be converted to the
        datetime object and added to the dictionary {node name: datetime}

        Args:
         - inf(str): path to input file

        KWargs:
         - verbose(int): how verbose should the output be

        Returns:
         - dic(dic): dictionary  {NodeName: Date as datetime object}
        """

        def str_to_date(instr):
            """
            Convert input string to datetime object.

            Args:
             - instr (str): input string. Accepts one of the formats:
             {YYYY.MM.DD, YYYY.MM, YYYY}.

            Returns:
             - date (datetime.datetime): parsed date object. If the parsing
             failed, None is returned
            """
            try:
                date = datetime.datetime.strptime(instr, "%Y.%m.%d")
            except ValueError:
                date = None
            if date is not None:
                return date

            try:
                date = datetime.datetime.strptime(instr, "%Y.%m")
            except ValueError:
                date = None

            if date is not None:
                return date

            try:
                date = datetime.datetime.strptime(instr, "%Y")
            except ValueError:
                date = None

            return date

        if 'verbose' in kwargs:
            verbose = kwargs['verbose']
        else:
            verbose = 10

        if verbose > 3:
            print ("Reading datetime info for the tree nodes...")
        with open(inf, 'r') as finf:
            all_ss = finf.readlines()
        if verbose > 5:
            print ("Loaded %d lines form dates file" % len(all_ss))
        try:
            dic = {s.split(',')[0]: str_to_date(s.split(',')[1].strip())
                   for s in all_ss if not s.startswith("#")}
            if verbose > 3:
                print ("Parsed data in %d lines of %d input, %d corrupted"
                       % (len(dic), len(all_ss), len(all_ss) - len(dic)))
            return dic
        except ValueError:
            # unable to read all dates, the file is corrupted - go one by one
            print ("Unable to perform parsing of the dates file, file is "
                   "corrupted. Return empty dictionary.")
            return {}

#  FIXME make consistent with metadata dictionary
def set_node_dates_from_dic(tree, dates_dic):
    """
    Read names of the leaves of the tree, mathc them with the provided dictionary
    and set the raw_date attribute to the nodes. If the dictionary has no entry
    for  a node, the node gets raw_date = None attribute.
    Args:
     - tree (TreeTime): instance of the tree time object with phylogenetic tree
     loaded.
     - dates_dic (dic): dictionary storing dates of the nodes as datetime.datetime
     object.
    Returns:
     - None, tree is being modified in-place
    """

    err_ = 0
    num_ = 0
    now = utils.numeric_date(datetime.datetime.now())
    for node in tree.tree.find_clades():

        if node.name is None or not node.name in dates_dic:
            node.numdate_given = None
            continue

        n_date = dates_dic[node.name] # assume the dictionary contains the numdate
        if not isinstance(n_date, float) and not isinstance(n_date, int): #  sanity check
            print ("Cannot set the numeric date tot the node. Float or int expected")
            continue

        try:

            if n_date > now:
                print ("Cannot set the date! the specified date is later "
                    " than today! cannot assign node date, skipping")
                node.numdate_given = None
                err_+=1
                continue
            else:
                node.numdate_given = n_date
                num_ += 1

        except:
            print ("Cannot assign date to the node: exception caught")
            node.numdate_given = None
            err_ += 1

    tu = (num_, err_)

    print ("Assigned dates to {0} nodes, {1} errors".format(*tu))

def set_node_dates_from_names(tree, date_func):
    """
    Read names of the leaves of the tree, extract the dates of the leaves from the
    names and asign the date to the nodes.
    Assumes that the dates are given in some human-readable format
    and are converted into the numericaldate (YYYY.F).
    After this function call, each node of
    the tree gets the numdate_given attribute. If the date was extracted from name
    successfully, the 'numdate_given' will be the days-before-present (int) value.
    Otherwise (either no node name, or date reading failed), the 'numdate_given' will be
    set to None.
    Args:
     - tree (TreeTime): instance of the tree time object with phylogenetic tree
     loaded.
     - date_func (callable): function to extract date and time from node name,
     should return float
    Returns:
     - None, tree is being modified in-place
    """
    #now = datetime.datetime.now()
    ## NOTE we do not rely on the datetime objects
    now = utils.numeric_date(datetime.datetime.now())
    for node in tree.tree.get_terminals():
        try:
            node_date = date_func(node.name)
        except:
            print ("Cannot extract numdate from the node name. Exception caugth.")
            node.numdate_given = None
            continue
        if node_date is None:
            #print ("Cannot parse the date from name: " + str(node.name) +
            #    " Setting node raw date to None")
            node.numdate_given = None # cannot extract the date from name - set None

        elif node_date > now:
            print ("Cannot set the date! the specified date is later "
                " than today")
            node.numdate_given = None
        else:
            node.numdate_given = node_date

    return

def read_metadata(tree, infile):
    if os.path.isfile(infile):
        try:
            import pandas
            df = pandas.read_csv(infile, index_col=0, sep=r'\s*,\s*')
            if df.index.name != "name" and df.index.name != "#name":
                print ("Cannot read metadata: first column should contain the leaves names")
                return

            potential_date_columns = []
            potential_numdate_columns = []

            for ci,col in enumerate(df.columns):
                if 'date' in col.lower():
                    try: #avoid date parsing when can be parsed as float
                        tmp = float(df.iloc[0,ci])
                        potential_numdate_columns.append((ci, col))
                    except: #otherwise add as potential date column
                        potential_date_columns.append((ci, col))


            # if a potential date column was found
            if len(potential_numdate_columns)>=1:
                # use the first numdate column
                #idx = potential_numdate_columns[0][0]
                name = potential_numdate_columns[0][1]
                # Use this column as numdate_given
                df.rename_axis({name:"numdate_given"}, axis=1, inplace=True)

            elif len(potential_date_columns)>=1:

                #try to parse the csv file with dates in the idx column:
                idx = potential_date_columns[0][0]
                name = potential_date_columns[0][1]

                # NOTE as the 0th column is the index, we should parse the dates
                # for the column idx + 1
                df = pandas.read_csv(infile, index_col=0, sep=r'\s*,\s*', parse_dates=[1+idx])
                #convert to numdate
                df[name] = map (utils.numeric_date, df[name])
                # use this column as the numeric date:
                df.rename_axis({name:"numdate_given"}, axis=1, inplace=True)

            else:
                print ("Metadata file has nothing which looks like a sampling date!")

            dic = df.to_dict(orient='index')
            tree.set_metadata(**dic)
        except:
            print ("Cannot read the metadata using the pandas library. "
                "pandas is outdated or missing., reading csv directly")
            import csv
            dic = {}
            with open(infile) as ifile:
                for row in csv.DictReader(ifile):
                    dic[row['name']]=row
                    if 'numdate_given' in row:
                        dic[row['name']]['numdate_given'] = float(dic[row['name']]['numdate_given'])
            tree.set_metadata(**dic)
    else:
        print("meta data file not found!")


def save_gtr_to_file(gtr, outfile):


    with open(outfile, 'w') as of:
        of.write("#GTR alphabet\n" + ','.join(gtr.alphabet)+'\n')

    with open (outfile,'a') as of:
        of.write("#Mutation rate:\n" + str(gtr.mu) + '\n')

    with open(outfile, 'a') as of:
        np.savetxt(of, gtr.mu * np.dot(gtr.Pi, gtr.W), header="Full GTR matrix", delimiter=",")

    with open(outfile, 'a') as of:
        np.savetxt(of, np.diag(gtr.Pi), delimiter=",", header="Equilibrium character composition")

    with open(outfile, 'a') as of:
        np.savetxt(of, gtr.W, delimiter=",", header="Flow rate matrix")

if __name__=='__main__':
    pass
