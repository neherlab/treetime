import Bio
from Bio import Phylo, Align, AlignIO
import numpy as np
import json, copy, datetime
from treetime import TreeTime
import utils 
import seq_utils
import pandas 
import pandas

def treetime_from_newick(infile):
    """
    Create TreeTime object and load phylogenetic tree from newick file
    Args:
     - infile(str): path to the newick file.
    Returns:
     - tanc(TreeTime): tree time object with phylogenetic tree set and required
     parameters assigned to the nodes.
    """
    tanc = TreeTime()
    tanc.tree = Phylo.read(infile, 'newick')
    tanc.set_additional_tree_params()
    return tanc

def save_newick_tree(outfile):
    pass

def save_timetree_results(tree, outfile_prefix):
    """
    First, it scans the tree and assigns the namesto every node with no name
    then, it saves the information as the csv table 
    """
    df = pandas.DataFrame(columns=["Given_date", "Initial_root_dist", "Inferred_date"])
    aln = Align.MultipleSeqAlignment([])
    
    for node in tree.tree.find_clades():
        
        if node.name is None:
            if node.up is None:
                node.name = "ROOT"
            else:
                node.name = "node_" + str(np.random.randint(1e8))
        
        #  alignment
        aln.append(Bio.SeqRecord.SeqRecord(Align.Seq(''.join(node.sequence)), node.name))
        #  metadata
        df.loc[node.name] = [node.numdate_given, node.dist2root_0, node.numdate]


    # save everything
    df.to_csv(outfile_prefix + ".meta.csv")
    Phylo.write(tree.tree, outfile_prefix + ".tree.nwk", "newick")
    AlignIO.write(aln, outfile_prefix + ".aln.fasta", "fasta")
    pass 

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
            l.sequence = dic_aln[l.name]
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
    now = utils.numdate(datetime.datetime.now())
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
    
    print ("Assigned ddates to {0} nodes, {1} errors".format(*tu))

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
    
    try:
        
        df = pandas.read_csv(infile, index_col=0)
        if df.index.name != "name" and df.index.name != "#name":
            print ("Cannot read metadata: first columns should be leaves name")
            return
        dic = df.to_dict(orient='index')
    except:

        print ("Cannot read the metadata using psndas library. "
            "pandas is outdated or missing. trying to read plain csv...")

    tree.set_metadata(**dic)    


if __name__=='__main__':
    pass
