from Bio import Phylo
import numpy as np
import json, copy, datetime
from tree_time import TreeTime
import seq_utils

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

def treetime_from_json(infile):
    raise NotImplementedError("this feature is currently under development")

def _read_json_tree(node, json_clade, data_keys, date_func):
        """
        recursive function to populate tree from the json strcture
        Args:
         - json_clade(dic): the data for the tree node represented as dictionary

         - data_keys(dic): dictionary to convert (some) data keys into the internal
         tree_time notification

         - date_func(callable): function to convert string date in the json into
         the datetime object of the tree node
        """

        clade_key = 'clade' if 'clade' not in data_keys else  data_keys['clade']
        if clade_key in json_clade:
            node.clade = json_clade[clade_key]

        name_key = 'name' if 'name' not in data_keys else  data_keys['name']
        if name_key in json_clade:
            node.name = json_clade[name_key]

        strain_key = 'strain' if 'strain' not in data_keys else data_keys['strain']
        if strain_key in json_clade:
            node.strain = json_clade[strain_key]

        f_key = 'branch_length' if 'branch_length' not in data_keys else data_keys['branch_length']

        if f_key in json_clade:
            node.branch_length = float(json_clade[f_key])


        f_key = 'xvalue' if 'xvalue' not in data_keys else data_keys['xvalue']

        if f_key in json_clade:
            node.xvalue = float(json_clade[f_key])

        f_key = 'yvalue' if 'yvalue' not in data_keys else data_keys['yvalue']
        if f_key in json_clade:
            node.yvalue = float(json_clade[f_key])

        f_key = 'days_before_present' if 'days_before_present' not in data_keys else data_keys['days_before_present']
        if f_key in json_clade:
            node.date=float(json_clade[f_key])

        f_key = 'seq' if 'seq' not in data_keys else data_keys['seq']
        if f_key in json_clade:
            node.sequence = np.array(list(json_clade[f_key]))

        f_key = 'yvalue' if 'yvalue' not in data_keys else data_keys['yvalue']
        if f_key in json_clade:
            node.yvalue = float(json_clade[f_key])

        f_key = 'logLH' if 'logLH' not in data_keys else data_keys['logLH']
        if f_key in json_clade:
            node.logLH = float(json_clade[f_key])

        if len(node.clades):
            json["children"] = []
            for ch in node.clades:
                json["children"].append(self.to_json(ch))

def from_json(cls, inf, json_keys={"branch_len":"xvalue"}, date_func=lambda x: None):
    """
    Load tree from json file.

    Args:
     - inf(str): pth to the input file
     - json_keys(dic): names of the parameters in the json file. The names
     for the following parameters should be specified:
     {"date": ..., "branch_len": ...}
     - date_func(callable): function to convert the date representation from
     json parameter string to datetime object (will be assigned as raw_date)
    Returns:
     - TreeTime object
    """
    with open (inf) as json_file:
        data = json.load(json_file)

    if len(data) < 1 or 'children' not in data:
        raise IOError("Wrong format of json file")
    t = Phylo.BaseTree.Tree()
    ttime = TreeTime()

    ttime.read_json_tree(t.root, data, json_keys, date_func)

    return ttime

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
    raise NotImplementedError("This feature is currently under development")

def set_node_dates_from_names(tree, date_func):
    """
    Read names of the leaves of the tree, extract the dates of the leaves from the
    names and asign the date to the nodes. After this function call, each node of
    the tree gets the raw_date attribute. If the date was extracted from name
    successfully, the raw_date will be the days-before-present (int) value.
    Otherwise (either no node name, or date reading failed), the raw_date will be
    set to None.
    Args:
     - tree (TreeTime): instance of the tree time object with phylogenetic tree
     loaded.
     - date_func (callable): function to extract date and time from node name,
     should return datetime object
    Returns:
     - None, tree is being modified in-place
    """
    now = datetime.datetime.now()
    for node in tree.tree.find_clades():
        try:
            node_date = date_func(node.name)
            if node_date is None:
                node.raw_date = None # cannot extract the date from name - set None
                continue
            days_before_present = (now - node_date).days
            if days_before_present < 0:
                print ("Cannot set the date! the specified date is later "
                    " than today")
                node.raw_date = None
                continue
            node.raw_date = days_before_present
        except:
            node.raw_date = None
    return

if __name__=='__main__':
    pass
