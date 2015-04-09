"""
Class, which contains methods to optimize branch lengths given the time
constraints set to leaves
"""
from tree_anc import TreeAnc
import numpy as np
from Bio import AlignIO
import datetime
from scipy import stats
import matplotlib.pyplot as plt


class DateConversion(object):
    """
    Small container class to store parameters to convert between branch length as it is used in ML computations and the dates of the nodes.
    It is assumed that the conversion formula is 'length = k*date + b'
    """
    def __init__(self):
        self.k = 0
        self.b = 0
        self.r = 0
        self.pi = 0
        self.sigma = 0

    def get_branch_len(self, date1, date2):
        """
        Compute branch length given the dates of the two nodes.

        Args:
         - date1 (int): date of the first node (days before present)
         - date2 (int): date of the second node (days before present)

        Returns:
         - branch length (double): Branch length, assuming that the dependence between the node date and the node depth in the the tree is linear.
        """
        return abs(date1-date2)*self.k

    def get_date(self, node):
        """
        Get the approximate date of the tree node, assuming that the dependence between the node date and the node depth int the tree is linear.

        Args:
         - node(Phylo.Tree.Clade): node of the tree. Must be from the TreeAnc class (or its derivative), to contain the necessary attributes (dist2root).

        """
        yr = (self.b - node.dist2root)/self.k
        if (yr < 0):
            print ("Warning: got the negative date! Returning the inverse.")
            yr = abs(yr)
        return yr



class TreeTime(TreeAnc):

    def __init__(self, tree, dates_file=""):
        super(TreeTime,self).__init__(tree)
        self.date2dist = None # we do not know anything about the conversion
        self.tree_file = ""

    @classmethod
    def from_files(cls, tree_file, aln_file, dates_file="", tree_fmt='newick', aln_fmt='fasta'):
        t = cls.from_file(tree_file, tree_fmt)
        t.tree_file = tree_file

        aln = AlignIO.read(aln_file, aln_fmt)
        t.set_seqs_to_leaves(aln)
        if dates_file != "":
            d_dic = t._read_dates_file(dates_file)
            t._set_dates_to_all_nodes(d_dic) # even with empty dictionary we need to set dates as empty objects
        return t

    def load_dates(self, dates_file):
        """
        Args:
         - dates_file(str): path to file containing timing info about the tree nodes. File should be in csv format: 'nodeName,year'
        """
        d_dic = self._read_dates_file(dates_file)
        if (len(d_dic)) < 1: return
        self._set_dates_to_all_nodes(d_dic)

    def _read_dates_file(self, inf, **kwargs):
        """
        Read dates from the file into python dictionary. The input file should be in csv format 'node name, date'. The date will be converted to the datetime object and added to the dictionary {node name: datetime}
        """

        if 'verbose' in kwargs:
            verbose = kwargs['verbose']
        else:
            verbose = 10

        if (verbose > 3):
            print ("Reading datetime info for the tree nodes...")
        with open(inf, 'r') as finf:
            ss = finf.readlines()
        if verbose > 5:
            print ("Loaded %d lines form dates file" %len(ss))
        try:
            d = {s.split(',')[0]: datetime.datetime(int(s.split(',')[1].strip()),1,1) for s in ss}
            if verbose > 3:
                print ("Parsed data in %d lines of %d input, %d corrupted" %(len(d), len(ss), len(ss) - len(d)))
            return d
        except:
            ## unable to read all dates, the file is corrupted - go one by one
            print ("Unable to perform parsing of the dates file, file is corrupted. Return empty dictionary.")
            return {}

    def _set_dates_to_all_nodes(self, dates_dic):
        """
        Set the time information to all nodes.
        Gets the datetime object of the nodes specified, calculate the time before present  (in days) for the nodes and sets this parameter (as int) to the node.date attribute.
        Args:
         - dates_dic(dic): dictionary with datetime informationfor nodes. Format: {node name: datetime object}

        """
        now = datetime.datetime.now()

        for node in self.tree.find_clades(order='preorder'):
            if dates_dic.has_key(node.name):
                days_before_present = (now - dates_dic[node.name]).days
                if (days_before_present < 0):
                    print ("Cannot set the date ")
                    continue
                node.date = days_before_present
            else:
                node.date = None
        # FIXME for now, we just reroot the to the most ancient node. Later, it should be done in a cleverer way.
        self.tree.root_with_outgroup(sorted(self.tree.get_terminals(), key = lambda x:x.date)[-1])
        self.tree.root.date = None

        #fix tree lengths, etc
        self._add_node_params()

    def dates2dist_compute(self):
        """
        Get the conversion coefficients between the dates and the branch lengths as they are used in ML computations. The conversion formula is assumed to be 'length = k*date + b'. For convenience, these coefficients as well as regression parameters are stored in the dates2dist object.

        Note: that tree must have dates set to all nodes before calling this function. (The latter is accomplished by calling load_dates func).
        """

        self.date2dist = DateConversion()

        d = []
        for n in self.tree.find_clades():
            if n.date is not None:
                d.append((n.date, n.dist2root ))
        d = np.array(d)

        self.date2dist.k,\
        self.date2dist.b,\
        self.date2dist.r,\
        self.date2dist.pi,\
        self.date2dist.sigma =  stats.linregress(d[:, 0], d[:, 1])

    def date2dist_plot(self):
        """
        Plot the dependence between the node depth in the tree and the given node date information.
        """

        d = []
        for n in self.tree.find_clades():
            if n.date is not None:
                d.append((n.date, n.dist2root ))
        d = np.array(d)

        if self.date2dist is None:
            self.date2dist = DateConversion()

            self.date2dist.k,\
            self.date2dist.b,\
            self.date2dist.r,\
            self.date2dist.pi,\
            self.date2dist.sigma =  stats.linregress(d[:, 0], d[:, 1])

        plt.plot(d[:,0], d[:,1], 'o',c='r', fillstyle='none', markersize=9, label='Data')


        plt.plot([0, d[:,0].max()],
            [self.date2dist.b,self.date2dist.b+self.date2dist.k*d[:,0].max()],
            lw=2, c='r', label='Linear regression')

        plt.grid()
        plt.ylabel("Distance from root (node depth)")
        plt.xlabel("Node date (days before present)")
        plt.title ("Dependence between the node depth\nand the given date time constraints of the node.\n Tree file: %s" % self.tree_file)
        plt.legend()











