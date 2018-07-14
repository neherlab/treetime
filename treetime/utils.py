from __future__ import division, print_function, absolute_import
import datetime
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy import stats
from scipy.ndimage import binary_dilation
from treetime import config as ttconf

def setup_TreeRegr(tt, covariation=True):
    from .treeregression import TreeRegression

    tip_value = lambda x:np.mean(x.numdate_given) if (x.is_terminal() and x.bad_branch is False) else None
    branch_value = lambda x:x.mutation_length
    if covariation:
        om = tt.one_mutation
        branch_variance = lambda x:(x.mutation_length+(ttconf.OVER_DISPERSION*om if x.is_terminal() else 0.0))*om
    else:
        branch_variance = lambda x:1.0 if x.is_terminal() else 0.0

    return TreeRegression(tt.tree, tip_value=tip_value,
                         branch_value=branch_value, branch_variance=branch_variance)


class DateConversion(object):
    """
    Small container class to store parameters to convert between branch length
    as it is used in ML computations and the dates of the nodes.
    It is assumed that the conversion formula is 'length = k*date + b'
    """
    def __init__(self):

        self.clock_rate = 0
        self.intercept = 0
        self.chisq = 0
        self.r_val = 0
        self.cov = 0
        self.sigma = 0

    def __str__(self):
        # TODO: fix the chi^2 vs r^2 output
        outstr = ('Root-Tip-Regression:\n --rate:\t%1.3e\n --chi^2:\t%1.2f\n\n --r^2:\t%1.2f\n'
                  %(self.clock_rate, self.chisq**2, self.r_val**2))
        return outstr


    @classmethod
    def from_regression(cls, clock_model):
        """
        Create the conversion object automatically from the tree

        Parameters
        ----------

         clock_model : dict
            dictionary as returned from TreeRegression with fields intercept and slope

        """
        dc = cls()
        dc.clock_rate = clock_model['slope']
        dc.intercept = clock_model['intercept']
        dc.chisq = clock_model['chisq']
        if 'cov' in clock_model:
            dc.cov = clock_model['cov']
        dc.r_val = clock_model['r_val']
        return dc

    def get_branch_len(self, date1, date2):
        """
        Compute branch length given the dates of the two nodes.

        Parameters
        -----------

         date1 : int
            date of the first node (days before present)

         date2 : int
            date of the second node (days before present)

        Returns:
        --------

         branch length : double
            Branch length, assuming that the dependence
            between the node date and the node depth in the the tree is linear.

        """
        return abs(date1 - date2) * self.clock_rate

    def get_time_before_present(self, numdate):
        """
        Convert the numeric date to the branch-len scale
        """
        return (numeric_date() - numdate) * abs(self.clock_rate)

    def to_years(self, abs_t):
        """
        Convert the time before present measured in branch length units to years

        """
        return abs_t / abs(self.clock_rate)

    def to_numdate(self, tbp):
        """
        Convert the numeric date to the branch-len scale
        """
        return numeric_date() - self.to_years(tbp)

    def numdate_from_dist2root(self, d2r):
        """
        estimate the numerical date based on the distance to root.
        -> crude dating of internal nodes
        """
        return (d2r-self.intercept)/self.clock_rate


def min_interp(interp_object):
    """
    Find the global minimum of a function represented as an interpolation object.
    """
    try:
        return interp_object.x[interp_object(interp_object.x).argmin()]
    except Exception as e:
        s = "Cannot find minimum of the interpolation object" + str(interp_object.x) + \
        "Minimal x: " + str(interp_object.x.min()) + "Maximal x: " + str(interp_object.x.max())
        raise e


def median_interp(interp_object):
    """
    Find the median of the function represented as an interpolation object.
    """
    new_grid = np.sort(np.concatenate([interp_object.x[:-1] + 0.1*ii*np.diff(interp_object.x)
                                       for ii in range(10)]).flatten())

    tmp_prop = np.exp(-(interp_object(new_grid)-interp_object.y.min()))
    tmp_cumsum = np.cumsum(0.5*(tmp_prop[1:]+tmp_prop[:-1])*np.diff(new_grid))
    median_index = min(len(tmp_cumsum)-3, max(2,np.searchsorted(tmp_cumsum, tmp_cumsum[-1]*0.5)+1))
    return new_grid[median_index]


def numeric_date(dt=None):
    """
    Convert datetime object to the numeric date.
    The numeric date format is YYYY.F, where F is the fraction of the year passed

    Parameters
    ----------
     dt:  datetime.datetime, None
        date of to be converted. if None, assume today

    """
    if dt is None:
        dt = datetime.datetime.now()

    try:
        res = dt.year + dt.timetuple().tm_yday / 365.25
    except:
        res = 0.0

    return res

def tree_layout(tree):
    leaf_count=0
    for ni,node in enumerate(tree.find_clades(order="postorder")):
        if node.is_terminal():
            leaf_count+=1
            node.ypos=leaf_count
        else:
            tmp = np.array([c.ypos for c in node])
            node.ypos=0.5*(np.max(tmp) + np.min(tmp))

def tree_inference(aln_fname, tree_fname, tmp_dir=None,
                   methods = ['iqtree', 'fasttree', 'raxml'], **kwargs):
    import os,shutil
    from Bio import Phylo
    if not os.path.isfile(aln_fname):
        print("alignment file does not exist")

    cwd = os.getcwd()
    if tmp_dir:
        if not os.path.isdir(tmp_dir):
            try:
                os.makedirs(tmp_dir)
            except OSError as e:
                print("Cannot create run_dir",e)
        aln_fname_base = os.path.basename(aln_fname)
        shutil.copyfile(aln_fname,os.path.join(tmp_dir, aln_fname_base))
        aln_fname = aln_fname_base
        os.chdir(tmp_dir)

    for method in methods:
        T = None
        try:
            if method.lower()=='iqtree':
                T = build_newick_iqtree(aln_fname)
            elif method.lower()=='fasttree':
                T = build_newick_fasttree(aln_fname, nuc=True)
            elif method.lower()=='raxml':
                T = build_newick_raxml(aln_fname)
            else:
                print("Method not supported",method)
            if T:
                break
        except:
            continue
    os.chdir(cwd)
    if T is None:
        print("tree building failed. tried", ", ".join(methods), "but none worked")
    else:
        Phylo.write(T, tree_fname, 'newick')


def build_newick_fasttree(aln_fname, nuc=True):
    import os
    from Bio import Phylo
    print("Building tree with fasttree")
    tree_cmd = ["fasttree"]
    if nuc: tree_cmd.append("-nt")

    tree_cmd.extend([aln_fname,"1>","tmp.nwk", "2>", "fasttree_stderr"])
    os.system(" ".join(tree_cmd))
    return Phylo.read("tmp.nwk", 'newick')


def build_newick_raxml(aln_fname, nthreads=2, raxml_bin="raxml", **kwargs):
    import shutil,os
    from Bio import Phylo, AlignIO
    AlignIO.write(AlignIO.read(aln_fname, 'fasta'),"temp.phyx", "phylip-relaxed")
    cmd = raxml_bin + " -f d -T " + str(nthreads) + " -m GTRCAT -c 25 -p 235813 -n tre -s temp.phyx"
    os.system(cmd)
    return Phylo.read('RAxML_bestTree.tre', "newick")


def build_newick_iqtree(aln_fname, nthreads=2, iqtree_bin="iqtree",
                        iqmodel="HKY",  **kwargs):
    import os
    from Bio import Phylo, AlignIO
    with open(aln_fname) as ifile:
        tmp_seqs = ifile.readlines()
    aln_file = "temp.fasta"
    with open(aln_file, 'w') as ofile:
        for line in tmp_seqs:
            ofile.write(line.replace('/', '_X_X_').replace('|','_Y_Y_'))

    fast_opts = [
        "-ninit", "2",
        "-n",     "2",
        "-me",    "0.05"
    ]

    call = ["iqtree", *fast_opts, "-nt", str(nthreads), "-s", aln_file, "-m", iqmodel,
            ">", "iqtree.log"]

    os.system(" ".join(call))
    T = Phylo.read(aln_file+".treefile", 'newick')
    for n in T.get_terminals():
        n.name = n.name.replace('_X_X_','/').replace('_Y_Y_','|')
    return T

if __name__ == '__main__':
    pass

