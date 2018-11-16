from __future__ import division, print_function, absolute_import
import os,sys
import datetime
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy import stats
from scipy.ndimage import binary_dilation
from treetime import config as ttconf

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
        self.cov = None
        self.sigma = 0
        self.valid_confidence = False

    def __str__(self):
        if self.cov is not None and self.valid_confidence:
            dslope = np.sqrt(self.cov[0,0])
            outstr = ('Root-Tip-Regression:\n --rate:\t%1.3e +/- %1.2e (one std-dev)\n --chi^2:\t%1.2f\n --r^2:  \t%1.2f\n'
                  %(self.clock_rate, dslope, self.chisq**2, self.r_val**2))
        else:
            outstr = ('Root-Tip-Regression:\n --rate:\t%1.3e\n --r^2:  \t%1.2f\n'
                  %(self.clock_rate, self.r_val**2))

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
        dc.chisq = clock_model['chisq'] if 'chisq' in clock_model else None
        dc.valid_confidence = clock_model['valid_confidence'] if 'valid_confidence' in clock_model else False
        if 'cov' in clock_model and dc.valid_confidence:
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
        res = None

    return res



def parse_dates(date_file):
    """
    parse dates from the arguments and return a dictionary mapping
    taxon names to numerical dates.

    Parameters
    ----------
    date_file : str
        name of file to parse meta data from

    Returns
    -------
    dict
        dictionary linking fields in a column interpreted as taxon name
        (first column that contains 'name', 'strain', 'accession')
        to a numerical date inferred from a column that contains 'date'.
        It will first try to parse the column as float, than via
        pandas.to_datetime and finally as ambiguous date such as 2018-05-XX
    """
    print("\nAttempting to parse dates...")
    dates = {}
    if not os.path.isfile(date_file):
        print("\n\tERROR: file %s does not exist, exiting..."%date_file)
        return dates
    # separator for the csv/tsv file. If csv, we'll strip extra whitespace around ','
    full_sep = '\t' if date_file.endswith('.tsv') else r'\s*,\s*'

    try:
        # read the metadata file into pandas dataframe.
        df = pd.read_csv(date_file, sep=full_sep, engine='python')
        # check the metadata has strain names in the first column
        # look for the column containing sampling dates
        # We assume that the dates might be given either in human-readable format
        # (e.g. ISO dates), or be already converted to the numeric format.
        potential_date_columns = []
        potential_numdate_columns = []
        potential_index_columns = []
        # Scan the dataframe columns and find ones which likely to store the
        # dates
        for ci,col in enumerate(df.columns):
            d = df.iloc[0,ci]
            # strip quotation marks
            if type(d)==str and d[0] in ['"', "'"] and d[-1] in ['"', "'"]:
                for i,tmp_d in enumerate(df.iloc[:,ci]):
                    df.iloc[i,ci] = tmp_d.strip(d[0])
            if 'date' in col.lower():
                potential_date_columns.append((ci, col))
            if any([x==col.lower() for x in ['name', 'strain', 'accession']]):
                potential_index_columns.append((ci, col))

        dates = {}
        # if a potential numeric date column was found, use it
        # (use the first, if there are more than one)
        if not len(potential_index_columns):
            print("ERROR: Cannot read metadata: need at least one column that contains the taxon labels."
                  " Looking for the first column that contains 'name', 'strain', or 'accession' in the header.", file=sys.stderr)
            return dates
        else:
            # use the first column that is either 'name', 'strain', 'accession'
            index_col = sorted(potential_index_columns)[0][1]
            print("\tUsing column '%s' as name. This needs match the taxon names in the tree!!"%index_col)

        if len(potential_date_columns)>=1:
            #try to parse the csv file with dates in the idx column:
            idx = potential_date_columns[0][0]
            col_name = potential_date_columns[0][1]
            print("\tUsing column '%s' as date."%col_name)
            for ri, row in df.iterrows():
                date_str = row.loc[col_name]
                k = row.loc[index_col]
                # try parsing as a float first
                try:
                    dates[k] = float(date_str)
                    continue
                except ValueError:
                    # try whether the date string can be parsed as [2002.2:2004.3]
                    # to indicate general ambiguous ranges
                    if date_str[0]=='[' and date_str[-1]==']' and len(date_str[1:-1].split(':'))==2:
                        try:
                            dates[k] = [float(x) for x in date_str[1:-1].split(':')]
                            continue
                        except ValueError:
                            pass
                    # try date format parsing 2017-08-12
                    try:
                        tmp_date = pd.to_datetime(date_str)
                        dates[k] = numeric_date(tmp_date)
                    except ValueError:  # try ambiguous date format parsing 2017-XX-XX
                        lower, upper = ambiguous_date_to_date_range(date_str, '%Y-%m-%d')
                        if lower is not None:
                            dates[k] = [numeric_date(x) for x in [lower, upper]]

        else:
            print("ERROR: Metadata file has no column which looks like a sampling date!", file=sys.stderr)

        if all(v is None for v in dates.values()):
            print("ERROR: Cannot parse dates correctly! Check date format.", file=sys.stderr)
            return {}
        return dates
    except:
        print("ERROR: Cannot read the metadata file!", file=sys.stderr)
        return {}


def ambiguous_date_to_date_range(mydate, fmt="%Y-%m-%d", min_max_year=None):
    """parse an abiguous date such as 2017-XX-XX to [2017,2017.999]

    Parameters
    ----------
    mydate : str
        date string to be parsed
    fmt : str
        format descriptor. default is %Y-%m-%d
    min_max_year : None, optional
        if date is completely unknown, use this as bounds.

    Returns
    -------
    tuple
        upper and lower bounds on the date. return (None, None) if errors
    """
    from datetime import datetime
    sep = fmt.split('%')[1][-1]
    min_date, max_date = {}, {}
    today = datetime.today().date()

    for val, field  in zip(mydate.split(sep), fmt.split(sep+'%')):
        f = 'year' if 'y' in field.lower() else ('day' if 'd' in field.lower() else 'month')
        if 'XX' in val:
            if f=='year':
                if min_max_year:
                    min_date[f]=min_max_year[0]
                    if len(min_max_year)>1:
                        max_date[f]=min_max_year[1]
                    elif len(min_max_year)==1:
                        max_date[f]=4000 #will be replaced by 'today' below.
                else:
                    return None, None
            elif f=='month':
                min_date[f]=1
                max_date[f]=12
            elif f=='day':
                min_date[f]=1
                max_date[f]=31
        else:
            try:
                min_date[f]=int(val)
                max_date[f]=int(val)
            except ValueError:
                print("Can't parse date string: "+mydate, file=sys.stderr)
                return None, None
    max_date['day'] = min(max_date['day'], 31 if max_date['month'] in [1,3,5,7,8,10,12]
                                           else 28 if max_date['month']==2 else 30)
    lower_bound = datetime(year=min_date['year'], month=min_date['month'], day=min_date['day']).date()
    upper_bound = datetime(year=max_date['year'], month=max_date['month'], day=max_date['day']).date()
    return (lower_bound, upper_bound if upper_bound<today else today)


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

    call = ["iqtree"] + fast_opts +["-nt", str(nthreads), "-s", aln_file, "-m", iqmodel,
            ">", "iqtree.log"]

    os.system(" ".join(call))
    T = Phylo.read(aln_file+".treefile", 'newick')
    for n in T.get_terminals():
        n.name = n.name.replace('_X_X_','/').replace('_Y_Y_','|')
    return T

if __name__ == '__main__':
    pass

