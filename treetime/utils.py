from __future__ import division, print_function
import numpy as np
from scipy.interpolate import interp1d
import config as ttconf
from scipy.integrate import quad
from scipy import stats
import datetime
from scipy.ndimage import binary_dilation

class DateConversion(object):
    """
    Small container class to store parameters to convert between branch length
    as it is used in ML computations and the dates of the nodes.
    It is assumed that the conversion formula is 'length = k*date + b'
    """
    def __init__(self):

        self.slope = 0
        self.intercept = 0
        self.r_val = 0
        self.p_val = 0
        self.sigma = 0

    def __str__(self):
        outstr = ('Root-Tip-Regression:\n --slope:\t%f\n --intercept:\t%f\n --R^2:\t\t%f\n'
                  %(self.slope, self.intercept, self.r_val))
        return outstr


    @classmethod
    def from_tree(cls, t, slope=None):
        """
        Create the conversion object automatically from the tree
        """
        dates = []
        for node in t.find_clades():
            if hasattr(node, "numdate_given" )
                and node.numdate_given is not None.
                dates.append((np.mean(node.numdate_given), node.dist2root))

        if len(dates) == 0:
            raise RuntimeError("Cannot proceed with the TreeTime computations: "
                "No date has been assigned to the terminal nodes!")
        dates = np.array(dates)
        dc = cls()

        if slope is None:
            if len(dates) < 3:
                raise(RuntimeError("There are to few dates set at the leaves of the tree."
                    " Cannot make the conversion function. Aborting."))
            # simple regression
            dc.slope,\
                dc.intercept,\
                dc.r_val,\
                dc.p_val,\
                dc.sigma = stats.linregress(dates[:, 0], dates[:, 1])
        else:
            # TODO this seems awkward
            dc.slope = slope # slope is given
            min_numdate_given = ttconf.BIG_NUMBER
            max_numdate_given = -ttconf.BIG_NUMBER
            max_diam = 0.0
            for node in t.get_terminals():
                # NOTE:  raw_date is time before present in years
                if hasattr(node, 'numdate_given') and node.numdate_given is not None:
                    if node.numdate_given < min_numdate_given:
                        min_numdate_given = node.numdate_given
                    if node.numdate_given > max_numdate_given:
                        max_numdate_given = node.numdate_given
                        max_diam = node.dist2root

            if max_numdate_given == -ttconf.BIG_NUMBER:
                print ("Warning! cannot set the minimal raw date. using today")
                max_numdate_given = 0.0

            if max_diam == 0.0:
                print ("Error! cannot set the intercept for the date2dist conversion!"
                    "Cannot read tree diameter")

            dc.intercept = max_diam - slope * max_numdate_given

        # set the root-mean-square deviation:
        dc.rms = np.sqrt(np.sum((dates[:, 1] - (dc.intercept + dc.slope * dates[:, 0]))**2) / dates.shape[0])
        return dc

    def get_branch_len(self, date1, date2):
        """
        Compute branch length given the dates of the two nodes.

        Args:
         - date1 (int): date of the first node (days before present)
         - date2 (int): date of the second node (days before present)

        Returns:
         - branch length (double): Branch length, assuming that the dependence
         between the node date and the node depth in the the tree is linear.
        """
        return abs(date1 - date2) * self.slope

    def get_time_before_present(self, numdate):
        """
        Convert the numeric date to the branch-len scale
        """
        return (numeric_date() - numdate) * abs(self.slope)

    def get_date(self, abs_t):
        """
        Get the approximate date of the tree node, assuming that the
        dependence between the node date and the node depth int the tree is
        linear.

        Args:
         - node(Phylo.Tree.Clade): node of the tree. Must be from the TreeAnc
         class (or its derivative), to contain the necessary attributes (
            dist2root).

        """
        days = abs_t / abs(self.slope)  #(self.intercept - abs_t) / self.slope
        if days < 0:
            print ("The inferred date of the node is later than today!")
            #print ("Warning: got the negative date! Returning the inverse.")
            #days = abs(days)
        return days


def min_interp(interp_object):
    """
    Find the global minimum of a function represented as an interpolation object.
    """
    try:
        return interp_object.x[interp_object(interp_object.x).argmin()]
    except Exception, e:
        s = "Cannot find minimum of tthe interpolation object" + str(interp_object.x) + \
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
    Args:
     - dt: (datetime.datetime) date of to be converted. if None, assume today
    """
    if dt is None:
        dt = datetime.datetime.now()

    try:
        res = dt.year + dt.timetuple().tm_yday / 365.25
    except:
        res = 0.0

    return res

if __name__ == '__main__':
    pass
