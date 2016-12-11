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
                  %(self.slope, self.intercept, self.r_val**2))
        return outstr


    @classmethod
    def from_tree(cls, t, slope=None):
        """
        Create the conversion object automatically from the tree
        """
        dates = []
        for node in t.find_clades():
            if hasattr(node, "numdate_given") and node.numdate_given is not None:
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
            dc.slope = slope # slope is given
            dc.intercept = np.mean(dates[:,1]) - slope * np.mean(dates[:,0])
            dc.r_val = np.corrcoef(dates[:,1], dates[:,0])[0,1]

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

    def to_years(self, abs_t):
        """
        Convert the time before present measured in branch length units to years

        """
        return abs_t / abs(self.slope)

    def to_numdate(self, tbp):
        """
        Convert the numeric date to the branch-len scale
        """
        return numeric_date() - self.to_years(tbp)



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
