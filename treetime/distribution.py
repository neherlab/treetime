from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.interpolate import interp1d
try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable
from copy import deepcopy as make_copy
from scipy.ndimage import binary_dilation
from .config import BIG_NUMBER, MIN_LOG, MIN_INTEGRATION_PEAK, TINY_NUMBER
from .utils import clip

class Distribution(object):
    """
    Class to implement the probability distribution. This class wraps the scipy
    linear interpolation object, and implements some additional operations,
    needed to manipulate distributions for tree nodes positions, branch lengths,
    etc.
    This class is callable, so it can be treated similarly to the scipy interpolation
    object.
    """

    @staticmethod
    def calc_fwhm(distribution, is_neg_log=True):
        """
        Assess the width of the probability distribution. This returns
        full-width-half-max
        """

        if isinstance(distribution, interp1d):

            if is_neg_log:
                ymin = distribution.y.min()
                log_prob = distribution.y-ymin
            else:
                log_prob = -np.log(distribution.y)
                log_prob -= log_prob.min()

            xvals = distribution.x

        elif isinstance(distribution, Distribution):
            # Distribution always stores neg log-prob with the peak value subtracted
            xvals = distribution._func.x
            log_prob = distribution._func.y
        else:
            raise TypeError("Error in computing the FWHM for the distribution. "
                " The input should be either Distribution or interpolation object")

        L = xvals.shape[0]
        # 0.69... is log(2), there is always one value for which this is true since
        # the minimum is subtracted
        tmp = np.where(log_prob < 0.693147)[0]
        x_l, x_u = tmp[0], tmp[-1]
        if L < 2:
            print ("Not enough points to compute FWHM: returning zero")
            return min(TINY_NUMBER, distribution.xmax - distribution.xmin)
        else:
            # need to guard against out-of-bounds errors
            return max(TINY_NUMBER, xvals[min(x_u+1,L-1)] - xvals[max(0,x_l-1)])


    @classmethod
    def delta_function(cls, x_pos, weight=1., min_width=MIN_INTEGRATION_PEAK):
        """
        Create delta function distribution.
        """

        distribution = cls(x_pos,0.,is_log=True, min_width=min_width)
        distribution.weight  = weight
        return distribution


    @classmethod
    def shifted_x(cls, dist, delta_x):
        return Distribution(dist.x+delta_x, dist.y, kind=dist.kind)


    @staticmethod
    def multiply(dists):
        '''
        multiplies a list of Distribution objects
        '''
        if  not all([isinstance(k, Distribution) for k in dists]):
            raise NotImplementedError("Can only multiply Distribution objects")

        n_delta = np.sum([k.is_delta for k in dists])
        min_width = np.max([k.min_width for k in dists])
        if n_delta>1:
            raise ArithmeticError("Cannot multiply more than one delta functions!")
        elif n_delta==1:
            delta_dist_ii = np.where([k.is_delta for k in dists])[0][0]
            delta_dist = dists[delta_dist_ii]
            new_xpos = delta_dist.peak_pos
            new_weight  = np.prod([k.prob(new_xpos) for k in dists if k!=delta_dist_ii]) * delta_dist.weight
            res = Distribution.delta_function(new_xpos, weight = new_weight,min_width=min_width)
        else:
            new_xmin = np.max([k.xmin for k in dists])
            new_xmax = np.min([k.xmax for k in dists])

            x_vals = np.unique(np.concatenate([k.x for k in dists]))
            x_vals = x_vals[(x_vals>new_xmin-TINY_NUMBER)&(x_vals<new_xmax+TINY_NUMBER)]
            y_vals = np.sum([k.__call__(x_vals) for k in dists], axis=0)
            peak = y_vals.min()
            ind = (y_vals-peak)<BIG_NUMBER/1000
            n_points = ind.sum()
            if n_points == 0:
                print ("ERROR in distribution multiplication: Distributions do not overlap")
                x_vals = [0,1]
                y_vals = [BIG_NUMBER,BIG_NUMBER]
                res = Distribution(x_vals, y_vals, is_log=True,
                                   min_width=min_width, kind='linear')
            elif n_points == 1:
                res = Distribution.delta_function(x_vals[0])
            else:
                res = Distribution(x_vals[ind], y_vals[ind], is_log=True,
                                   min_width=min_width, kind='linear', assume_sorted=True)

        return res


    def __init__(self, x, y, is_log=True, min_width = MIN_INTEGRATION_PEAK,
                 kind='linear', assume_sorted=False):

        """
        Create Distribution instance
        """

        self.min_width = min_width
        if isinstance(x, Iterable) and isinstance (y, Iterable):

            self._delta = False # NOTE in classmethod this value is set explicitly to True.
            # first, prepare x, y values
            if assume_sorted:
                xvals, yvals = x,y
            else:
                xvals, yvals = np.array(sorted(zip(x,y))).T
            if not is_log:
                yvals = -np.log(yvals)
            # just for safety
            yvals[np.isnan(yvals)] = BIG_NUMBER
            # set the properties
            self._kind=kind
            # remember range
            self._xmin, self._xmax = xvals[0], xvals[-1]
            self._support = self._xmax - self._xmin
            # extract peak
            self._peak_idx = yvals.argmin()
            self._peak_val = yvals.min()
            self._peak_pos = xvals[self._peak_idx]
            yvals -= self._peak_val
            self._ymax = yvals.max()
            # store the interpolation object
            self._func= interp1d(xvals, yvals, kind=kind, fill_value=BIG_NUMBER,
                                 bounds_error=False, assume_sorted=True)
            self._fwhm = Distribution.calc_fwhm(self)

        elif np.isscalar(x):
            assert (np.isscalar(y) or y is None)
            self._delta = True
            self._peak_pos = x
            self._fwhm = 0
            if y is None:
                self._peak_val = np.inf
            else:
                self._peak_val = y

            self._xmin, self._xmax = x, x
            self._support = 0.
            self._func = lambda x : (x==self.peak_pos)*self.peak_val
        else:
            raise TypeError("Cannot create Distribution: "
                "Input arguments should be scalars or iterables!")


    @property
    def is_delta(self):
        return self._delta

    @property
    def kind(self):
        return self._kind

    @property
    def peak_val(self):
        return self._peak_val

    @property
    def peak_pos(self):
        return self._peak_pos

    @property
    def peak_idx(self):
        return self._peak_idx

    @property
    def support(self):
        return self._support

    @property
    def fwhm(self):
        return self._fwhm

    @property
    def x(self):
        if self.is_delta:
            return [self._peak_pos]
        else:
            return self._func.x

    @property
    def y(self):
        if self.is_delta:
            print("THIS SHOULDN'T BE CALLED ON A DELTA FUNCTION")
            return [self.weight]
        else:
            return self._peak_val + self._func.y

    @property
    def xmin(self):
        return self._xmin

    @property
    def xmax(self):
        return self._xmax


    def __call__(self, x):
        if isinstance(x, Iterable):
            valid_idxs = (x > self._xmin-TINY_NUMBER) & (x < self._xmax+TINY_NUMBER)
            res = np.full(np.shape(x), BIG_NUMBER+self.peak_val, dtype=float)
            tmp_x = x[valid_idxs]
            res[valid_idxs] = self._peak_val + self._func(clip(tmp_x, self._xmin+TINY_NUMBER, self._xmax-TINY_NUMBER))
            return res

        elif np.isreal(x):
            if x < self._xmin or x > self._xmax:
                return BIG_NUMBER+self.peak_val
            # x is within interpolation range
            elif self._delta == True:
                return self._peak_val
            else:
                return self._peak_val + self._func(x)
        else:
            raise TypeError("Wrong type: should be float or array")


    def __mul__(self, other):
        return Distribution.multiply((self, other))


    def _adjust_grid(self, rel_tol=0.01, yc=10):
        updated = True
        n_iter=0
        while len(self.y)>200 and updated and n_iter<5:
            interp_err = 2*self.y[1:-1] - self.y[2:] - self.y[:-2]
            ind = np.ones_like(self.y, dtype=bool)
            dy = self.y-self.peak_val
            prune = interp_err[::2] > rel_tol*(1+ (dy[1:-1:2]/yc)**4)
            ind[1:-1:2] = prune
            if np.mean(prune)<1.0:
                self._func.y = self._func.y[ind]
                self._func.x = self._func.x[ind]
                updated=True
                n_iter+=1
            else:
                updated=False
                n_iter+=1

        self._peak_idx = self.__call__(self._func.x).argmin()
        self._peak_pos = self._func.x[self._peak_idx]
        self._peak_val = self.__call__(self.peak_pos)


    def prob(self,x):
        return np.exp(-1 * self.__call__(x))

    def prob_relative(self,x):
        return np.exp(-1 * (self.__call__(x)-self.peak_val))

    def x_rescale(self, factor):
        self._func.x*=factor
        self._peak_pos*=factor
        if factor>=0:
            self._xmin*=factor
            self._xmax*=factor
        else:
            tmp = self.xmin
            self._xmin = factor*self.xmax
            self._xmax = factor*tmp
            self._func.x = self._func.x[::-1]
            self._func.y = self._func.y[::-1]


    def integrate(self, return_log=False ,**kwargs):
        if self.is_delta:
            return self.weight
        else:
            integral_result = self.integrate_simpson(**kwargs)
            if return_log:
                if integral_result==0:
                    return -self.peak_val - BIG_NUMBER
                else:
                    return -self.peak_val + max(-BIG_NUMBER, np.log(integral_result))
            else:
                return np.exp(-self.peak_val)*integral_result

    def integrate_trapez(self, a=None, b=None,n=None):
        mult = 0.5
        if a>b:
            b,a = a,b
            mult=-0.5

        x = np.linspace(a,b,n)
        dx = np.diff(x)
        y = self.prob_relative(x)
        return mult*np.sum(dx*(y[:-1] + y[1:]))


    def integrate_simpson(self, a=None,b=None,n=None):
        if n % 2 == 0:
            n += 1
        mult = 1.0/6
        dpeak = max(10*self.fwhm, self.min_width)
        threshold = np.array([a,self.peak_pos-dpeak, self.peak_pos+dpeak,b])
        threshold = threshold[(threshold>=a)&(threshold<=b)]
        threshold.sort()
        res = []
        for lw, up in zip(threshold[:-1], threshold[1:]):
            x = np.linspace(lw,up,n)
            dx = np.diff(x[::2])
            y = self.prob_relative(x)
            res.append(mult*(dx[0]*y[0]+ np.sum(4*dx*y[1:-1:2])
                    + np.sum((dx[:-1]+dx[1:])*y[2:-1:2]) + dx[-1]*y[-1]))

        return np.sum(res)

