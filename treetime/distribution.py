import numpy as np
from scipy.interpolate import interp1d
from config import WIDTH_DELTA,BIG_NUMBER,MIN_LOG


class Distribution(object):
    """docstring for Distribution"""
    def __init__(self, x, y, is_log=False):
        super(Distribution, self).__init__()
        xvals, yvals = np.array(sorted(zip(x,y))).T
        if not is_log:
            yvals = -np.log(yvals)

        # remember range
        self.xmin, self.xmax = xvals[0], xvals[-1]
        self.support = self.xmax - self.xmin

        # extract peak
        self.peak_val = yvals.min()
        yvals -= self.peak_val

    def make_padded_interpolator(self, x, y):
        tmp_scale = self.support+1
        dx = self.support*WIDTH_DELTA
        tmp_x = np.concatenate(([self.xmin-tmp_scale*BIG_NUMBER, self.xmin-dx],
                                x, [self.xmin-dx, self.xmax+tmp_scale*BIG_NUMBER]))
        tmp_y = np.concatenate(([-MIN_LOG, -MIN_LOG], y, [-MIN_LOG, -MIN_LOG])
        self.func= interp1d(tmp_x, tmp_y, kind='linear')

    def prob(self,x):
        return np.exp(-self.neg_log(x))

    def neg_log(self,x):
        return self.peak_val + self.func(x)

