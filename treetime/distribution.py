import numpy as np
from scipy.interpolate import interp1d
from config import WIDTH_DELTA,BIG_NUMBER,MIN_LOG


class Distribution(object):
    """docstring for Distribution"""
    def __init__(self, x, y=1.0, is_log=False, kind='linear'):
        super(Distribution, self).__init__()
        if np.isscalar(x):
            self.kind='delta'
            self.delta_pos = x
            self.weight = y
        else:
            self.kind=kind
            xvals, yvals = np.array(sorted(zip(x,y))).T
            if not is_log:
                yvals = -np.log(yvals)

            # remember range
            self.xmin, self.xmax = xvals[0], xvals[-1]
            self.support = self.xmax - self.xmin

            # extract peak
            self.peak_val = yvals.min()
            yvals -= self.peak_val
            self.ymax = yvals.max()
            self.func= interp1d(xvals, yvals, kind=kind, fill_value=BIG_NUMBER, bounds_error=False)

    def prob(self,x):
        if self.kind=='delta':
            return x==self.delta_pos
        else:
            return np.exp(-self.neg_log(x))

    def neg_log(self,x):
        if self.kind=='delta':
            return (x!=self.delta_pos)*BIG_NUMBER
        else:
            return self.peak_val + self.func(x)

    def x_rescale(self, factor):
        self.func.x*=factor
        self.xmax*=factor
        self.xmax*=factor

def multiply_distributions(dists):
    if any([d.kind=='delta' for d in dists]):
        print("can't multiply delta functions")
    else:
        new_min = np.min([d.xmin for d in dists])
        new_max = np.min([d.xmax for d in dists])
        x_vals = np.concatenate([d.func.x for x in dists])
        x_vals = xvals[(xvals>=new_min)&(xvals<=new_max)]
        y_vals = np.sum([d.neg_log(xvals) for d in dists], axis=0)





