from __future__ import division, print_function
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


    def integrate_trapez(self, a, b,n):
        mult = 0.5
        if a>b:
            b,a = a,b
            mult=-0.5

        x = np.linspace(a,b,n)
        dx = np.diff(x)
        y = self.prob(x)
        return mult*np.sum(dx*(y[:-1] + y[1:]))


    def integrate_simpson(self, a, b, n):
        mult = 1.0/6
        if a>b:
            b,a = a,b
            mult=-1.0/6
        x = np.linspace(a,b,n)
        dx = np.diff(x[::2])
        y = self.prob(x)
        print(x.shape, dx.shape)
        return mult*(dx[0]*y[0]+ np.sum(4*dx*y[1:-1:2]) + np.sum((dx[:-1]+dx[1:])*y[2:-1:2]) + dx[-1]*y[-1])


    def integrate_piecewise(self, a, b):
        mult = 1.0
        if a>b:
            b,a = a,b
            mult=-1.0
        ind = (self.func.x>a)&(self.func.x<b)
        x = np.concatenate(([a], self.func.x[ind], [b]))
        dx = np.diff(x)
        y = self.neg_log(x)
        yexp = self.prob(x)
        slope = -np.diff(y)/dx

        return mult*np.sum(yexp[:-1]*np.abs((1-np.exp(slope*dx))/slope))

    def integrate_piecewise_gaussian(self, a, b):
        from scipy.special import erf
        from scipy.interpolate import spleval
        mult = 0.5*np.sqrt(np.pi)
        if a>b:
            b,a = a,b
            mult=-mult
        ind = (self.func.x>a)&(self.func.x<b)
        x = np.concatenate(([a], self.func.x[ind], [b]))
        dx = np.diff(x)
        y = self.neg_log(x)
        yexp = self.prob(x)
        cslope =    np.squeeze(spleval(self.func._spline,x,deriv=1))[1:-1]
        curvature = 0.5*np.squeeze(spleval(self.func._spline,x,deriv=2))[1:-1]
        slope=cslope
        #slope = np.diff(y)/dx
        #curvature = np.diff(slope)/(dx[:1]+dx[1:])*2
        #cslope = (y[2:]-y[:-2])/(dx[:1]+dx[1:])*2

        y0 = yexp[1:-1]*np.exp(cslope**2/4.0/curvature)
        t1 = np.sqrt(curvature)*(cslope/2.0/curvature-dx[:-1]*0.5)
        t2 = np.sqrt(curvature)*(cslope/2.0/curvature +dx[1:]*0.5)
        ind = (curvature>0)&(np.abs(cslope**2/4.0/curvature)<200)
        #print(y,'X', t1[ind], t2[ind])
        #print(curvature[ind], slope[ind])
        i1 = np.sum(y0[ind]*(erf(t2[ind]) - erf(t1[ind]))/np.sqrt(curvature[ind]))

        return mult*i1



class quadratic_interpolator(object):
    """docstring for quadratic_interpolator"""
    def __init__(self, x,y):
        super(quadratic_interpolator, self).__init__()
        if len(x)%2==0:
            x = np.concatenate((x[:2],[x[-1]]))
            y = np.concatenate((y[:2],[y[-1]]))
        self.x = x[::2]
        self.xa = x[1::2]
        self.y = y[::2]
        self.ya = y[1::2]
        self.set_up()

    def set_up(self):
        x1 = self.x[:-1]
        x2 = self.xa
        x3 = self.x[1:]
        x1sq = x1**2
        x2sq = x2**2
        x3sq = x3**2
        dx12 = x1 - x2
        dx23 = x2 - x3
        y1 = self.y[:-1]
        y2 = self.ya
        y3 = self.y[1:]
        dy12 = y1 - y2
        dy23 = y2 - y3

        a = (dy12/dx12 - dy23/dx23)/((x1sq-x2sq)/dx12 - (x2sq-x3sq)/dx23)
        b = (dy12 - a*(x1sq-x2sq))/dx12
        c = y1 - a*x1sq - b*x1

        self.a = a
        self.b = b
        self.c = c

    def __call__(self, x):
        ii = np.minimum(np.maximum(self.x.searchsorted(x)-1,0), self.a.shape[0]-1)
        return self.a[ii]*x**2 + self.b[ii]*x + self.c[ii]



def multiply_distributions(dists):
    if any([d.kind=='delta' for d in dists]):
        print("can't multiply delta functions")
        return None
    else:
        new_xmin = np.max([d.xmin for d in dists])
        new_xmax = np.min([d.xmax for d in dists])
        x_vals = np.concatenate([d.func.x for x in dists])
        x_vals = xvals[(xvals>=new_xmin)&(xvals<=new_xmax)]
        y_vals = np.sum([d.neg_log(xvals) for d in dists], axis=0)
        return Distribution(x_vals, y_vals, is_log=True, kind='linear')


if __name__=="__main__":
    from matplotlib import pyplot as plt
    plt.ion()

    def f(x):
        #return (x-5)**2+np.abs(x)**3
        return (x**2-5)**2 #(x-5)**2+np.abs(x)**3
    def g(x):
        return (x-4)**2*(x**(1.0/3)-5)**2

    plot=False
    error = {}
    for kind in ['linear', 'quadratic', 'cubic', 'Q']:
        error[kind]=[[],[]]
    npoints = [11,21] #,31,41, 51,75,101]
    for ex, func in [[0,f],[1,g]]:
        for npoint in npoints:
            if ex==0:
                xnew = np.linspace(-5,15,1000)
                x = np.linspace(0,10,npoint)
            elif ex==1:
                xnew = np.linspace(0,150,1000)
                x = np.linspace(0,9.0,npoint)**3

            if plot:
                plt.figure()
                plt.plot(x, np.exp(-func(x)),'-o', label = 'data')
                plt.plot(xnew, np.exp(-func(xnew)),'-',label='true')
            for kind in ['linear', 'quadratic', 'cubic']:
                try:
                    dist = Distribution(x, func(x), kind=kind, is_log=True)
                    if plot: plt.plot(xnew, dist.prob(xnew), label=kind)
                    E = np.mean((np.exp(-func(xnew))-dist.prob(xnew))[(xnew>dist.xmin) & (xnew<dist.xmax)]**2)
                    print(kind,npoint, E)
                except:
                    E=np.nan
                error[kind][ex].append(E)
            try:
                distQ = quadratic_interpolator(x, func(x))
                if plot: plt.plot(xnew[(xnew>dist.xmin) & (xnew<dist.xmax)], np.exp(-distQ(xnew))[(xnew>dist.xmin) & (xnew<dist.xmax)], label='Q')
                E = np.mean((np.exp(-func(xnew))-np.exp(-distQ(xnew)))[(xnew>dist.xmin) & (xnew<dist.xmax)]**2)
                print('Q',npoint, E)
            except:
                E=np.nan
            error['Q'][ex].append(E)
            if plot:
                plt.yscale('log')
                plt.legend()

    for ex in [0,1]:
        plt.figure()
        for k in error:
            plt.plot(npoints, error[k][ex],'-o', label=k)
        plt.yscale('log')
        plt.legend()

    integration_error = {'trapez':[], 'simpson':[], 'piecewise':[]}
    npoints = [11,21,31, 41, 51,75,101, 201, 501, 1001]
    xnew = np.linspace(-5,15,1000)
    x = np.linspace(0,10,3000)
    dist = Distribution(x, g(x), kind='linear', is_log=True)
    for npoint in npoints:
        integration_error['trapez'].append(dist.integrate_trapez(0,10,npoint))
        integration_error['simpson'].append(dist.integrate_simpson(0,10,npoint))
        xtmp = np.linspace(0,10,min(100,npoint))
        disttmp = Distribution(xtmp, g(xtmp), kind='cubic', is_log=True)
        integration_error['piecewise'].append(disttmp.integrate_piecewise_gaussian(0,10))

    plt.figure()
    base_line = integration_error['simpson'][-1]
    plt.plot(npoints, np.abs(integration_error['trapez']-base_line), label='trapez')
    plt.plot(npoints, np.abs(integration_error['simpson']-base_line), label='simpson')
    plt.plot(npoints, np.abs(integration_error['piecewise']-base_line), label='piecewise')
    plt.xlabel('npoints')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
