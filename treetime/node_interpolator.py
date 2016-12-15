import numpy as np
from distribution import Distribution
import numpy  as np
import config as ttconf

def _create_initial_grid(node_dist, branch_dist):
    pass

def _convolution_integrand(t_val, f, g, inverse_time=None, return_log=False):
    '''
    evaluates int_tau f(t+tau)g(tau) or int_tau f(t-tau)g(tau) if inverse time is TRUE
    '''
    if inverse_time is None:
        raise Exception("Inverse time argument must be set!")

    # determine integration boundaries:
    if inverse_time:
        ## tau>g.xmin and t-tau<f.xmax
        tau_min = max(t_val - f.xmax, g.xmin)
        ## tau<g.xmax and t-tau>f.xmin
        tau_max = min(t_val - f.xmin, g.xmax)
    else:
        ## tau>g.xmin and t+tau>f.xmin
        tau_min = max(f.xmin-t_val, g.xmin)
        ## tau<g.xmax and t+tau<f.xmax
        tau_max = min(f.xmax-t_val, g.xmax)
        #print(tau_min, tau_max)


    if tau_max <= tau_min + ttconf.TINY_NUMBER:
        if return_log:
            return ttconf.BIG_NUMBER
        else:
            return 0.0 #  functions do not overlap

    else:
        # create the tau-grid for the interpolation object in the overlap region
        if inverse_time:
            tau = np.unique(np.concatenate((g.x, t_val-f.x,[tau_min,tau_max])))
        else:
            tau = np.unique(np.concatenate((g.x, f.x-t_val,[tau_min,tau_max])))
        tau = tau[(tau>tau_min-ttconf.TINY_NUMBER)&(tau<tau_max+ttconf.TINY_NUMBER)]
        if len(tau)<10:
            tau = np.linspace(tau_min, tau_max, 10)

        if inverse_time: # add negative logarithms
            tnode = t_val - tau
            fg = f(tnode) + g(tau, tnode=tnode)
        else:
            fg = f(t_val + tau) + g(tau, tnode=t_val)

        # create the interpolation object on this grid
        FG = Distribution(tau, fg, is_log=True, kind='linear')
        return FG



def _max_of_integrand(t_val, f, g, inverse_time=None, return_log=False):

    # return log is always True
    FG = _convolution_integrand(t_val, f, g, inverse_time, return_log=True)

    if FG == ttconf.BIG_NUMBER:
        res = ttconf.BIG_NUMBER, 0

    else:

        # TODO use minimize_scalar here !!
        X = FG.x[FG.y.argmin()]
        Y = FG.y.min()
        res =  Y, X

        #import matplotlib.pyplot as plt
        #plt.figure(11)
        #plt.plot(FG.x, FG.y)
        #plt.vlines(res[1], -1e3, Y)

    if not return_log:
        res[0] = np.log(res[0])


    return res

def _evaluate_convolution(t_val, f, g,  n_integral = 100, inverse_time=None, return_log=False):


    FG = _convolution_integrand(t_val, f, g, inverse_time, return_log)

    #integrate the interpolation object, return log, make neg_log
        #print('FG:',FG.xmin, FG.xmax, FG(FG.xmin), FG(FG.xmax))
    if (return_log and FG == ttconf.BIG_NUMBER) or \
        (not return_log and FG == 0.0): # distributions do not overlap
        res = ttconf.BIG_NUMBER # we integrate log funcitons
    else:
        res = -FG.integrate(a=FG.xmin, b=FG.xmax, n=n_integral, return_log=True)

    if return_log:
        return res, -1
    else:
        return np.exp(-res), -1


class NodeInterpolator (Distribution):

    @classmethod
    def convolve(cls, node_interp, branch_interp, max_or_integral='integral',
                 n_integral=100, inverse_time=True, rel_tol=0.05, yc=10):

        '''
        calculate H(t) = \int_tau f(t-tau)g(tau) if inverse_time=True
                  H(t) = \int_tau f(t+tau)g(tau) if inverse_time=False

        This function determines the time points of the grid of the result to
        ensure an accurate approximation.
        '''

        if max_or_integral not in ['max', 'integral']:
            raise Exception("Max_or_integral expected to be 'max' or 'integral', got "
                            + str(max_or_integral)  + " instead.")

        def conv_in_point(time_point):

            if max_or_integral == 'integral': # compute integral of the convolution
                return _evaluate_convolution(time_point, node_interp, branch_interp,
                                               n_integral=n_integral, return_log=True,
                                               inverse_time = inverse_time)

            else: # compute max of the convolution
                return _max_of_integrand(time_point, node_interp, branch_interp,
                                               return_log=True, inverse_time = inverse_time)

        # estimate peak and width
        joint_fwhm  = (node_interp.fwhm + branch_interp.fwhm)
        min_fwhm  = min(node_interp.fwhm, branch_interp.fwhm)
        # determine support of the resulting convolution
        # in order to be positive, the flipped support of f, shifted by t and g need to overlap
        if inverse_time:
            new_peak_pos = node_interp.peak_pos + branch_interp.peak_pos
            tmin = node_interp.xmin+branch_interp.xmin
            tmax = node_interp.xmax+branch_interp.xmax
        else:
            new_peak_pos = node_interp.peak_pos - branch_interp.peak_pos
            tmin = node_interp.xmin - branch_interp.xmax
            tmax = node_interp.xmax - branch_interp.xmin

        # make initial node grid consisting of linearly spaced points around
        # the center and quadratically spaced points at either end
        n_grid_points = ttconf.NODE_GRID_SIZE
        n = n_grid_points/3
        center_width = 3*joint_fwhm
        grid_center = new_peak_pos + np.linspace(-1, 1, n)*center_width

        # add the right and left grid if it is needed
        right_range = (tmax - grid_center[-1])
        if right_range>4*center_width:
            grid_right = grid_center[-1] + right_range*(np.linspace(0, 1, n)**2.0)
        elif right_range>0: # use linear grid the right_range is comparable to center_width
            grid_right = grid_center[-1] + right_range*np.linspace(0,1, int(min(n,1+0.5*n*right_range/center_width)))
        else:
            grid_right =[]

        left_range = grid_center[0]-tmin
        if left_range>4*center_width:
            grid_left = tmin + left_range*(np.linspace(0, 1, n)**2.0)
        elif left_range>0:
            grid_left = tmin + left_range*np.linspace(0,1, int(min(n,1+0.5*n*left_range/center_width)))
        else:
            grid_left =[]

        # make grid and calculate convolution
        t_grid_0 = np.unique(np.concatenate([grid_left[:-1], grid_center, grid_right[1:], [tmin, tmax]]))
        t_grid_0 = t_grid_0[(t_grid_0 > tmin-ttconf.TINY_NUMBER) & (t_grid_0 < tmax+ttconf.TINY_NUMBER)]

        # res0 - the values of the convolution (integral or max)
        # t_0 - the value, at which the res0 achieves maximum (for 'max' type  of convolution)
        # for 'itegral' type t_0 is -1, without any meaning
        res_0, t_0 = np.array([conv_in_point(t_val) for t_val in t_grid_0]).T

        # refine grid as necessary and add new points
        # calculate interpolation error at all internal points [2:-2] bc end points are sometime off scale
        interp_error = np.abs(res_0[3:-1]+res_0[1:-3]-2*res_0[2:-2])
        # determine the number of extra points needed, criterion depends on distance from peak dy
        dy = (res_0[2:-2]-res_0.min())
        dx = np.diff(t_grid_0)
        refine_factor = np.minimum(np.minimum(np.array(np.floor(np.sqrt(interp_error/(rel_tol*(1+(dy/yc)**4)))), dtype=int),
                                   np.array(100*(dx[1:-2]+dx[2:-1])/min_fwhm, dtype=int)), 10)

        insert_point_idx = np.zeros(interp_error.shape[0]+1, dtype=int)
        insert_point_idx[1:] = refine_factor
        insert_point_idx[:-1] += refine_factor
        # add additional points if there are any to add

        if np.sum(insert_point_idx):
            add_x = np.concatenate([np.linspace(t1,t2,n+2)[1:-1] for t1,t2,n in
                               zip(t_grid_0[1:-2], t_grid_0[2:-1], insert_point_idx) if n>0])
            # calculate convolution at these points
            add_y, add_t = np.array([conv_in_point(t_val) for t_val in add_x]).T

            t_grid_0 = np.concatenate((t_grid_0, add_x))
            res_0 = np.concatenate ((res_0, add_y))
            t_0 = np.concatenate ((t_0, add_t))

        # instantiate the new interpolation object and return
        res_y = cls(t_grid_0, res_0, is_log=True, kind='linear')

        # the interpolation object,  which is used to store the value of the
        # grid, which maximizes the convolution (for 'max' option),
        # or flat -1 distribution (for 'integral' option)
        res_t = Distribution(t_grid_0, t_0, is_log=True, kind='linear')

        # the convolution and the values of the 'tau' grid, which maximize the
        # convolution values in the 'max' regime

        #import matplotlib.pyplot as plt
        #plt.figure(1)
        #plt.plot(res_y.x, res_y.y)
        #plt.figure(2)
        #plt.plot(res_y.x, t_0)
        #import ipdb; ipdb.set_trace()
        return res_y, res_t

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    plt.ion()
    plt.show()

    print ("NodeInterpolator test")

    for n in [2, 5, 10, 100, 500]:

        xl = np.linspace(0,1,n)
        xr = np.linspace(1, 2, n)
        yl = 1e-10 + 10*xl
        yr = 1e-10 + -10*xr + 20

        x = np.concatenate([xl, xr])
        y = np.concatenate([yl, yr])

        d1 = Distribution(x, y,is_log=False)
        d2 = Distribution(x, y,is_log=False)

        ni = NodeInterpolator.convolve(d1, d2, n_integral=1000)

        plt.figure(1)
        plt.plot(d1.x, np.exp(-d1.y), 'o--', label="# tau points: " + str(n))
        plt.figure(2)
        plt.plot(ni.x, ni.y, 'o-' , label="# tau points: " + str(n))
        plt.figure(3)
        plt.plot(ni.x, np.exp(-ni.y), 'o-', label="# tau points: " + str(n))




    plt.figure(1)
    plt.legend()
    plt.figure(2)
    plt.legend()
    plt.figure(3)
    plt.legend()





