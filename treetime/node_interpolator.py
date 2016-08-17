import numpy as np
from distribution import Distribution
import numpy  as np
import config as ttconf

def _create_initial_grid(node_dist, branch_dist):
    pass

def _convolution_in_point(t_val,f, g,  n_integral = 100, inverse_time=None, return_log=False):
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


    if tau_max <= tau_min :
        if return_log:
            return ttconf.BIG_NUMBER
        else:
            return 0.0 #  functions do not overlap

    else:
        # create the tau-grid for the interpolation object in the overlap region
        # TODO: make clever grid

        # create the interpolation object on this grid
        if inverse_time: # add negative logarithms
            tau = np.unique(np.concatenate((g.x, t_val-f.x)))
        else:
            tau = np.unique(np.concatenate((g.x, f.x-t_val)))
        tau = tau[(tau>=tau_min)&(tau<tau_max)]
        if len(tau)<200:
            tau = np.linspace(tau_min, tau_max, 150)

        if inverse_time: # add negative logarithms
            fg = f(t_val - tau) + g(tau)
        else:
            fg = f(t_val + tau) + g(tau)

        # TODO: break into segments: peak and tails
        FG = Distribution(tau, fg, is_log=True, kind='linear')
        #integrate the interpolation object, return log, make neg_log
        res = -FG.integrate(a=FG.xmin, b=FG.xmax, n=n_integral, return_log=True)

        if return_log:
            return res
        else:
            return np.exp(-res)

class NodeInterpolator (Distribution):

    @classmethod
    def convolve(cls, node_interp, branch_interp, n_integral=1000, inverse_time=True):

        '''
        calculate H(t) = \int_tau f(t-tau)g(tau) if inverse_time=True
                  H(t) = \int_tau f(t+tau)g(tau) if inverse_time=False

        This function determines the time points of the grid of the result to
        ensure an accurate approximation.
        '''

        import matplotlib.pyplot as plt

        # create coarse grid (5 points)
        joint_fwhm  = (node_interp.fwhm + branch_interp.fwhm)
        new_peak_pos = node_interp.peak_pos + branch_interp.peak_pos

        # determine support of the resulting convolution
        # in order to be positive, the flipped support of f, shifted by t and g need to overlap
        if inverse_time:
            tmin = node_interp.xmin+branch_interp.xmin
            tmax = node_interp.xmax+branch_interp.xmax
        else:
            tmin = node_interp.xmin - branch_interp.xmax
            tmax = node_interp.xmax - branch_interp.xmin

        # make initial node grid
        n_grid_points = ttconf.NODE_GRID_SIZE
        n = n_grid_points/3
        center_width = 3*joint_fwhm
        grid_center = new_peak_pos + np.linspace(-1, 1, n)*center_width

        dx = grid_center[1]-grid_center[0]
        right_range = (tmax - grid_center[-1])

        if right_range>center_width:
            grid_right = grid_center[-1] + right_range*(np.linspace(0, 1, n)**2.0)
        elif right_range>0:
            grid_right = grid_center[-1] + right_range*np.linspace(0,1,n)
        else:
            grid_right =[]

        left_range = grid_center[0]-tmin
        if left_range>center_width:
            grid_left = tmin - left_range*(np.linspace(0, 1, n)**2.0)
        elif left_range>0:
            grid_left = tmin + left_range*np.linspace(0,1,n)
        else:
            grid_left =[]

        t_grid_0 = np.concatenate([grid_left[:-1], grid_center, grid_right[1:]])
        t_grid_0 = t_grid_0[(t_grid_0 > tmin) & (t_grid_0 < tmax)]
        res_0 = np.array([_convolution_in_point(t_val, node_interp, branch_interp,
                                               n_integral=n_integral, return_log=True,
                                               inverse_time = inverse_time)
                        for t_val in t_grid_0])

        # (determine the threshold error value)
        step = 0
        print ("Refine grid for node convolution, step: " + str(step))
        while step < 1:
            step += 1
            interp_error = np.abs(res_0[3:-1]+res_0[1:-3]-2*res_0[2:-2])
            dy = (res_0[2:-2]-res_0.min())
            yc = 10
            rel_tol = 0.01
            refine_factor = np.array(np.floor(np.sqrt(interp_error/(rel_tol*(1+(dy/yc)**4)))), dtype=int)
            insert_point_idx = np.zeros(interp_error.shape[0]+1, dtype=int)
            insert_point_idx[1:] = refine_factor
            insert_point_idx[:-1] += refine_factor
            if np.sum(insert_point_idx):
                add_x = np.concatenate([np.linspace(t1,t2,n+2)[1:-1] for t1,t2,n in
                                   zip(t_grid_0[1:-2], t_grid_0[2:-1], insert_point_idx) if n>0])
            else:
                break
            # additional points
            add_y = np.array([_convolution_in_point(t_val, node_interp, branch_interp,
                                                    n_integral=1000, return_log=True,
                                                    inverse_time = inverse_time)
                             for t_val in add_x])

            n_x = np.concatenate((t_grid_0, add_x)) #[insert_point_idx]))
            n_y = np.concatenate ((res_0, add_y)) #[insert_point_idx]))
            n_x, n_y = np.array(sorted(zip(n_x, n_y))).T
            t_grid_0 = n_x
            res_0 = n_y
        print ("final grid: " , len(res_0))

        res = cls(t_grid_0, res_0, is_log=True, kind='linear')

        return res

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





