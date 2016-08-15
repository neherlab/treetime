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
        if len(tau)<50:
            tau = np.linspace(tau_min, tau_max, 50)

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
        joint_fwhm  = 0.5 * (node_interp.fwhm+ branch_interp.fwhm)
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
        grid_left =  tmin + (new_peak_pos-tmin) * (1 - np.linspace(1, 0.0, n_grid_points/2)**2.0)
        grid_right = new_peak_pos + (tmax-new_peak_pos)*(np.linspace(0, 1, n_grid_points/2)**2)

        t_grid_0 = np.concatenate([grid_left, grid_right[1:]])
        res_0 = np.array([_convolution_in_point(t_val, node_interp, branch_interp,
                                               n_integral=n_integral, return_log=True,
                                               inverse_time = inverse_time)
                        for t_val in t_grid_0])

        #res_0 = Distribution(t_grid_0, res_0)

        # refine grid
        # TODO estimate error and stop when error is lower than some value
        # (determine the threshold error value)

        refine_idx = np.ones_like(t_grid_0, dtype=bool)
        refine_idx[-1]=False

        step = 0
        while step < 25:

            step += 1

            print ("Refine grid for node convolution, step: " + str(step))
            if refine_idx.sum() < 2:
                break

            if step >= 10:
                import ipdb; ipdb.set_trace()

            # additional points
            add_x = 0.5*(t_grid_0[refine_idx] + t_grid_0[np.roll(refine_idx, 1)])
            add_y = np.array([_convolution_in_point(t_val, node_interp, branch_interp,
                                               n_integral=1000, return_log=True,
                                               inverse_time = inverse_time)
                        for t_val in add_x])

            tan_0 = abs((res_0[np.roll(refine_idx, 1)] - res_0[refine_idx]) / \
                    (t_grid_0[np.roll(refine_idx, 1)] - t_grid_0[refine_idx]))

            tan_1 = abs((add_y - res_0[refine_idx]) / (add_x - t_grid_0[refine_idx]))

            tan_delta = abs((tan_1 - tan_0) / (1 + tan_1 * tan_0))
            insert_point_idx = tan_delta > 1e-3

            plt.plot(add_x, tan_delta, label = str(step))


            # new grid and make new refine idxs:
            n_x = np.concatenate((t_grid_0, add_x[insert_point_idx]))
            n_y = np.concatenate ((res_0, add_y[insert_point_idx]))

            try:
                # For the next refinement, only leave those indexes,
                # where we have inserted additional points
                refine_idx[np.where(refine_idx)] = insert_point_idx
            except:
                import ipdb; ipdb.set_trace()

            n_idx = np.concatenate((refine_idx,
                np.ones_like(add_y[insert_point_idx], dtype=bool)))

            # sort and prepare for the next iteration:
            n_x, n_y, n_idx = np.array(sorted(zip(n_x, n_y, n_idx))).T

            print ("Inserted " + str(insert_point_idx.sum()) + " points to the grid."
                "Points left to refine: " + str(n_idx.sum()) )

            refine_idx = n_idx.astype (bool)
            refine_idx[-1] = False
            t_grid_0 = n_x
            res_0 = n_y

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





