import numpy as np
from distribution import Distribution
import numpy  as np
from config import BIG_NUMBER

def _create_initial_grid(node_dist, branch_dist):
    pass

def _convolution_in_point(node_interp, branch_interp, t_val, n_integral = 100, inverse_time=None, return_log=False):

    if inverse_time is None:
        raise Exception("Inverse time argument must be set!")

    # functions overlap region
    tau_min = np.max((t_val - node_interp.xmax, branch_interp.xmin))
    tau_max = np.min((t_val - node_interp.xmin, branch_interp.xmax))

    if tau_max <= tau_min :
        res = 0.0 #  functions do not overlap

    else:
        # create the tau-grid for the interpolation object in the overlap region
        tau = np.linspace(tau_min, tau_max, n_integral)
        # create the interpolation object on this grid
        if inverse_time:
            fg = node_interp(t_val - tau) + branch_interp(tau)
        else:
            fg = node_interp(t_val + tau) + branch_interp(tau)

        FG = Distribution(tau, fg)
        #integrate the interpolation object
        res = FG.integrate(a=FG.xmin, b=FG.xmax, n=n_integral)

    if not return_log:
        return res

    elif res <= 0.:
        return BIG_NUMBER

    else:
        return np.log(res)

class NodeInterpolator (Distribution):

    @classmethod
    def convolve(cls, node_interp, branch_interp, n_integral=100, inverse_time=True):

        # create coarse grid (5 points)
        import ipdb; ipdb.set_trace()
        joint_fwhm  = 0.5 * (node_interp.fwhm+ branch_interp.fwhm)
        new_peak_pos = node_interp.peak_pos + branch_interp.peak_pos

        initial_times = np.unique([2 * branch_interp.xmin,
                                   new_peak_pos - joint_fwhm,
                                   new_peak_pos,
                                   new_peak_pos + joint_fwhm,
                                   2 * branch_interp.xmax])


        res = np.ones_like(initial_times)

        for t_idx, t_val in enumerate(initial_times):

            res[t_idx] = _convolution_in_point(t_val, node_interp, branch_interp)

        #res = -np.log(res)
        res = cls(initial_times, res, is_log=False)
        # ad
        return res
        # compute the F(t) = integral[f(t-tau) g(tau) dtau] on these 5 points
        # while the curvature is less than the threshold, insert the points in between:

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





