import numpy as np
from . import config as ttconf
from .distribution import Distribution
from .utils import clip
from .config import FFT_FWHM_GRID_SIZE

def _convolution_integrand(t_val, f, g,
                           inverse_time=None, return_log=False):
    '''
    Evaluates int_tau f(t+tau)*g(tau) or int_tau f(t-tau)g(tau) if inverse time is TRUE

    Parameters
    -----------

     t_val : double
        Time point

     f : Interpolation object
        First multiplier in convolution

     g : Interpolation object
        Second multiplier in convolution

     inverse_time : bool, None
        time direction. If True, then the f(t-tau)*g(tau) is calculated, otherwise,
        f(t+tau)*g(tau)

     return_log : bool
        If True, the logarithm will be returned


    Returns
    -------

     FG : Distribution
        The function to be integrated as Distribution object (interpolator)

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


    if tau_max <= tau_min:
        if return_log:
            return ttconf.BIG_NUMBER
        else:
            return 0.0 #  functions do not overlap

    else:
        # create the tau-grid for the interpolation object in the overlap region
        if inverse_time:
            tau = np.concatenate((g.x, t_val-f.x,[tau_min,tau_max]))
        else:
            tau = np.concatenate((g.x, f.x-t_val,[tau_min,tau_max]))
        tau = np.unique(clip(tau, tau_min-ttconf.TINY_NUMBER, tau_max+ttconf.TINY_NUMBER))
        if len(tau)<10:
            tau = np.linspace(tau_min, tau_max, 10)

        if inverse_time: # add negative logarithms
            tnode = t_val - tau
            fg = f(tnode) + g(tau)
        else:
            fg = f(t_val + tau) + g(tau)

        # create the interpolation object on this grid
        FG = Distribution(tau, fg, is_log=True, min_width = np.max([f.min_width, g.min_width]),
                          kind='linear', assume_sorted=True)
        return FG



def _max_of_integrand(t_val, f, g, inverse_time=None, return_log=False):

    '''
    Evaluates max_tau f(t+tau)*g(tau) or max_tau f(t-tau)g(tau) if inverse time is TRUE

    Parameters
    -----------

     t_val : double
        Time point

     f : Interpolation object
        First multiplier in convolution

     g : Interpolation object
        Second multiplier in convolution

     inverse_time : bool, None
        time direction. If True, then the f(t-tau)*g(tau) is calculated, otherwise,
        f(t+tau)*g(tau)

     return_log : bool
        If True, the logarithm will be returned


    Returns
    -------

     FG : Distribution
        The function to be integrated as Distribution object (interpolator)

    '''
    # return log is always True
    FG = _convolution_integrand(t_val, f, g, inverse_time, return_log=True)

    if FG == ttconf.BIG_NUMBER:
        res = [ttconf.BIG_NUMBER, 0]

    else:
        X = FG.x[FG.y.argmin()]
        Y = FG.y.min()
        res =  [Y, X]

    if not return_log:
        res[0] = np.exp(res[0])

    return res

def _evaluate_convolution(t_val, f, g,  n_integral = 100, inverse_time=None, return_log=False):
    """
    Calculate convolution F(t) = int { f(tau)g(t-tau) } dtau
    """

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
    """
    Node's position distribution function. This class extends the distribution
    class ind implements the convolution constructor.
    """

    @classmethod
    def convolve_fft(cls, node_interp, branch_interp, fft_grid_size=FFT_FWHM_GRID_SIZE, inverse_time=True):

        dt = max(branch_interp.one_mutation*0.005, min(node_interp.fwhm, branch_interp.fwhm)/fft_grid_size)
        ratio = node_interp.fwhm/branch_interp.fwhm
        if ratio < 1/fft_grid_size and 4*dt > node_interp.fwhm:
            ## node distribution is much narrower than the branch distribution, proceed as if node distribution is
            ## a delta distribution
            log_scale_node_interp = node_interp.integrate(return_log=True, a=node_interp.xmin,b=node_interp.xmax,n=max(100, len(node_interp.x))) #probability of node distribution 
            if inverse_time:
                x = branch_interp.x + node_interp._peak_pos
                dist = Distribution(x, branch_interp(x - node_interp._peak_pos) - log_scale_node_interp, min_width=max(node_interp.min_width, branch_interp.min_width), is_log=True)  
            else:
                x = - branch_interp.x + node_interp._peak_pos
                dist = Distribution(x, branch_interp(branch_interp.x) - log_scale_node_interp, min_width=max(node_interp.min_width, branch_interp.min_width), is_log=True)
            return dist
        elif ratio > fft_grid_size and 4*dt > branch_interp.fwhm:
            raise ValueError("ERROR: Unexpected behavior: branch distribution is much narrower than the node distribution.")
        else:
            b_effsupport = branch_interp.effective_support
            n_effsupport = node_interp.effective_support

            tmax = 2*max(b_effsupport[1]-b_effsupport[0], n_effsupport[1]-n_effsupport[0])

            Tb = np.arange(b_effsupport[0], b_effsupport[0] + tmax + dt, dt)
            if inverse_time:
                Tn = np.arange(n_effsupport[0], n_effsupport[0] + tmax + dt, dt)
                Tmin = node_interp.xmin
                Tmax = ttconf.MAX_BRANCH_LENGTH
            else:
                Tn = np.arange(n_effsupport[1] - tmax, n_effsupport[1] + dt, dt)
                Tmin = -ttconf.MAX_BRANCH_LENGTH
                Tmax = node_interp.xmax

            raw_len = len(Tb)
            fft_len = 2*raw_len

            fftb = branch_interp.fft(Tb, n=fft_len)
            fftn = node_interp.fft(Tn, n=fft_len, inverse_time=inverse_time)
            if inverse_time:
                fft_res = np.fft.irfft(fftb*fftn, fft_len)[:raw_len]
                Tres = Tn + Tb[0]
            else:
                fft_res = np.fft.irfft(fftb*fftn, fft_len)[::-1]
                fft_res = fft_res[raw_len:]
                Tres = Tn - Tb[0]

            # determine region in which we can trust the FFT convolution and avoid
            # inaccuracies due to machine precision. 1e-13 seems robust
            ind = fft_res>fft_res.max()*1e-13
            res = -np.log(fft_res[ind]) + branch_interp.peak_val + node_interp.peak_val - np.log(dt)
            Tres_cropped = Tres[ind]

            # extrapolate the tails exponentially: use margin last data points
            margin = np.minimum(3, Tres_cropped.shape[0]//3)
            if margin<1 or len(res)==0:
                import ipdb; ipdb.set_trace()
            else:
                left_slope = (res[margin]-res[0])/(Tres_cropped[margin]-Tres_cropped[0])
                right_slope = (res[-1]-res[-margin-1])/(Tres_cropped[-1]-Tres_cropped[-margin-1])

            # only extrapolate on the left when the slope is negative and we are not on the boundary
            if Tmin<Tres_cropped[0] and left_slope<0:
                Tleft = np.linspace(Tmin, Tres_cropped[0],10)[:-1]
                res_left = res[0] + left_slope*(Tleft - Tres_cropped[0])
            else:
                Tleft, res_left = [], []

            # only extrapolate on the right when the slope is positive and we are not on the boundary
            if Tres_cropped[-1]<Tmax and right_slope>0:
                Tright = np.linspace(Tres_cropped[-1], Tmax,10)[1:]
                res_right = res[-1] + right_slope*(Tright - Tres_cropped[-1])
            else: #otherwise
                Tright, res_right = [], []

            # instantiate the new interpolation object and return
            return cls(np.concatenate((Tleft,Tres_cropped,Tright)),
                    np.concatenate((res_left, res, res_right)),
                    is_log=True, kind='linear', assume_sorted=True)

    @classmethod
    def convolve(cls, node_interp, branch_interp, max_or_integral='integral',
                 n_grid_points = ttconf.NODE_GRID_SIZE, n_integral=ttconf.N_INTEGRAL,
                 inverse_time=True, rel_tol=0.05, yc=10):

        r'''
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
        n = n_grid_points//3
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


        if tmin>-1:
            grid_zero_left = tmin + (tmax-tmin)*np.linspace(0,0.01,11)**2
        else:
            grid_zero_left = [tmin]
        if tmax<1:
            grid_zero_right = tmax - (tmax-tmin)*np.linspace(0,0.01,11)**2
        else:
            grid_zero_right = [tmax]

        # make grid and calculate convolution
        t_grid_0 = np.unique(np.concatenate([grid_zero_left, grid_left[:-1], grid_center, grid_right[1:], grid_zero_right]))
        t_grid_0 = t_grid_0[(t_grid_0 > tmin-ttconf.TINY_NUMBER) & (t_grid_0 < tmax+ttconf.TINY_NUMBER)]

        # res0 - the values of the convolution (integral or max)
        # t_0  - the value, at which the res0 achieves maximum
        #        (when determining the maximum of the integrand, otherwise meaningless)
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

        # the interpolation object, which is used to store the value of the
        # grid, which maximizes the convolution (for 'max' option),
        # or flat -1 distribution (for 'integral' option)
        # this grid is the optimal branch length
        res_t = Distribution(t_grid_0, t_0, is_log=True,
                             min_width=node_interp.min_width, kind='linear')

        return res_y, res_t
