"""
methods to calculate merger models for a time tree
"""
import numpy as np
import scipy.special as sf
from scipy.interpolate import interp1d
try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable
from . import config as ttconf
from .utils import clip


class Coalescent(object):
    """
    Class for adding the Kingman coalescent model to the tree time inference, this is optional.
    The coalescent model is based on the idea that certain tree structures are more likely given a specific population structure
    and this likelihood can be added to divergence time inference. The coalescent depends on the effective population size (:math:`Tc`)
    and the number of lineages at any point in time :math:`k(t)`.
    """

    def __init__(self, tree, Tc=0.001, logger=None, date2dist=None, n_branches_posterior=False):
        '''
        Initialize :math:`k(t)` and :math:`Tc` functions

        Parameters
        -----------
            Tc:    float
                time scale of coalescence / effective population size

            n_branches_posterior:  boolean
                True if the uncertainty of the divergence time estimates should be taken into consideration when calculating
                the number of lineages function :math:`k(t)`, False if current divergence times should be seen as fixed (default).
                Using the uncertainty should make :math:`k(t)` more smooth.

        '''

        super(Coalescent, self).__init__()
        self.tree = tree
        self.n_branches_posterior = n_branches_posterior
        self.calc_branch_count(posterior=n_branches_posterior)
        self.set_Tc(Tc)
        self.date2dist = date2dist
        if logger is None:
            def f(*args):
                print(*args)
            self.logger = f
        else:
            self.logger = logger

    def set_Tc(self, Tc, T=None):
        '''
        initialize the merger model with a coalescent time and calculate the integral of the merger rate

        Parameters
        ----------
            Tc:  float, iterable
                float if the coalescence rate is constant, if this should be approximated by a piece-wise linear
                merger rate (skyline method) an iterable with another argument T of the same shape is required
            T:  array
                an array of same shape as Tc that specifies the time pivots corresponding to Tc
                note that this array is ordered past to present corresponding to
                decreasing 'time before present' values
        '''
        if isinstance(Tc, Iterable):
            if len(Tc)==len(T):
                x = np.concatenate(([ttconf.BIG_NUMBER], T, [-ttconf.BIG_NUMBER]))
                y = np.concatenate(([Tc[0]], Tc, [Tc[-1]]))
                self.Tc = interp1d(x,y)
            else:
                self.logger("need Tc values and Timepoints of equal length",2,warn=True)
                self.Tc = interp1d([-ttconf.BIG_NUMBER, ttconf.BIG_NUMBER], [1e-5, 1e-5])
        else:
            self.Tc = interp1d([-ttconf.BIG_NUMBER, ttconf.BIG_NUMBER],
                               [Tc+ttconf.TINY_NUMBER, Tc+ttconf.TINY_NUMBER])
        self.calc_integral_merger_rate()



    def calc_branch_count(self, posterior=False):
        '''
        Calculates an interpolation object that maps time to the number of concurrent branches in the tree:
        :math:`k(t)`

        Parameters
        ----------
            posterior:  boolean
                If False use current best estimate of divergence times, else use posterior distributions of divergence times
                (If the marginal posterior time distribution of a node has been calculated this is used or
                approximated using the joint posterior time distribution)
        '''
        ## Divide merger events into either smooth merger events where a posterior likelihood distribution is known or
        ## delta events where either a date constraint for that node exists or the likelihood distribution is unknown.
        ## For delta distributions the corresponding nbranches step function can be calculated faster as the nodes can be
        ## sorted by time and mergers added or subtracted from the previous time, for smooth distributions when a new merger
        ## event occurs the previous distribution must be evaluated at the corresponding position.

        self.tree_events = sorted([(n.time_before_present, len(n.clades)-1)
                        for n in self.tree.find_clades() if not n.bad_branch],
                        key=lambda x:-x[0])

        tree_delta_events = []
        tree_smooth_events = []

        if not posterior:
            tree_delta_events = self.tree_events
        else:
            y_power = np.array([-8, -4, -3, -2, 0, 2, 3, 4, 8])
            y_points= np.exp(y_power)/(1 + np.exp(y_power))
            for n in self.tree.find_clades():
                cdf_function=None
                # use cdf function if exists and not from a delta function
                if hasattr(n, 'marginal_inverse_cdf') and not n.marginal_pos_LH.is_delta:
                    cdf_function=n.marginal_inverse_cdf
                elif hasattr(n, 'joint_inverse_cdf') and (n.date_constraint is None or not n.date_constraint.is_delta):
                    cdf_function=n.joint_inverse_cdf

                if cdf_function is not None:
                    x_vals = np.concatenate([[-ttconf.BIG_NUMBER], cdf_function(y_points), [ttconf.BIG_NUMBER]])
                    y_vals = np.concatenate([ [(len(n.clades)-1),(len(n.clades)-1)], (1-y_points[1:-1]), [0,0]])
                    tree_smooth_events +=  [interp1d(x_vals, y_vals, kind="linear")]
                else:
                    tree_delta = [(n.time_before_present, len(n.clades)-1)]
                    tree_delta_events += tree_delta
            tree_delta_events= sorted(tree_delta_events, key=lambda x:-x[0])

        if tree_delta_events:
            # collapse multiple events at one time point into sum of changes
            from collections import defaultdict
            dn_branch = defaultdict(int)
            for (t, dn) in tree_delta_events:
                dn_branch[t]+=dn
            unique_mergers = np.array(sorted(dn_branch.items(), key = lambda x:-x[0]))

            # calculate the branch count at each point summing the delta branch counts
            nbranches_discrete = [[ttconf.BIG_NUMBER, 1], [unique_mergers[0,0]+ttconf.TINY_NUMBER, 1]]
            for ti, (t, dn) in enumerate(unique_mergers[:-1]):
                new_n = nbranches_discrete[-1][1]+dn
                next_t = unique_mergers[ti+1,0]+ttconf.TINY_NUMBER
                nbranches_discrete.append([t, new_n])
                nbranches_discrete.append([next_t, new_n])

            new_n += unique_mergers[-1,1]
            nbranches_discrete.append([unique_mergers[ti+1,0], new_n])
            nbranches_discrete.append([-ttconf.BIG_NUMBER, new_n])
            nbranches_discrete=np.array(nbranches_discrete)
            nbranches_discrete = interp1d(nbranches_discrete[:,0], nbranches_discrete[:,1], kind='linear')

        if tree_smooth_events:
            # add all smooth events by evaluating at all unique x points
            x_tot = np.unique(np.concatenate([t.x for t in tree_smooth_events]))
            y_tot = np.array([t(x_tot) for t in tree_smooth_events]).sum(axis=0)
            nbranches_smooth = interp1d(x_tot, y_tot +1, kind='linear')
            if tree_delta_events:
                # join smooth and delta merger events into one distribution object
                x_tot = np.unique(np.concatenate([nbranches_discrete.x, nbranches_smooth.x]))
                y_tot = nbranches_discrete(x_tot) + nbranches_smooth(x_tot)
                # if both delta and smooth event objects exist must remove the initial starting value so not double
                self.nbranches = interp1d(x_tot, y_tot -1, kind='linear')
            else:
                self.nbranches = nbranches_smooth
        else:
            self.nbranches = nbranches_discrete

        self.tree_events = np.array(self.tree_events)

    def calc_integral_merger_rate(self):
        '''
        calculates the integral :math:`int_0^t (k(t')-1)/2Tc(t') dt` and stores it as
        self.integral_merger_rate. This differences of this quantity evaluated at
        different times points are the cost of a branch.
        '''
        # integrate the piecewise constant branch count function.
        tvals = np.unique(self.nbranches.x[1:-1])
        rate = self.branch_merger_rate(tvals)
        avg_rate = 0.5*(rate[1:] + rate[:-1])
        cost = np.concatenate(([0],np.cumsum(np.diff(tvals)*avg_rate)))
        # make interpolation objects for the branch count and its integral
        # the latter is scaled by 0.5/Tc
        # need to add extra point at very large time before present to
        # prevent 'out of interpolation range' errors
        self.integral_merger_rate = interp1d(np.concatenate(([-ttconf.BIG_NUMBER], tvals,[ttconf.BIG_NUMBER])),
                                  np.concatenate(([cost[0]], cost,[cost[-1]])), kind='linear')

    def branch_merger_rate(self, t):
        '''
        rate at which one particular branch merges with any other branch at time t,
        in the Kingman model this is: :math:`\kappa(t) = (k(t)-1)/(2Tc(t))`
        '''
        # note that we always have a positive merger rate by capping the
        # number of branches at 0.5 from below. in these regions, the
        # function should only be called if the tree changes.
        return 0.5*np.maximum(0.5,self.nbranches(t)-1.0)/self.Tc(t)

    def total_merger_rate(self, t):
        '''
        rate at which any branch merges with any other branch at time t,
        in the Kingman model this is: :math:`\lambda(t) = k(t)(k(t)-1)/(2Tc(t))`
        '''
        # note that we always have a positive merger rate by capping the
        # number of branches at 0.5 from below. in these regions, the
        # function should only be called if the tree changes.
        nlineages = np.maximum(0.5,self.nbranches(t)-1.0)
        return 0.5*nlineages*(nlineages+1)/self.Tc(t)


    def cost(self, t_node, branch_length, multiplicity=2.0):
        '''
        returns the cost associated with a branch starting with divergence time t_node (:math:`t_n`)
        having a branch length :math:`\\tau`.
        This is equivalent to the probability of there being no merger on that branch and a merger at the end of the branch,
        calculated in the negative log
        :math:`-log(\lambda(t_n+ \\tau)^{(m-1)/m}) + \int_{t_n}^{t_n+ \\tau} \kappa(t) dt`, where m is the multiplicity

        Parameters
        ----------
            t_node: float
                time of the node
            branch_length: float
                branch length, determines when this branch merges with sister
            multiplicity: int
                2 if merger is binary, higher if this is a polytomy
        '''
        merger_time = t_node + np.maximum(0,branch_length)
        return self.integral_merger_rate(merger_time) - self.integral_merger_rate(t_node)\
                 - np.log(self.total_merger_rate(merger_time))*(multiplicity-1.0)/multiplicity


    def node_contribution(self, node, t, multiplicity=None):
        '''
        returns the contribution of node at time t to cost of merging branch that node is parent of
        '''
        from treetime.node_interpolator import NodeInterpolator
        if multiplicity is None:
            multiplicity = len(node.clades)
        # the number of mergers is 'number of children' - 1
        multiplicity -= 1.0
        y = (self.integral_merger_rate(t) - np.log(self.total_merger_rate(t)))*multiplicity
        return NodeInterpolator(t, y, is_log=True)


    def total_LH(self):
        LH = 0.0 #np.log(self.total_merger_rate([node.time_before_present for node in self.tree.get_nonterminals()])).sum()
        for node in self.tree.find_clades():
            if node.up:
                LH -= self.cost(node.time_before_present, node.branch_length)
        return LH


    def optimize_Tc(self):
        '''
        determines the coalescent time scale Tc that optimizes the coalescent likelihood of the tree
        (product of the cost of coalescence of all nodes)
        '''
        from scipy.optimize import minimize_scalar
        initial_Tc = self.Tc
        def cost(logTc):
            self.set_Tc(np.exp(logTc))
            return -self.total_LH()

        sol = minimize_scalar(cost, bracket=[-20.0, 2.0], method='brent')
        if "success" in sol and sol["success"]:
            self.set_Tc(np.exp(sol['x']))
        else:
            self.logger("merger_models:optimize_Tc: optimization of coalescent time scale failed: " + str(sol), 0, warn=True)
            self.set_Tc(initial_Tc.y, T=initial_Tc.x)


    def optimize_skyline(self, n_points=20, stiffness=2.0, method = 'SLSQP',
                         tol=0.03, regularization=10.0, **kwarks):
        '''
        optimize the trajectory of the clock rate 1./T_c to maximize the
        coalescent likelihood, this is the product of the cost of coalescence of all nodes

        Parameters
        ----------
            n_points: int
                number of pivots of the Tc interpolation object
            stiffness: float
                penalty for rapid changes in log(Tc)
            methods: str
                method used to optimize, see documentation of scipy.optimize.minimize for options
            tol: float
                optimization tolerance
            regularization: float
                cost of moving log(Tc) outsize of the range [-100,0]
                merger rate is measured in branch length units, no
                plausible rates should ever be outside this window
        '''
        self.logger("Coalescent:optimize_skyline:... current LH: %f"%self.total_LH(),2)
        from scipy.optimize import minimize
        initial_Tc = self.Tc
        tvals = np.linspace(self.tree_events[0,0], self.tree_events[-1,0], n_points)
        def cost(logTc):
            # cap log Tc to avoid under or overflow and nan in logs
            self.set_Tc(np.exp(clip(logTc, -200, 100)), tvals)
            neglogLH = -self.total_LH() + stiffness*np.sum(np.diff(logTc)**2) \
                       + np.sum((logTc>0)*logTc)*regularization\
                       - np.sum((logTc<-100)*logTc)*regularization
            return neglogLH

        sol = minimize(cost, np.ones_like(tvals)*np.log(self.Tc.y.mean()), method=method, tol=tol)
        if "success" in sol and sol["success"]:
            dlogTc = 0.1
            opt_logTc = sol['x']
            dcost = []
            for ii in range(len(opt_logTc)):
                tmp = opt_logTc.copy()
                tmp[ii]+=dlogTc
                cost_plus = cost(tmp)
                tmp[ii]-=2*dlogTc
                cost_minus = cost(tmp)
                dcost.append([cost_minus, cost_plus])

            dcost = np.array(dcost)
            optimal_cost = cost(opt_logTc)
            self.confidence = dlogTc/np.sqrt(np.abs(2*optimal_cost - dcost[:,0] - dcost[:,1]))
            self.logger("Coalescent:optimize_skyline:...done. new LH: %f"%self.total_LH(),2)
        else:
            self.set_Tc(initial_Tc.y, T=initial_Tc.x)
            self.confidence = [np.nan for i in initial_Tc.x]
            self.logger("Coalescent:optimize_skyline:...failed:"+str(sol),0, warn=True)


    def skyline_empirical(self, gen=1.0, n_points = 20):
        '''
        returns the skyline, i.e., an estimate of the inverse rate of coalesence.
        Here, the skyline is estimated from a sliding window average of the observed
        mergers, i.e., without reference to the coalescence likelihood.

        Parameters
        ----------
            gen: float
                number of generations per year
            n_points: int
        '''
        merger_times = np.array(self.tree_events[self.tree_events[:,1]>0, 0])
        nlineages = self.nbranches(merger_times -ttconf.TINY_NUMBER)
        expected_merger_density = nlineages*(nlineages-1)*0.5

        nmergers = len(merger_times)
        et = merger_times
        ev = 1.0/expected_merger_density
        # reduce the window size if there are few events in the tree
        if 2*n_points>len(expected_merger_density):
            n_points = len(ev)//4

        # smoothes with a sliding window over data points
        avg = np.sum(ev)/np.abs(et[0]-et[-1])
        dt = et[0]-et[-1]
        mid_points = np.concatenate(([et[0]-0.5*(et[1]-et[0])],
                                      0.5*(et[1:] + et[:-1]),
                                     [et[-1]+0.5*(et[-1]-et[-2])]))

        # this smoothes the ratio of expected and observed merger rate
        # epsilon is added to avoid division by 0 and to normalize Tc
        epsilon= (1/n_points)*dt/nmergers
        self.Tc_inv = interp1d(mid_points[n_points:-n_points],
                        [np.sum(ev[(et>=l)&(et<u)])/(u-l+epsilon)
                        for u,l in zip(mid_points[:-2*n_points],mid_points[2*n_points:])])

        return interp1d(self.date2dist.to_numdate(self.Tc_inv.x), gen/self.date2dist.clock_rate/self.Tc_inv.y)


    def skyline_inferred(self, gen=1.0, confidence=False):
        '''
        return the skyline, i.e., an estimate of the inverse rate of coalesence.
        This function merely returns the merger rate self.Tc that was set or
        estimated by other means. If it was determined using self.optimize_skyline,
        the returned skyline will maximize the coalescent likelihood.

        Parameters
        ----------
            gen: float
                number of generations per year. Unit of time is branch length,
                hence this needs to be the inverse substitution rate per generation
            confidence: boolean, float
                False, or number of standard deviations of confidence intervals
        '''
        if len(self.Tc.x)<=2:
            print("no skyline has been inferred, returning constant population size")
            return gen/self.date2dist.clock_rate*self.Tc.y[-1]

        skyline = interp1d(self.date2dist.to_numdate(self.Tc.x[1:-1]), gen/self.date2dist.clock_rate*self.Tc.y[1:-1])
        if confidence and hasattr(self, 'confidence'):
            conf = [skyline.y*np.exp(-confidence*self.confidence), skyline.y*np.exp(confidence*self.confidence)]
            return skyline, conf
        else:
            return skyline, None





