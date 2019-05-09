import numpy as np
from Bio import Phylo

tavgii, davgii, tsqii, dtavgii, dsqii, sii = 0,1,2,3,4,5

def base_regression(Q, slope=None):
    """
    this function calculates the regression coefficients for a
    given vector containing the averages of tip and branch
    quantities.

    Parameters
    ----------
    Q : numpy.array
        vector with
    slope : None, optional
        Description

    Returns
    -------
    TYPE
        Description
    """
    if slope is None:
        if (Q[tsqii] - Q[tavgii]**2/Q[sii])>0:
            slope = (Q[dtavgii] - Q[tavgii]*Q[davgii]/Q[sii]) \
                /(Q[tsqii] - Q[tavgii]**2/Q[sii])
        else:
            raise ValueError("No variation in sampling dates! Please specify your clock rate explicitly.")
        only_intercept=False
    else:
        only_intercept=True

    intercept = (Q[davgii] - Q[tavgii]*slope)/Q[sii]
    if (Q[tsqii] - Q[tavgii]**2/Q[sii])>0:
        chisq = 0.5*(Q[dsqii] - Q[davgii]**2/Q[sii] - (Q[dtavgii] - Q[davgii]*Q[tavgii]/Q[sii])**2/(Q[tsqii] - Q[tavgii]**2/Q[sii]))
    else:
        chisq = 0.5*(Q[dsqii] - Q[davgii]**2/Q[sii])

    if only_intercept:
        return {'slope':slope, 'intercept':intercept,
                'chisq': chisq}

    estimator_hessian = np.array([[Q[tsqii], Q[tavgii]], [Q[tavgii], Q[sii]]])

    return {'slope':slope, 'intercept':intercept,
            'chisq':chisq, 'hessian':estimator_hessian,
            'cov':np.linalg.inv(estimator_hessian)}


class TreeRegression(object):
    """TreeRegression
    This class implements an efficient regression method
    for quantity associated with tips and one that changes
    in an additive manner along the branches of the tree,
    e.g. the distance to the root. This implemented
    algorithm take into account the correlation structure
    of the data under the assumptions that variance increase
    linearly along branches as well.
    """
    def __init__(self, tree_in, tip_value = None,
                 branch_value = None, branch_variance = None):
        """
        Parameters
        ----------
         T : (Bio.Phylo.tree)
            Tree for which the covariances and regression
            are to be calculated.

         tip_value : (callable)
            function that for each tip returns the value to
            be used in the regression.

         branch_value : (callable)
            function that for each node of the tree returns
            the contribution of this branch to the value of
            the subtending tips.

         variance_function : (callable)
            function that for each node of the tree returns
            the accumulated variance
        """
        super(TreeRegression, self).__init__()
        self.tree = tree_in
        # prep tree
        for li, l in enumerate(self.tree.get_terminals()):
            l._ii = np.array([li])
        total_bl = 0
        for n in self.tree.get_nonterminals(order='postorder'):
            n._ii = np.concatenate([c._ii for c in n])
            n._ii.sort()
            for c in n:
                c.up=n
                total_bl+=c.branch_length
        self.tree.root.up=None
        self.N = self.tree.root._ii.shape[0]

        if tip_value is None:
            self.tip_value = lambda x:np.mean(x.numdate) if x.is_terminal() else None
        else:
            self.tip_value = tip_value

        if branch_value is None:
            self.branch_value = lambda x:x.branch_length
        else:
            self.branch_value = branch_value

        if branch_variance is None:
            # provide a default equal to the branch_length (Poisson) and add
            # a tenth of the average branch length to avoid numerical instabilities and division by 0.
            self.branch_variance = lambda x:x.branch_length + 0.05*total_bl/self.N
        else:
            self.branch_variance = branch_variance


    def Cov(self):
        """
        calculate the covariance matrix of the tips assuming variance
        has accumulated along branches of the tree accoriding to the
        the provided
        Returns
        -------

         M : (np.array)
            covariance matrix with tips arranged standard transersal order.
        """
        # accumulate the covariance matrix by adding 'squares'
        M = np.zeros((self.N, self.N))
        for n in self.tree.find_clades():
            if n == self.tree.root:
                continue
            M[np.meshgrid(n._ii, n._ii)] += self.branch_variance(n)
        return M


    def CovInv(self):
        """
        Inverse of the covariance matrix

        Returns
        -------

         H : (np.array)
            inverse of the covariance matrix.
        """
        self.recurse(full_matrix=True)
        return self.tree.root.cinv


    def recurse(self, full_matrix=False):
        """
        recursion to calculate inverse covariance matrix

        Parameters
        ----------
        full_matrix : bool, optional
            if True, the entire inverse matrix is calculated. otherwise, only the weighing vector.
        """
        for n in self.tree.get_nonterminals(order='postorder'):
            n_leaves = len(n._ii)
            if full_matrix: M = np.zeros((n_leaves, n_leaves), dtype=float)
            r = np.zeros(n_leaves, dtype=float)
            c_count = 0
            for c in n:
                ssq = self.branch_variance(c)
                nc = len(c._ii)
                if c.is_terminal():
                    if full_matrix:
                        M[c_count, c_count] = 1.0/ssq
                    r[c_count] = 1.0/ssq
                else:
                    if full_matrix:
                        M[c_count:c_count+nc, c_count:c_count+nc] = c.cinv - ssq*np.outer(c.r,c.r)/(1+ssq*c.s)
                    r[c_count:c_count+nc] = c.r/(1+ssq*c.s)
                c_count += nc

            if full_matrix: n.cinv = M
            n.r = r #M.sum(axis=1)
            n.s = n.r.sum()


    def _calculate_averages(self):
        """
        calculate the weighted sums of the tip and branch values and
        their second moments.
        """
        for n in self.tree.get_nonterminals(order='postorder'):
            Q = np.zeros(6, dtype=float)
            for c in n:
                tv = self.tip_value(c)
                bv = self.branch_value(c)
                var = self.branch_variance(c)
                Q += self.propagate_averages(c, tv, bv, var)
            n.Q=Q

        for n in self.tree.find_clades(order='preorder'):
            O = np.zeros(6, dtype=float)
            if n==self.tree.root:
                n.Qtot = n.Q
                continue

            for c in n.up:
                if c==n:
                    continue

                tv = self.tip_value(c)
                bv = self.branch_value(c)
                var = self.branch_variance(c)
                O += self.propagate_averages(c, tv, bv, var)

            if n.up!=self.tree.root:
                c = n.up
                tv = self.tip_value(c)
                bv = self.branch_value(c)
                var = self.branch_variance(c)
                O += self.propagate_averages(c, tv, bv, var, outgroup=True)
            n.O = O

            if not n.is_terminal():
                tv = self.tip_value(n)
                bv = self.branch_value(n)
                var = self.branch_variance(n)
                n.Qtot = n.Q + self.propagate_averages(n, tv, bv, var, outgroup=True)


    def propagate_averages(self, n, tv, bv, var, outgroup=False):
        """
        This function implements the propagation of the means,
        variance, and covariances along a branch. It operates
        both towards the root and tips.

        Parameters
        ----------
         n : (node)
            the branch connecting this node to its parent is used
            for propagation
         tv : (float)
            tip value. Only required if not is terminal
         bl : (float)
            branch value. The increment of the tree associated quantity'

         var : (float)
            the variance increment along the branch

        Returns
        -------
         Q : (np.array)
            a vector of length 6 containing the updated quantities
        """
        if n.is_terminal() and outgroup==False:
            if tv is None or np.isinf(tv) or np.isnan(tv):
                res = np.array([0, 0, 0, 0, 0, 0])
            elif var==0:
                res = np.array([np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])
            else:
                res = np.array([
                    tv/var,
                    bv/var,
                    tv**2/var,
                    bv*tv/var,
                    bv**2/var,
                    1.0/var], dtype=float)
        else:
            tmpQ = n.O if outgroup else n.Q
            denom = 1.0/(1+var*tmpQ[sii])
            res = np.array([
                tmpQ[tavgii]*denom,
                (tmpQ[davgii] + bv*tmpQ[sii])*denom,
                tmpQ[tsqii] - var*tmpQ[tavgii]**2*denom,
                tmpQ[dtavgii] + tmpQ[tavgii]*bv - var*tmpQ[tavgii]*(tmpQ[davgii] + bv*tmpQ[sii])*denom,
                tmpQ[dsqii] + 2*bv*tmpQ[davgii] + bv**2*tmpQ[sii] - var*(tmpQ[davgii]**2 + 2*bv*tmpQ[davgii]*tmpQ[sii] + bv**2*tmpQ[sii]**2)*denom,
                tmpQ[sii]*denom]
            )

        return res

    def explained_variance(self):
        """calculate standard explained variance

        Returns
        -------
        float
            r-value of the root-to-tip distance and time.
            independent of regression model, but dependent on root choice
        """
        self.tree.root._v=0
        for n in self.tree.get_nonterminals(order='preorder'):
            for c in n:
                c._v = n._v + self.branch_value(c)
        raw = np.array([(self.tip_value(n), n._v) for n in self.tree.get_terminals()
                         if self.tip_value(n) is not None])
        return np.corrcoef(raw.T)[0,1]


    def regression(self, slope=None):
        """regress tip values against branch values

        Parameters
        ----------
        slope : None, optional
            if given, the slope isn't optimized

        Returns
        -------
        dict
            regression parameters
        """
        self._calculate_averages()

        clock_model = base_regression(self.tree.root.Q, slope=slope)
        clock_model['r_val'] = self.explained_variance()

        return clock_model



    def find_best_root(self, force_positive=True, slope=None):
        """
        determine the position on the tree that minimizes the bilinear
        product of the inverse covariance and the data vectors.

        Returns
        -------
         best_root : (dict)
            dictionary with the node, the fraction `x` at which the branch
            is to be split, and the regression parameters
        """
        self._calculate_averages()
        best_root = {"chisq": np.inf}
        for n in self.tree.find_clades():
            if n==self.tree.root:
                continue

            tv = self.tip_value(n)
            bv = self.branch_value(n)
            var = self.branch_variance(n)
            x, chisq = self._optimal_root_along_branch(n, tv, bv, var, slope=slope)
            if (chisq<best_root["chisq"]):
                tmpQ = self.propagate_averages(n, tv, bv*x, var*x) \
                     + self.propagate_averages(n, tv, bv*(1-x), var*(1-x), outgroup=True)
                reg = base_regression(tmpQ, slope=slope)
                if reg["slope"]>=0 or (force_positive==False):
                    best_root = {"node":n, "split":x}
                    best_root.update(reg)

        if 'node' not in best_root:
            print("TreeRegression.find_best_root: No valid root found!", force_positive)
            return None

        if 'hessian' in best_root:
            # calculate differentials with respect to x
            deriv = []
            n = best_root["node"]
            tv = self.tip_value(n)
            bv = self.branch_value(n)
            var = self.branch_variance(n)
            for dx in [-0.001, 0.001]:
                y = min(1.0, max(0.0, best_root["split"]+dx))
                tmpQ = self.propagate_averages(n, tv, bv*y, var*y) \
                     + self.propagate_averages(n, tv, bv*(1-y), var*(1-y), outgroup=True)
                reg = base_regression(tmpQ, slope=slope)
                deriv.append([y,reg['chisq'], tmpQ[tavgii], tmpQ[davgii]])

            estimator_hessian = np.zeros((3,3))
            estimator_hessian[:2,:2] = best_root['hessian']
            estimator_hessian[2,2] = (deriv[0][1] + deriv[1][1] - 2.0*best_root['chisq'])/(deriv[0][0] - deriv[1][0])**2
            # estimator_hessian[2,0] = (deriv[0][2] - deriv[1][2])/(deriv[0][0] - deriv[1][0])
            # estimator_hessian[2,1] = (deriv[0][3] - deriv[1][3])/(deriv[0][0] - deriv[1][0])
            estimator_hessian[0,2] = estimator_hessian[2,0]
            estimator_hessian[1,2] = estimator_hessian[2,1]
            best_root['hessian'] = estimator_hessian
            best_root['cov'] = np.linalg.inv(estimator_hessian)

        return best_root


    def _optimal_root_along_branch(self, n, tv, bv, var, slope=None):
        from scipy.optimize import minimize_scalar
        def chisq(x):
            tmpQ = self.propagate_averages(n, tv, bv*x, var*x) \
                 + self.propagate_averages(n, tv, bv*(1-x), var*(1-x), outgroup=True)
            return base_regression(tmpQ, slope=slope)['chisq']

        chisq_prox = np.inf if n.is_terminal() else base_regression(n.Qtot, slope=slope)['chisq']
        chisq_dist = np.inf if n==self.tree.root else base_regression(n.up.Qtot, slope=slope)['chisq']

        grid = np.linspace(0.001,0.999,6)
        chisq_grid = np.array([chisq(x) for x in grid])
        min_chisq = chisq_grid.min()
        if chisq_prox<=min_chisq:
            return 0.0, chisq_prox
        elif chisq_dist<=min_chisq:
            return 1.0, chisq_dist
        else:
            ii = np.argmin(chisq_grid)
            bounds = (0 if ii==0 else grid[ii-1], 1.0 if ii==len(grid)-1 else grid[ii+1])
            sol = minimize_scalar(chisq, bounds=bounds, method="bounded")
            if sol["success"]:
                return sol['x'], sol['fun']
            else:
                return np.nan, np.inf


    def optimal_reroot(self, force_positive=True, slope=None):
        """
        determine the best root and reroot the tree to this value.
        Note that this can change the parent child relations of the tree
        and values associated with branches rather than nodes
        (e.g. confidence) might need to be re-evaluated afterwards

        Parameters
        ----------
        force_positive : bool, optional
            if True, the search for a root will only consider positive rate estimates

        slope : float, optional
            if given, it will find the optimal root given a fixed rate. If slope==0, this
            corresponds to minimal root-to-tip variance rooting (min_dev)

        Returns
        -------
        dict
            regression parameters
        """
        best_root = self.find_best_root(force_positive=force_positive, slope=slope)
        best_node = best_root["node"]

        x = best_root["split"]
        if x<1e-5:
            new_node = best_node
        elif x>1.0-1e-5:
            new_node = best_node.up
        else:
            # create new node in the branch and root the tree to it
            new_node = Phylo.BaseTree.Clade()

            # insert the new node in the middle of the branch
            # by simple re-wiring the links on the both sides of the branch
            # and fix the branch lengths
            new_node.branch_length = best_node.branch_length*(1-x)
            new_node.up = best_node.up
            new_node.clades = [best_node]
            new_node.up.clades = [k if k!=best_node else new_node
                                  for k in best_node.up.clades]

            best_node.branch_length *= x
            best_node.up = new_node

        new_node.rtt_regression = best_root
        self.tree.root_with_outgroup(new_node)

        self.tree.ladderize()
        for n in self.tree.get_nonterminals(order='postorder'):
            for c in n:
                c.up=n
        return best_root


    def clock_plot(self, add_internal=False, ax=None, regression=None,
                   confidence=True, n_sigma = 2, fs=14):
        """Plot root-to-tip distance vs time as a basic time-tree diagnostic

        Parameters
        ----------
        add_internal : bool, optional
            add internal nodes. this will only work if the tree has been dated already
        ax : None, optional
            an matplotlib axis to plot into. if non provided, a new figure is opened
        regression : None, optional
            a dict containing parameters of a root-to-tip vs time regression as
            returned by the function base_regression
        confidence : bool, optional
            add confidence area to the regression line
        n_sigma : int, optional
            number of standard deviations for the confidence area.
        fs : int, optional
            fontsize

        """
        import matplotlib.pyplot as plt
        if ax is None:
            plt.figure()
            ax=plt.subplot(111)

        self.tree.root._v=0
        for n in self.tree.get_nonterminals(order='preorder'):
            for c in n:
                c._v = n._v + self.branch_value(c)

        tips = self.tree.get_terminals()
        internal = self.tree.get_nonterminals()

        # get values of terminals
        xi = np.array([self.tip_value(n) for n in tips])
        yi = np.array([n._v for n in tips])
        ind = np.array([n.bad_branch  if hasattr(n, 'bad_branch') else False for n in tips])
        if add_internal:
            xi_int = np.array([n.numdate for n in internal])
            yi_int = np.array([n._v for n in internal])
            ind_int = np.array([n.bad_branch  if hasattr(n, 'bad_branch') else False  for n in internal])

        if regression:
            # plot regression line
            t_mrca = -regression['intercept']/regression['slope']
            if add_internal:
                time_span = np.max(xi_int[~ind_int]) - np.min(xi_int[~ind_int])
                x_vals = np.array([max(np.min(xi_int[~ind_int]), t_mrca) - 0.1*time_span, np.max(xi[~ind])+0.05*time_span])
            else:
                time_span = np.max(xi[~ind]) - np.min(xi[~ind])
                x_vals = np.array([max(np.min(xi[~ind]), t_mrca) - 0.1*time_span, np.max(xi[~ind]+0.05*time_span)])

            # plot confidence interval
            if confidence and 'cov' in regression:
                x_vals = np.linspace(x_vals[0], x_vals[1], 100)
                y_vals = regression['slope']*x_vals + regression['intercept']
                dev = n_sigma*np.array([np.sqrt(regression['cov'][:2,:2].dot(np.array([x, 1])).dot(np.array([x,1]))) for x in x_vals])
                dev_slope = n_sigma*np.sqrt(regression['cov'][0,0])
                ax.fill_between(x_vals, y_vals-dev, y_vals+dev, alpha=0.2)
                dp = np.array([regression['intercept']/regression['slope']**2,-1./regression['slope']])
                dev_rtt = n_sigma*np.sqrt(regression['cov'][:2,:2].dot(dp).dot(dp))

            else:
                dev_rtt = None
                dev_slope = None

            ax.plot(x_vals, regression['slope']*x_vals + regression['intercept'],
                    label = r"$y=\alpha + \beta t$"+"\n"+
                            r"$\beta=$%1.2e"%(regression["slope"])
                            + ("+/- %1.e"%dev_slope if dev_slope else "") +
                            "\nroot date: %1.1f"%(-regression['intercept']/regression['slope']) +
                            ("+/- %1.2f"%dev_rtt if dev_rtt else ""))


        ax.scatter(xi[~ind], yi[~ind], label=("tips" if add_internal else None))
        if ind.sum():
            try:
                # note: this is treetime specific
                tmp_x = np.array([np.mean(n.raw_date_constraint) if n.raw_date_constraint else None
                                  for n in self.tree.get_terminals()])
                ax.scatter(tmp_x[ind], yi[ind], label="ignored tips", c='r')
            except:
                pass
        if add_internal:
            ax.scatter(xi_int[~ind_int], yi_int[~ind_int], label="internal nodes")

        ax.set_ylabel('root-to-tip distance', fontsize=fs)
        ax.set_xlabel('date', fontsize=fs)
        ax.ticklabel_format(useOffset=False)
        ax.tick_params(labelsize=fs*0.8)
        ax.set_ylim([0, 1.1*np.max(yi)])
        plt.tight_layout()
        plt.legend(fontsize=fs*0.8)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import time
    plt.ion()
    # tree_file = '../data/H3N2_NA_allyears_NA.20.nwk'
    # date_file = '../data/H3N2_NA_allyears_NA.20.metadata.csv'
    tree_file = '../data/ebola.nwk'
    date_file = '../data/ebola.metadata.csv'

    T = Phylo.read(tree_file, 'newick')
    #T.root_with_outgroup('A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416')

    dates = {}
    with open(date_file) as ifile:
        ifile.readline()
        for line in ifile:
            if line[0]!='#':
                fields = line.strip().split(',')
                dates[fields[0]] = float(fields[1])


    for l in T.get_terminals():
        l.numdate = dates[l.name]
    branch_variance = lambda x:(x.branch_length+(0.0005 if x.is_terminal() else 0.0))/19000.0
    #branch_variance = lambda x:(x.branch_length+(0.005 if x.is_terminal() else 0.0))/1700.0
    #branch_variance = lambda x:1.0 if x.is_terminal() else 0.0
    tstart = time.time()
    mtc = TreeRegression(T, branch_variance = branch_variance)
    print(time.time()-tstart)
    reg = mtc.optimal_reroot()
    print(time.time()-tstart)
    print(reg)

    plt.figure()
    ti = []
    rtt = []
    T.root.rtt=0
    for n in T.get_nonterminals(order='preorder'):
        for c in n:
            c.rtt = n.rtt + c.branch_length
    for l in T.get_terminals():
        ti.append(l.numdate)
        rtt.append(l.rtt)

    ti = np.array(ti)
    rtt = np.array(rtt)
    plt.plot(ti, rtt)
    plt.plot(ti, reg["slope"]*ti + reg["intercept"])
    plt.show()

    Phylo.draw(T)
