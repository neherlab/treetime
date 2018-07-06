import numpy as np
from Bio import Phylo

tavgii, davgii, tsqii, dtavgii, dsqii, sii = 0,1,2,3,4,5

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
    def __init__(self, T, tip_value = None,
                 branch_value = None, branch_variance = None):
        """
        Parameters
        ----------
         T : (Bio.Phylo.Tree)
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
        self.T = T
        # prep tree
        self.N = self.T.count_terminals()
        total_bl = 0
        for n in self.T.get_nonterminals(order='postorder'):
            for c in n:
                c.up=n
                total_bl+=c.branch_length
        self.T.root.up=None


        if tip_value is None:
            self.tip_value = lambda x:x.numdate if x.is_terminal() else None
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

        # gather tip indices
        for li, l in enumerate(self.T.get_terminals()):
            l._ii = np.array([li])
        for n in self.T.get_nonterminals(order='postorder'):
            n._ii = np.concatenate([c._ii for c in n])
            n._ii.sort()
            for c in n:
                c.up=n
        self.T.root.up=None

        # accumulate the covariance matrix by adding 'squares'
        M = np.zeros((self.N, self.N))
        for n in self.T.find_clades():
            if n == self.T.root:
                continue
            M[np.meshgrid(n._ii, n._ii)] += self.branch_variance(n)
        return M

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
        res = np.zeros(6, dtype=float)
        if n.is_terminal() and outgroup==False:
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

    def CovInv(self):
        """
        Inverse of the covariance matrix

        Returns
        -------

         H : (np.array)
            inverse of the covariance matrix.
        """
        self.recurse(full_matrix=True)
        return self.T.root.cinv

    def recurse(self, full_matrix=False):
        """
        recursion to calculate inverse covariance matrix
        """
        for n in self.T.get_nonterminals(order='postorder'):
            if full_matrix: M = np.zeros((len(n.ii), len(n.ii)))
            r = np.zeros((len(n._ii)))
            c_count = 0
            for c in n:
                ssq = self.branch_variance(c)
                if c.is_terminal():
                    if full_matrix:
                        M[c_count:c_count+nc, c_count:c_count+nc] = 1.0/ssq
                    r[c_count:c_count+nc] = 1.0/ssq
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
        for n in self.T.get_nonterminals(order='postorder'):
            Q = np.zeros(6, dtype=float)
            for c in n:
                tv = self.tip_value(c)
                bv = self.branch_value(c)
                var = self.branch_variance(c)
                Q+=self.propagate_averages(c, tv, bv, var)
            n.Q=Q

        for n in self.T.find_clades(order='preorder'):
            O = np.zeros(6, dtype=float)
            if n==self.T.root:
                continue

            for c in n.up:
                if c==n:
                    continue

                tv = self.tip_value(c)
                bv = self.branch_value(c)
                var = self.branch_variance(c)
                O += self.propagate_averages(c, tv, bv, var)

            if n.up!=self.T.root:
                c = n.up
                tv = self.tip_value(c)
                bv = self.branch_value(c)
                var = self.branch_variance(c)
                O += self.propagate_averages(c, tv, bv, var, outgroup=True)
            n.O = O

    def _regression(self, Q):
        """
        this function calculates the regression coefficients for a
        given vector containing the averages of tip and branch
        quantities.
        """
        slope = (Q[dtavgii] - Q[tavgii]*Q[davgii]/Q[sii]) \
                    /(Q[tsqii] - Q[tavgii]**2/Q[sii])
        intercept = (Q[davgii] - Q[tavgii]*slope)/Q[sii]
        chisq = 0.5*(Q[dsqii] - Q[davgii]**2/Q[sii]
                    - (Q[dtavgii] - Q[davgii]*Q[tavgii]/Q[sii])**2/(Q[tsqii]
                    - Q[tavgii]**2/Q[sii]))

        estimator_hessian = np.array([[Q[tsqii], Q[tavgii]], [Q[tavgii], Q[sii]]])

        return {'slope':slope, 'intercept':intercept,
                'chisq':chisq, 'hessian':estimator_hessian,
                'cov':np.linalg.inv(estimator_hessian)}


    def regression(self):
        self._calculate_averages()
        return self._regression(T.root.Q)


    def _optimal_root_along_branch(self, n, tv, bv, var):
        from scipy.optimize import minimize_scalar
        def chisq(x):
            tmpQ = self.propagate_averages(n, tv, bv*x, var*x) \
                 + self.propagate_averages(n, tv, bv*(1-x), var*(1-x), outgroup=True)
            return self._regression(tmpQ)['chisq']

        sol = minimize_scalar(chisq, bounds=(0,1), method="bounded")
        if sol["success"]:
            return sol['x'], sol['fun']
        else:
            return np.nan, np.inf

    def find_best_root(self):
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
        for n in self.T.find_clades():
            if n==self.T.root:
                continue

            tv = self.tip_value(n)
            bv = self.branch_value(n)
            var = self.branch_variance(n)
            x, chisq = self._optimal_root_along_branch(n, tv, bv, var)

            if (chisq<best_root["chisq"]):
                tmpQ = self.propagate_averages(n, tv, bv*x, var*x) \
                     + self.propagate_averages(n, tv, bv*(1-x), var*(1-x), outgroup=True)
                reg = self._regression(tmpQ)
                if reg["slope"]>0:
                    best_root = {"node":n, "split":x}
                    best_root.update(reg)

        # calculate differentials with respect to x
        deriv = []
        for dx in [-0.001, 0.001]:
            y = min(1.0, max(0.0, best_root["split"]+dx))
            tmpQ = self.propagate_averages(n, tv, bv*y, var*y) \
                 + self.propagate_averages(n, tv, bv*(1-y), var*(1-y), outgroup=True)
            reg = self._regression(tmpQ)
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


    def optimal_reroot(self):
        """
        determine the best root and reroot the tree to this value.
        Note that this can change the parent child relations of the tree
        and values associated with branches rather than nodes
        (e.g. confidence) might need to be re-evaluated afterwards
        """
        best_root = self.find_best_root()
        best_node = best_root["node"]
        x = best_root["split"]

        # create new node in the branch and root the tree to it
        new_node = Phylo.BaseTree.Clade()

        # insert the new node in the middle of the branch
        # by simple re-wiring the links on the both sides of the branch
        # and fix the branch lengths
        new_node.branch_length = best_node.branch_length*(1-x)
        new_node.up = best_node.up
        new_node.clades = [best_node]
        new_node.up.clades = [k if k != best_node else new_node
                              for k in best_node.up.clades]

        best_node.branch_length*(1-x)
        best_node.up = new_node
        new_node.rtt_regression = best_root
        self.T.root_with_outgroup(new_node)
        self.T.ladderize()
        return best_root



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    plt.ion()
    T = Phylo.read('../data/H3N2_NA_allyears_NA.20.nwk', 'newick')
    T.root_with_outgroup('A/Canterbury/58/2000|CY009150|09/05/2000|New_Zealand||H3N2/8-1416')

    dates = {}
    with open('../data/H3N2_NA_allyears_NA.20.metadata.csv') as ifile:
        ifile.readline()
        for line in ifile:
            if line[0]!='#':
                fields = line.strip().split(',')
                dates[fields[0]] = float(fields[1])


    for l in T.get_terminals():
        l.numdate = dates[l.name]
    mtc = TreeRegression(T, branch_variance = lambda x:(x.branch_length+0.005)/1700.0)
    reg = mtc.optimal_reroot()

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
