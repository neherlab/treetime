import numpy as np

def residual_filter(tt, n_iqd):
    terminals = tt.tree.get_terminals()
    clock_rate = tt.clock_model['slope']
    icpt = tt.clock_model['intercept']
    res = {}
    for node in terminals:
        if hasattr(node, 'raw_date_constraint') and  (node.raw_date_constraint is not None):
            res[node] = node.dist2root - clock_rate*np.mean(node.raw_date_constraint) - icpt

    residuals = np.array(list(res.values()))
    iqd = np.percentile(residuals,75) - np.percentile(residuals,25)
    bad_branch_count = 0
    for node,r in res.items():
        if abs(r)>n_iqd*iqd and node.up.up is not None:
            tt.logger('TreeTime.ClockFilter: marking %s as outlier, residual %f interquartile distances'%(node.name,r/iqd), 3, warn=True)
            node.bad_branch=True
            bad_branch_count += 1
        else:
            node.bad_branch=False

    return bad_branch_count
