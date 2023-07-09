import numpy as np
import pandas as pd

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
    outliers = {}
    for node,r in res.items():
        if abs(r)>n_iqd*iqd and node.up.up is not None:
            node.bad_branch=True
            outliers[node.name] = {'tau':(node.dist2root - icpt)/clock_rate,  'avg_date': np.mean(node.raw_date_constraint),
                                'exact_date': node.raw_date_constraint if type(node) is float else None,
                                'residual': r/iqd}
        else:
            node.bad_branch=False

    tt.outliers=None
    if len(outliers):
        outlier_df = pd.DataFrame(outliers).T.loc[:,['avg_date', 'tau', 'residual']]\
                                .rename(columns={'avg_date':'given_date', 'tau':'apparent_date'})
        tt.logger("Clock_filter.residual_filter marked the following outliers:", 2, warn=True)
        if tt.verbose>=2:
            print(outlier_df)
        tt.outliers = outlier_df
    return len(outliers)

def local_filter(tt, z_score_threshold):
    tt.logger(f"TreeTime.ClockFilter: starting local_outlier_detection", 2)

    node_info = collect_node_info(tt)

    node_info, z_scale = calculate_node_timings(tt, node_info)
    tt.logger(f"TreeTime.ClockFilter: z-scale {z_scale:1.2f}", 2)

    outliers = flag_outliers(tt, node_info, z_score_threshold, z_scale)

    for n in tt.tree.get_terminals():
        if n.name in outliers:
            n.bad_branch = True

    tt.outliers=None
    if len(outliers):
        outlier_df = pd.DataFrame(outliers).T.loc[:,['avg_date', 'tau', 'z', 'diagnosis']]\
                                .rename(columns={'avg_date':'given_date', 'tau':'apparent_date', 'z':'z-score'})
        tt.logger("Clock_filter.local_filter marked the following outliers", 2, warn=True)
        if tt.verbose>=2:
            print(outlier_df)
        tt.outliers = outlier_df
    return len(outliers)


def flag_outliers(tt, node_info, z_score_threshold, z_scale):
    def add_outlier_info(z, n, n_info, parent_tau, mu):
        n_info['z'] = z
        diagnosis=''
        # muts = n_info["nmuts"] if n.is_terminal() else 0.0
        # parent_tau = node_info[n.up.name]['tau'] if n.is_terminal() else n_info['tau']
        if z<0:
            if np.abs(n_info['avg_date']-parent_tau) > n_info["nmuts"]/mu:
                diagnosis='date_too_early'
            else:
                diagnosis = 'excess_mutations'
        else:
            diagnosis = 'date_too_late'
        n_info['diagnosis'] = diagnosis
        return n_info

    outliers = {}
    mu = tt.clock_model['slope']*tt.data.full_length
    for n in tt.tree.get_terminals():
        if n.up.up is None:
            continue # do not label children of the root as bad -- typically a problem with root choice that will be fixed anyway
        n_info = node_info[n.name]
        parent_tau = node_info[n.up.name]['tau']
        if n_info['exact_date']:
            z = (n_info['avg_date'] - n_info['tau'])/z_scale
            if np.abs(z) > z_score_threshold:
                outliers[n.name] = add_outlier_info(z, n, n_info, parent_tau, mu)
        elif n.raw_date_constraint and len(n.raw_date_constraint):
            zs = [(n_info['tau'] - x)/z_scale for x in n.raw_date_constraint]
            if zs[0]*zs[1]>0 and np.min(np.abs(zs))>z_score_threshold:
                outliers[n.name] = add_outlier_info(zs[0] if np.abs(zs[0])<np.abs(zs[1]) else zs[1],
                                                    n, n_info, parent_tau, mu)

    return outliers

def calculate_node_timings(tt, node_info, eps=0.2):
    mu = tt.clock_model['slope']*tt.data.full_length
    sigma_sq = (3/mu)**2
    tt.logger(f"Clockfilter.calculate_node_timings: mu={mu:1.3e}/y, sigma={3/mu:1.3e}y", 2)
    for n in tt.tree.find_clades(order='postorder'):
        p = node_info[n.name]
        if not p['exact_date'] or p['skip']:
            continue

        if n.is_terminal():
            prefactor = (p["observations"]/sigma_sq + mu**2/(p["nmuts"]+eps))
            p["a"] = (p["avg_date"]/sigma_sq + mu*p["nmuts"]/(p["nmuts"]+eps))/prefactor
        else:
            children = [node_info[c.name] for c in n if (not node_info[c.name]['skip']) and node_info[c.name]['exact_date']]
            if n==tt.tree.root:
                tmp_children_1 = mu*np.sum([(mu*c["a"]-c["nmuts"])/(eps+c["nmuts"]) for c in children])
                tmp_children_2 = mu**2*np.sum([(1-c["b"])/(eps+c["nmuts"]) for c in children])
                prefactor = (p["observations"]/sigma_sq + tmp_children_2)
                p["a"] = (p["observations"]*p["avg_date"]/sigma_sq + tmp_children_1)/prefactor
            else:
                tmp_children_1 = mu*np.sum([(mu*c["a"]-c["nmuts"])/(eps+c["nmuts"]) for c in children])
                tmp_children_2 = mu**2*np.sum([(1-c["b"])/(eps+c["nmuts"]) for c in children])
                prefactor = (p["observations"]/sigma_sq + mu**2/(p["nmuts"]+eps) + tmp_children_2)
                p["a"] = (p["observations"]*p["avg_date"]/sigma_sq + mu*p["nmuts"]/(p["nmuts"]+eps)+tmp_children_1)/prefactor
        p["b"] = mu**2/(p["nmuts"]+eps)/prefactor

    node_info[tt.tree.root.name]["tau"] = node_info[tt.tree.root.name]["a"]

    ## need to deal with tips without exact dates below.
    dev = []
    for n in tt.tree.get_nonterminals(order='preorder'):
        p = node_info[n.name]
        for c in n:
            c_info = node_info[c.name]
            if c_info['skip']:
                c_info['tau']=p['tau']
            else:
                if c_info['exact_date']:
                    c_info["tau"] = c_info["a"] + c_info["b"]*p["tau"]
                else:
                    c_info["tau"] = p["tau"] + c_info['nmuts']/mu
            if c.is_terminal() and c_info['exact_date']:
                dev.append(c_info['avg_date']-c_info['tau'])

    return node_info, np.std(dev)


def collect_node_info(tt, percentile_for_exact_date=90):
    node_info = {}
    aln = tt.aln or False
    if aln and (not tt.sequence_reconstruction):
        tt.infer_ancestral_sequences(infer_gtr=False)
    L = tt.data.full_length

    date_uncertainty = [np.abs(n.raw_date_constraint[1]-n.raw_date_constraint[0])
                            if type(n.raw_date_constraint)!=float else 0.0
                        for n in tt.tree.get_terminals()
                            if n.raw_date_constraint is not None]
    from scipy.stats import scoreatpercentile
    uncertainty_cutoff = scoreatpercentile(date_uncertainty, percentile_for_exact_date)*1.01

    for n in tt.tree.get_nonterminals(order='postorder'):
        parent = {"dates": [], "tips": {}, "skip":False}
        exact_dates = 0
        for c in n:
            if c.is_terminal():
                child = {'skip':False}
                child["nmuts"] = len([m for m in c.mutations if m[-1] in 'ACGT']) if aln \
                                      else np.round(c.branch_length*L)
                if c.raw_date_constraint is None:
                    child['exact_date'] = False
                elif type(c.raw_date_constraint)==float:
                    child['exact_date'] = True
                else:
                    child['exact_date'] = np.abs(c.raw_date_constraint[1]-c.raw_date_constraint[0])<=uncertainty_cutoff

                if child['exact_date']:
                    exact_dates += 1
                    if child["nmuts"]==0:
                        child['skip'] = True
                        parent["tips"][c.name]={'date': np.mean(c.raw_date_constraint),
                                                'exact_date':child['exact_date']}
                    else:
                        child['skip'] = False
                        child['observations'] = 1
                if c.raw_date_constraint is not None:
                    child["avg_date"] = np.mean(c.raw_date_constraint)
                node_info[c.name] = child
            else:
                if node_info[c.name]['exact_date']:
                    exact_dates += 1

        parent['exact_date'] = exact_dates>0

        parent["nmuts"] = len([m for m in n.mutations if m[-1] in 'ACGT']) if aln else np.round(n.branch_length*L)
        d = [v['date'] for v in parent['tips'].values() if v['exact_date']]
        parent["observations"] = len(d)
        parent["avg_date"] = np.mean(d) if len(d) else 0.0
        node_info[n.name] = parent

    return node_info
