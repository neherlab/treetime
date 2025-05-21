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

    prepare_gauss_model(tt)
    genome_wide_mu = tt.clock_model['slope']*tt.data.full_length
    initial_params = {'mu':genome_wide_mu, 'sigma':3/genome_wide_mu, 'eps':0.5, 'Tc':None}
    gauss_params = iterate_gauss_model(tt, initial_params)
    tt.logger(f"TreeTime.ClockFilter: z-scale {gauss_params['z_scale']:1.2f}", 2)

    flag_outliers(tt, z_score_threshold, gauss_params['z_scale'])

    outliers = []
    for n in tt.tree.get_terminals():
        if n.gauss_model['outlier']:
            n.bad_branch = True
            outliers.append({'tau':n.gauss_model['tau'], 'avg_date':np.mean(n.gauss_model['observations']), 'obs':n.gauss_model['observations'],
                            'z':n.gauss_model['z'], 'diagnosis':n.gauss_model['diagnosis'], 'name':n.name})

    tt.outliers=None
    if len(outliers):
        outlier_df = pd.DataFrame(outliers).loc[:,['name', 'avg_date', 'tau', 'z', 'diagnosis']]\
                                .rename(columns={'avg_date':'given_date', 'tau':'apparent_date', 'z':'z-score'})
        outlier_df.set_index('name', inplace=True)
        tt.logger("Clock_filter.local_filter marked the following outliers", 2, warn=True)
        if tt.verbose>=2:
            print(outlier_df)
        tt.outliers = outlier_df
    return len(outliers)


def flag_outliers(tt, z_score_threshold, z_scale):
    def add_outlier_info(z, n, parent_tau, mu):
        gm = n.gauss_model
        gm['z'] = z
        diagnosis=''
        # muts = n_info["nmuts"] if n.is_terminal() else 0.0
        # parent_tau = node_info[n.up.name]['tau'] if n.is_terminal() else n_info['tau']
        if z<0:
            if np.abs(np.mean(gm['observations']-parent_tau)) > gm["n_muts"]/mu:
                diagnosis='date_too_early'
            else:
                diagnosis = 'excess_mutations'
        else:
            diagnosis = 'date_too_late'
        gm['diagnosis'] = diagnosis
        if diagnosis:
            gm['outlier'] = True

    mu = tt.clock_model['slope']*tt.data.full_length
    for n in tt.tree.get_terminals():
        gm = n.gauss_model
        if n.up.up is None:
            continue # do not label children of the root as bad -- typically a problem with root choice that will be fixed anyway
        parent_tau = n.up.gauss_model['tau']
        if n.bad_branch:
            add_outlier_info(0, n, parent_tau, mu)
        elif gm['descendent_observations']>0:
            z = (np.mean(gm['observations']) - gm['tau'])/z_scale
            if np.abs(z) > z_score_threshold:
                add_outlier_info(z, n, parent_tau, mu)
        elif n.raw_date_constraint and len(n.raw_date_constraint):
            zs = [(gm['tau'] - x)/z_scale for x in n.raw_date_constraint]
            if zs[0]*zs[1]>0 and np.min(np.abs(zs))>z_score_threshold:
                add_outlier_info(zs[0] if np.abs(zs[0])<np.abs(zs[1]) else zs[1],
                                  n, parent_tau, mu)


def iterate_gauss_model(tt, params):
    from scipy.interpolate import interp1d
    mu = params['mu']
    sigma_sq = params['sigma']**2
    eps = params['eps']
    Tc = params['Tc']
    with_coal = (Tc is not None) and ("lineages" in params)

    # backward pass
    for n in tt.tree.find_clades(order='postorder'):
        gm  = n.gauss_model
        if gm['descendent_observations']==0:
            continue
        if with_coal:
            n_lineages = params['lineages'](gm['tau'])
            coal= n_lineages/Tc*(len(n.clades)-1) + 1/Tc
        else:
            coal = 0

        if n.is_terminal():
            prefactor = (gm["n_observations"]/sigma_sq + mu**2/(gm["n_muts"]+eps))
            gm["a"] = (np.sum(gm["observations"])/sigma_sq + mu*gm["n_muts"]/(gm["n_muts"]+eps) + coal)/prefactor
        else:
            valid_children = [c.gauss_model for c in n.clades if c.gauss_model['descendent_observations']>0]
            tmp_children_1 = mu*np.sum([(mu*c["a"]-c["n_muts"])/(eps+c["n_muts"]) for c in valid_children])
            tmp_children_2 = mu**2*np.sum([(1-c["b"])/(eps+c["n_muts"]) for c in valid_children])
            if n==tt.tree.root:
                prefactor = (gm["n_observations"]/sigma_sq + tmp_children_2)
                gm["a"] = (np.sum(gm['observations'])/sigma_sq + tmp_children_1 + coal)/prefactor
            else:
                prefactor = (gm["n_observations"]/sigma_sq + mu**2/(gm["n_muts"]+eps) + tmp_children_2)
                gm["a"] = (np.sum(gm['observations'])/sigma_sq + mu*gm["n_muts"]/(gm["n_muts"]+eps)+tmp_children_1 + coal)/prefactor
        if tt.tree.root!=n:
            gm["b"] = mu**2/(gm["n_muts"]+eps)/prefactor

    # evaluation of tau at the root
    gm_r = tt.tree.root.gauss_model
    gm_r["tau"] = gm_r["a"]
    n_lineage_list = [(gm_r["tau"], len(tt.tree.root.clades))]
    # forward pass
    tsq = 0
    dt = 0
    dev = []
    logL = np.log(2*np.pi*sigma_sq)*0.5*gm_r["descendent_observations"]
    for n in tt.tree.get_nonterminals(order='preorder'):
        gm = n.gauss_model
        for c in n:
            gm_c = c.gauss_model
            if gm_c['descendent_observations']==0:
                gm_c['tau']=gm['tau'] + gm_c['n_muts']/mu
            else:
                gm_c["tau"] = gm_c["a"] + gm_c["b"]*gm["tau"]
                n_lineage_list.append((gm_c["tau"], len(c.clades)-1))
                tsq += (gm_c["tau"] - gm["tau"])**2/(2*(gm_c["n_muts"]+eps))
                dt += (gm_c["tau"] - gm["tau"])*gm_c["n_muts"]/(2*(gm_c["n_muts"]+eps))
                if gm_c['n_observations']:
                    logL += 0.5*np.sum((gm_c['observations'] - gm_c['tau'])**2)/sigma_sq
                    dev.extend(gm_c['observations'] - gm_c['tau'])
                logL += 0.5*((gm_c["tau"] - gm["tau"])*mu - gm_c["n_muts"])**2/(gm_c["n_muts"]+eps)
        if with_coal:
            n_lineages = params['lineages'](gm['tau'])
            logL -= np.log(n_lineages/Tc)*(len(n.clades)-1)

    n_lineage_list.sort()
    n_lineage_list[-1] = (n_lineage_list[-1][0], 0)
    return {'lineages': interp1d([x[0] for x in n_lineage_list], np.maximum(1, np.cumsum([x[1] for x in n_lineage_list])), kind='previous', bounds_error=False, fill_value=1),
            'mu':dt/tsq, 'dmu':1/tsq**0.5, 'logL':logL, "sigma": params['sigma'], 'Tc':Tc, 'eps':eps, 'z_scale':np.std(dev)}


def prepare_gauss_model(tt, uncertainty_cutoff=0.2):
    aln = tt.aln or False
    L = tt.data.full_length
    for n in tt.tree.find_clades(order='postorder'):
        n.gauss_model = {'a':0, 'b':0, 'tau':0, 'z':0, 'outlier':False}
        n.gauss_model['n_muts'] = len([m for m in n.mutations if m[-1] in 'ACGT']) if aln \
                                      else np.round(n.branch_length*L)
        obs = []
        if n.is_terminal() and (n.raw_date_constraint is not None) and (n.bad_branch is False):
            if isinstance(n.raw_date_constraint, float):
                obs.append(n.raw_date_constraint)
            else: # interval
                if n.raw_date_constraint[1]-n.raw_date_constraint[0]<=uncertainty_cutoff:
                    obs.append(np.mean(n.raw_date_constraint))

        n.gauss_model['observations'] = np.array(obs)
        n.gauss_model['n_observations'] = len(obs)
        # total number of explicit observation of any downstream genotype that isn't filtered
        n.gauss_model['descendent_observations'] = n.gauss_model['n_observations'] + np.sum([c.gauss_model['descendent_observations'] for c in n.clades])


