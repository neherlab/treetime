from treetime import MissingDataError,UnknownMethodError,NotReadyError
import numpy as np

def initialize_basetree(self, root=None, infer_gtr=True, relaxed_clock=None, n_iqd = None,
            resolve_polytomies=True, max_iter=0, Tc=None, fixed_clock_rate=None,
            time_marginal=False, sequence_marginal=False, branch_length_mode='auto',
            vary_rate=False, use_covariation=False, tracelog_file=None, **kwargs):
    """"
    Start of treetime.run function
    """

    # register the specified covaration mode
    self.use_covariation = use_covariation or (vary_rate and (not type(vary_rate)==float))

    if (self.tree is None) or (self.aln is None and self.data.full_length is None):
        raise MissingDataError("TreeTime.run: ERROR, alignment or tree are missing")
    if self.aln is None:
        branch_length_mode='input'

    self._set_branch_length_mode(branch_length_mode)

    # determine how to reconstruct and sample sequences
    seq_kwargs = {"marginal_sequences":sequence_marginal or (self.branch_length_mode=='marginal'),
                    "branch_length_mode": self.branch_length_mode,
                    "sample_from_profile":"root",
                    "reconstruct_tip_states":kwargs.get("reconstruct_tip_states", False)}

    tt_kwargs = {'clock_rate':fixed_clock_rate, 'time_marginal':False}
    tt_kwargs.update(kwargs)

    seq_LH = 0
    if "fixed_pi" in kwargs:
        seq_kwargs["fixed_pi"] = kwargs["fixed_pi"]
    if "do_marginal" in kwargs:
        time_marginal=kwargs["do_marginal"]

    # initially, infer ancestral sequences and infer gtr model if desired
    if self.branch_length_mode=='input':
        if self.aln:
            self.infer_ancestral_sequences(infer_gtr=infer_gtr, marginal=seq_kwargs["marginal_sequences"], **seq_kwargs)
            self.prune_short_branches()
    else:
        print("else executed")
        self.optimize_tree(infer_gtr=infer_gtr,
                            max_iter=1, prune_short=True, **seq_kwargs)
    avg_root_to_tip = np.mean([x.dist2root for x in self.tree.get_terminals()])

    # optionally reroot the tree either by oldest, best regression or with a specific leaf
    if n_iqd or root=='clock_filter':
        if "plot_rtt" in kwargs and kwargs["plot_rtt"]:
            plot_rtt=True
        else:
            plot_rtt=False
        reroot_mechanism = 'least-squares' if root=='clock_filter' else root
        self.clock_filter(reroot=reroot_mechanism, n_iqd=n_iqd, plot=plot_rtt, fixed_clock_rate=fixed_clock_rate)
        print("clock filter set")
    elif root is not None:
        self.reroot(root=root, clock_rate=fixed_clock_rate)
        print("root")

    if self.branch_length_mode=='input':
        if self.aln:
            self.infer_ancestral_sequences(**seq_kwargs)
            print("infer ancestral seq")
    else:
        self.optimize_tree(max_iter=1, prune_short=False,**seq_kwargs)
        print("optimize tree")

    # infer time tree and optionally resolve polytomies
    self.logger("###TreeTime.run: INITIAL ROUND",0)
    self.make_time_tree(**tt_kwargs)

    if self.aln:
        seq_LH = self.tree.sequence_marginal_LH if seq_kwargs['marginal_sequences'] else self.tree.sequence_joint_LH
        print("self Aln")
    self.LH =[[seq_LH, self.tree.positional_joint_LH, 0]]

    if root is not None and max_iter:
        new_root = self.reroot(root='least-squares' if root=='clock_filter' else root, clock_rate=fixed_clock_rate)
        self.logger("###TreeTime.run: rerunning timetree after rerooting",0)
        self.make_time_tree(**tt_kwargs)

    # iteratively reconstruct ancestral sequences and re-infer
    # time tree to ensure convergence.
    niter = 0
    ndiff = 0

    # Initialize the tracelog dict attribute
    self.trace_run = []
    self.trace_run.append(self.tracelog_run(niter=0, ndiff=0, n_resolved=0,
                            time_marginal = tt_kwargs['time_marginal'],
                            sequence_marginal = seq_kwargs['marginal_sequences'], Tc=None, tracelog=tracelog_file))
    need_new_time_tree=False
    while niter < max_iter:
        self.logger("###TreeTime.run: ITERATION %d out of %d iterations"%(niter+1,max_iter),0)
        # add coalescent prior
        tmpTc=None
        if Tc:
            if Tc=='skyline' and niter<max_iter-1:
                tmpTc='const'
            else:
                tmpTc=Tc
            self.add_coalescent_model(tmpTc, **kwargs)
            need_new_time_tree = True

        # # estimate a relaxed molecular clock
        # if relaxed_clock:
        #     print("relaxed_clock", relaxed_clock)
        #     self.relaxed_clock(**relaxed_clock)
        #     need_new_time_tree = True

        # n_resolved=0
        # if resolve_polytomies:
        #     # if polytomies are found, rerun the entire procedure
        #     n_resolved = self.resolve_polytomies()
        #     if n_resolved:
        #         self.prepare_tree()
        #         if self.branch_length_mode!='input': # otherwise reoptimize branch length while preserving branches without mutations
        #             self.optimize_tree(prune_short=False, max_iter=0, **seq_kwargs)

        #         need_new_time_tree = True

        niter+=1


def add_coalescent(self, FFT, root=None, infer_gtr=True, relaxed_clock=None, n_iqd = None,
            resolve_polytomies=True, max_iter=0, Tc=None, fixed_clock_rate=None,
            time_marginal=False, sequence_marginal=False, branch_length_mode='auto',
            vary_rate=False, use_covariation=False, tracelog_file=None, **kwargs):
    
    niter = 0
    tmpTc=None
    if Tc:
        if Tc=='skyline' and niter<max_iter-1:
            tmpTc='const'
        else:
            tmpTc=Tc
        #self.add_coalescent_model(tmpTc, **kwargs)
        #need_new_time_tree = True
        if FFT:
            from treetime_fft.merger_models import Coalescent
        else:
            from treetime.merger_models import Coalescent
        self.logger('TreeTime.run: adding coalescent prior with Tc='+str(Tc),1)
        self.merger_model = Coalescent(self.tree,
                                       date2dist=self.date2dist, logger=self.logger)

        if Tc=='skyline': # restrict skyline model optimization to last iteration
            self.merger_model.optimize_skyline(**kwargs)
            self.logger("optimized a skyline ", 2)
        else:
            if Tc in ['opt', 'const']:
                self.merger_model.optimize_Tc()
                self.logger("optimized Tc to %f"%self.merger_model.Tc.y[0], 2)
            else:
                try:
                    self.merger_model.set_Tc(Tc)
                except:
                    self.logger("setting of coalescent time scale failed", 1, warn=True)


def multiply(x_vals, dists):
    '''
    multiplies a list of Distribution objects, returning the output at the points given by x_vals
    '''
    from treetime_fft.distribution import Distribution

    if  not all([isinstance(k, Distribution) for k in dists]):
        from treetime.distribution import Distribution as Distribution
        if  not all([isinstance(k, Distribution) for k in dists]):
            raise NotImplementedError("Can only multiply Distribution objects")

    n_delta = np.sum([k.is_delta for k in dists])
    min_width = np.max([k.min_width for k in dists])
    if n_delta>1:
        raise ArithmeticError("Cannot multiply more than one delta functions!")
    elif n_delta==1:
        delta_dist_ii = np.where([k.is_delta for k in dists])[0][0]
        delta_dist = dists[delta_dist_ii]
        new_xpos = delta_dist.peak_pos
        new_weight  = np.prod([k.prob(new_xpos) for k in dists if k!=delta_dist_ii]) * delta_dist.weight
        res = Distribution.delta_function(new_xpos, weight = new_weight,min_width=min_width)
    else:
        #new_xmin = np.max([k.xmin for k in dists])
        #new_xmax = np.min([k.xmax for k in dists])

        #x_vals = np.unique(np.concatenate([k.x for k in dists]))
        x_vals = x_vals
        #print(len(x_vals))
        #x_vals = x_vals[(x_vals>new_xmin-TINY_NUMBER)&(x_vals<new_xmax+TINY_NUMBER)]
        y_vals = np.sum([k.__call__(x_vals) for k in dists], axis=0)
        peak = y_vals.min()
        #ind = (y_vals-peak)<BIG_NUMBER/1000
        #n_points = ind.sum()
        n_points = y_vals.sum()
        if n_points == 0:
            print ("ERROR in distribution multiplication: Distributions do not overlap")
            # x_vals = [0,1]
            # y_vals = [BIG_NUMBER,BIG_NUMBER]
            # res = Distribution(x_vals, y_vals, is_log=True,
            #                     min_width=min_width, kind='linear')
        elif n_points == 1:
            res = Distribution.delta_function(x_vals[0])
        else:
            #res = Distribution(x_vals[ind], y_vals[ind], is_log=True,
            #                    min_width=min_width, kind='linear', assume_sorted=True)
            res = Distribution(x_vals, y_vals, is_log=True,
                                min_width=min_width, kind='linear', assume_sorted=True)

    return res