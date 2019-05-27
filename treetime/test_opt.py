    def optimize_tree_marginal_new(self, damping=0.5):
        L = self.data.compressed_length
        n_states = self.gtr.alphabet.shape[0]
        # propagate leaves --> root, set the marginal-likelihood messages
        for node in self.tree.find_clades(order='postorder'): #leaves -> root
            if node.up is None and len(node.clades)==2:
                continue

            profiles = [c.marginal_subtree_LH for c in node] + [node.marginal_outgroup_LH]
            bls = [c.branch_length for c in nodes] + [node.branch_length]
            new_bls = self.optimize_star(profiles,bls, last_is_root=node.up is None)

            # regardless of what was before, set the profile to ones
            tmp_log_subtree_LH = np.zeros((L,n_states), dtype=float)
            node.marginal_subtree_LH_prefactor = np.zeros(L, dtype=float)
            for ch in ci,node.clades:
                ch.branch_length = new_bls[ci]
                ch.marginal_log_Lx = self.gtr.propagate_profile(ch.marginal_subtree_LH,
                                                                ch.branch_length, return_log=True)
                tmp_log_subtree_LH += ch.marginal_log_Lx
                node.marginal_subtree_LH_prefactor += ch.marginal_subtree_LH_prefactor

            node.marginal_subtree_LH, offset = normalize_profile(tmp_log_subtree_LH, log=True)
            node.marginal_subtree_LH_prefactor += offset # and store log-prefactor

            if node.up:
                node.marginal_log_Lx = self.gtr.propagate_profile(node.marginal_subtree_LH,
                                                node.branch_length, return_log=True) # raw prob to transfer prob up
                tmp_msg_from_parent = self.gtr.evolve(node.marginal_outgroup_LH,
                                                 self._branch_length_to_gtr(node), return_log=False)
                node.marginal_profile, pre = normalize_profile(node.marginal_subtree_LH * tmp_msg_from_parent, return_offset=False)
            else:
                node.marginal_profile, pre = normalize_profile(node.marginal_subtree_LH * node.marginal_outgroup_LH, return_offset=False)

        root=self.tree.root
        print(len(root.clades))
        if len(root.clades)==2:
            tmp_log_subtree_LH = np.zeros((L,n_states), dtype=float)
            root.marginal_subtree_LH_prefactor = np.zeros(L, dtype=float)
            old_bl = root.clades[0].branch_length + root.clades[1]
            bl = self.gtr.optimal_t_compressed((root.clades[0].marginal_subtree_LH*root.marginal_outgroup_LH,
                                                root.clades[1].marginal_subtree_LH), multiplicity=self.data.multiplicity,
                                                profiles=True, tol=1e-8)
            for ch in root:
                ch.branch_length *= ((1-damping)*old_bl + damping*bl)/old_bl
                ch.marginal_log_Lx = self.gtr.propagate_profile(ch.marginal_subtree_LH,
                                            ch.branch_length, return_log=True) # raw prob to transfer prob up
                tmp_log_subtree_LH += ch.marginal_log_Lx
                root.marginal_subtree_LH_prefactor += ch.marginal_subtree_LH_prefactor

            root.marginal_subtree_LH, offset = normalize_profile(tmp_log_subtree_LH, log=True)
            root.marginal_subtree_LH_prefactor += offset # and store log-prefactor


        self.total_LH_and_root_sequence(assign_sequence=False)
        self.preorder_traversal_marginal(assign_sequence=False, reconstruct_leaves=False)




    def optimize_tree_marginal_new2(self, n_iter_internal=2, damping=0.5):
        L = self.data.compressed_length
        n_states = self.gtr.alphabet.shape[0]
        # propagate leaves --> root, set the marginal-likelihood messages
        for node in self.tree.get_nonterminals(order='postorder'): #leaves -> root
            if node.up is None and len(node.clades)==2:
                continue
            # regardless of what was before, set the profile to ones
            for ii in range(n_iter_internal):
                damp = damping**(1+ii)
                tmp_log_subtree_LH = np.zeros((L,n_states), dtype=float)
                node.marginal_subtree_LH_prefactor = np.zeros(L, dtype=float)
                for ch in node.clades:
                    outgroup = np.exp(np.log(np.maximum(ttconf.TINY_NUMBER, node.marginal_profile)) - ch.marginal_log_Lx)

                    bl = self.gtr.optimal_t_compressed((ch.marginal_subtree_LH, outgroup), multiplicity=self.data.multiplicity, profiles=True, tol=1e-8)
                    new_bl = (1-damp)*bl + damp*ch.branch_length
                    ch.branch_length=new_bl
                    ch.marginal_log_Lx = self.gtr.propagate_profile(ch.marginal_subtree_LH,
                                                new_bl, return_log=True) # raw prob to transfer prob up
                    tmp_log_subtree_LH += ch.marginal_log_Lx
                    node.marginal_subtree_LH_prefactor += ch.marginal_subtree_LH_prefactor

                node.marginal_subtree_LH, offset = normalize_profile(tmp_log_subtree_LH, log=True)
                node.marginal_subtree_LH_prefactor += offset # and store log-prefactor

                if node.up:
                    bl = self.gtr.optimal_t_compressed((node.marginal_subtree_LH, node.marginal_outgroup_LH), multiplicity=self.data.multiplicity, profiles=True, tol=1e-8)
                    new_bl = (1-damp)*bl + damp*node.branch_length
                    node.branch_length=new_bl
                    node.marginal_log_Lx = self.gtr.propagate_profile(node.marginal_subtree_LH,
                                                    new_bl, return_log=True) # raw prob to transfer prob up
                    node.marginal_outgroup_LH, pre = normalize_profile(np.log(np.maximum(ttconf.TINY_NUMBER, node.up.marginal_profile)) - node.marginal_log_Lx,
                                                 log=True, return_offset=False)

                    tmp_msg_from_parent = self.gtr.evolve(node.marginal_outgroup_LH,
                                                     self._branch_length_to_gtr(node), return_log=False)
                    node.marginal_profile, pre = normalize_profile(node.marginal_subtree_LH * tmp_msg_from_parent, return_offset=False)
                else:
                    node.marginal_profile, pre = normalize_profile(node.marginal_subtree_LH * node.marginal_outgroup_LH, return_offset=False)


        root=self.tree.root
        print(len(root.clades))
        if len(root.clades)==2:
            tmp_log_subtree_LH = np.zeros((L,n_states), dtype=float)
            root.marginal_subtree_LH_prefactor = np.zeros(L, dtype=float)
            old_bl = root.clades[0].branch_length + root.clades[1]
            bl = self.gtr.optimal_t_compressed((root.clades[0].marginal_subtree_LH*root.marginal_outgroup_LH,
                                                root.clades[1].marginal_subtree_LH), multiplicity=self.data.multiplicity,
                                                profiles=True, tol=1e-8)
            for ch in root:
                ch.branch_length *= bl/old_bl
                ch.marginal_log_Lx = self.gtr.propagate_profile(ch.marginal_subtree_LH,
                                            ch.branch_length, return_log=True) # raw prob to transfer prob up
                tmp_log_subtree_LH += ch.marginal_log_Lx
                root.marginal_subtree_LH_prefactor += ch.marginal_subtree_LH_prefactor

            root.marginal_subtree_LH, offset = normalize_profile(tmp_log_subtree_LH, log=True)
            root.marginal_subtree_LH_prefactor += offset # and store log-prefactor


        self.total_LH_and_root_sequence(assign_sequence=False)
        self.preorder_traversal_marginal(assign_sequence=False, reconstruct_leaves=False)
