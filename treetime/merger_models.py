"""
methods to calculate merger models for a time tree
"""
from __future__ import print_function, division
import numpy as np
from Bio import AlignIO, Phylo


def coalescent(tree, Tc=None):
	'''
	assigns coalescent merger rates to all branches in the tree
	'''
	# determine the time intervals with constant merger rates
	branch_times = np.sort(np.unique([n.time_before_present for n in tree.find_clades()]))
	branch_counts = np.zeros_like(branch_times, dtype=int)
	for n in tree.find_clades():
		if n.up  is None:
			continue
		lw = branch_times.searchsorted(n.time_before_present)
		up = branch_times.searchsorted(n.up.time_before_present)
		branch_counts[lw:up]+=1

	# calculate the merger rates in each interval
	merger_rates = (branch_counts-1.0)*0.5/Tc

	# assign those rates to all nodes in the tree
	for n in tree.find_clades():
		if n.up  is None:
			continue
		lw = branch_times.searchsorted(n.time_before_present)
		up = branch_times.searchsorted(n.up.time_before_present)
		if lw==up:
			lw = min(lw,up-1)
		n.branch_length_interpolator.merger_rate = merger_rates[lw:up].mean()

def traveling_wave(tree, Tc=None, tau=None):
	'''
	assigns coalescent merger rates to all branches in the tree
	'''
	#
	if tau is None:
		tau = Tc/8.0   # 8 is roughly the factor between total coalescence and expansion rate in realistic populations
	# determine the msg to parents
	for n in tree.find_clades(order='postorder'):
		n._polarizer_to_parent = np.sum([c._polarizer_to_parent for c in n.clades])
		n._polarizer_to_parent*= np.exp(-n.branch_length/tau)
		n._polarizer_to_parent+=(1-np.exp(-n.branch_length/tau))*tau

	# determine the msg to children
	tree.root._polarizer_from_parent = 0.0
	for n in tree.get_nonterminals(order='preorder'):
		tmp = np.sum([c._polarizer_to_parent for c in n.clades]) + n._polarizer_from_parent
		for c in n.clades:
			c._polarizer_from_parent = tmp-c._polarizer_to_parent
			c._polarizer_from_parent*= np.exp(-c.branch_length/tau)
			c._polarizer_from_parent+=(1-np.exp(-c.branch_length/tau))*tau

	# determine the msg to parents
	for n in tree.find_clades(order='postorder'):
		n.lbi = n._polarizer_from_parent + np.sum([c._polarizer_to_parent for c in n.clades])

	# assign those rates to all nodes in the tree
	for n in tree.find_clades():
		n.merger_rate = n.lbi/tau/Tc

