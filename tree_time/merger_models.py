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
	branch_times = np.sort(np.unique([n.abs_t for n in tree.find_clades()]))
	branch_counts = np.zeros_like(branch_times, dtype=int)
	for n in tree.find_clades():
		if n.up  is None:
			continue
		lw = branch_times.searchsorted(n.abs_t)
		up = branch_times.searchsorted(n.up.abs_t)
		branch_counts[lw:up]+=1

	# calculate the merger rates in each interval
	merger_rates = branch_counts*(branch_counts-1.0)*0.5/Tc

	# assign those rates to all nodes in the tree
	for n in tree.find_clades():
		if n.up  is None:
			continue
		lw = branch_times.searchsorted(n.abs_t)
		up = branch_times.searchsorted(n.up.abs_t)
		if lw==up:
			lw = min(lw,up-1)
		n.merger_rate = merger_rates[lw:up].mean()

