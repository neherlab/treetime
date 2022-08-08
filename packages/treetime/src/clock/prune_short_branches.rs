/// If the branch length is less than the minimal value, remove the branch
/// from the tree.
///
/// **Requires** ancestral sequence reconstruction
pub fn prune_short_branches() {
  // for node in self.tree.find_clades():
  //     if node.up is None or node.is_terminal():
  //         continue
  //     # probability of the two seqs separated by zero time is not zero
  //     if  ((node.branch_length<0.1*self.one_mutation) and
  //          (self.gtr.prob_t(node.up._cseq, node._cseq, 0.0,
  //                           pattern_multiplicity=self.data.multiplicity(mask=node.mask)) > 0.1)):
  //         # re-assign the node children directly to its parent
  //         node.up.clades = [k for k in node.up.clades if k != node] + node.clades
  //         for clade in node.clades:
  //             clade.up = node.up

  unimplemented!("prune_short_branches: not yet implemented")
}
