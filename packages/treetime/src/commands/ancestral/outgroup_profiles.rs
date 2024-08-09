use crate::commands::ancestral::anc_graph::Node;
use crate::commands::ancestral::subtree_profiles::EPS;
use crate::gtr::gtr::GTR;
use maplit::btreemap;
use ndarray_stats::QuantileExt;

pub fn outgroup_profiles(node: &mut Node, parent: Option<&Node>, seq: &[char], logLH: &mut f64, gtr: &GTR) {
  match parent {
    // Root state gets different treatment
    None => {
      calculate_root_state(node, logLH, gtr);
    }
    // Calculate outgroup_profile and profile for non-root nodes
    Some(parent) => {
      let mut variable_pos = node.subtree_profile_variable.keys().copied().collect::<Vec<_>>();
      let parent_variable_pos = parent.profile_variable.keys();
      variable_pos.extend(parent_variable_pos);

      node.outgroup_profile_variable = btreemap! {};
      node.profile_variable = btreemap! {};
      for pos in variable_pos {
        if let Some(nuc) = seq.get(pos) {
          // For both the subtree and the parent profile, use the variable profile if available, otherwise use the fixed profile
          let stp = node
            .subtree_profile_variable
            .get(&pos)
            .unwrap_or_else(|| &node.subtree_profile_fixed[nuc]);

          let pp = parent
            .profile_variable
            .get(&pos)
            .unwrap_or_else(|| &parent.profile_fixed[nuc]);

          let vec = pp / node.expQt.dot(stp); // this is numerically tricky, need to guard against division by 0
          let vec_norm = vec.sum();
          let outgroup_profile_variable = vec / vec_norm;
          node.outgroup_profile_variable.insert(pos, outgroup_profile_variable);

          // If uncertainty is high, keep this position
          let vec = stp * node.expQt.t().dot(&node.outgroup_profile_variable[&pos]);
          let vec_norm = vec.sum();
          if (vec.max().unwrap() < &((1.0 - EPS) * vec_norm)) || (*nuc != gtr.alphabet.char(vec.argmax().unwrap())) {
            node.profile_variable.insert(pos, vec / vec_norm);
          }
        }
      }

      // Report for fixed positions for each nucleotide
      node.profile_fixed = btreemap! {};
      node.outgroup_profile_fixed = btreemap! {};
      for nuc in "ACGT".chars() {
        let vec = &parent.profile_fixed[&nuc] / node.expQt.dot(&node.subtree_profile_fixed[&nuc]);
        let vec_norm = vec.sum();
        node.outgroup_profile_fixed.insert(nuc, vec / vec_norm);
        node.profile_fixed.insert(
          nuc,
          &node.subtree_profile_fixed[&nuc] * node.expQt.t().dot(&node.outgroup_profile_fixed[&nuc]),
        );
      }
    }
  }
}

fn calculate_root_state(root: &mut Node, logLH: &mut f64, gtr: &GTR) {
  // Multiply the `message_to_parent` at the root with the equilibrium probabilities
  root.profile_variable = btreemap! {};
  let mut inert_nucs = root.nuc_composition.clone();

  // // Variable positions
  // for (&pos, vec) in &root.subtree_profile_variable {
  //   root.profile_variable.insert(pos, vec * &gtr.pi);
  //   let vec_norm = root.profile_variable[&pos].sum();
  //   *root.profile_variable.entry(pos).or_default() /= vec_norm;
  //   *logLH += vec_norm.ln();
  //   if let Some(nuc) = root.seq.get(pos) {
  //     if gtr.alphabet.contains(*nuc) {
  //       let count = inert_nucs.entry(*nuc).or_default();
  //       *count = count.saturating_sub(1);
  //     }
  //   }
  // }
  // 
  // // Fixed positions
  // root.profile_fixed = btreemap! {};
  // for &nuc in gtr.alphabet.chars() {
  //   if let Some(subtree_profile_fixed) = root.subtree_profile_fixed.get(&nuc) {
  //     root.profile_fixed.insert(nuc, subtree_profile_fixed * &gtr.pi);
  //     if let Some(profile_fixed) = root.profile_fixed.get_mut(&nuc) {
  //       let vec_norm = profile_fixed.sum();
  //       *profile_fixed /= vec_norm;
  //       *logLH += (inert_nucs[&nuc] as f64) * vec_norm.ln();
  //     }
  //   }
  // }
}
