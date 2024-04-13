use crate::commands::ancestral::anc_graph::{AncestralGraph, Edge, Node};
use crate::graph::graph::SafeNode;
use crate::gtr::gtr::GTR;
use crate::seq::find_mixed_sites::MixedSite;
use crate::utils::ndarray::{product_axis, stack_owned};
use itertools::Itertools;
use maplit::btreemap;
use ndarray::{Array2, Axis};
use ndarray_stats::QuantileExt;
use std::collections::BTreeSet;

// threshold to keep position among the variable ones
pub const EPS: f64 = 1e-6;

// Function to calculate the likelihood for subtree profiles
pub fn subtree_profiles(
  graph: &AncestralGraph,
  node: &mut Node,
  edge: Option<&Edge>,
  children: &[SafeNode<Node>], // TODO: Should we lock all these child locks in advance?
  seq: &[char],
  gtr: &GTR,
  logLH: &mut f64,
) {
  let prof_nuc = &gtr.alphabet.profile_map;

  // GTR matrix associated with this branch length. Using 0 length for the root saves an extra calculation below
  let t = edge.map(|edge| edge.weight).unwrap_or_default();
  node.expQt = gtr.expQt(t).t().to_owned(); // TODO: might make sense to save this on the edge

  // We have calculated the total nucleotide composition in the sequence representation.
  // From this, we will subtract positions that are tracked as variable positions -- hence need to copy
  let mut fixed_nuc_count = node.nuc_composition.clone();

  // Each node will get a vector with the probability distribution of non-variable positions.
  // This vector should always be peaked around the focal nucleotide and quantifies the uncertainty around it.
  node.subtree_profile_fixed = btreemap! {};

  if node.is_leaf() {
    // For terminal nodes, we deem mixed sites and sites that have mutations in the branch to the parent as variable
    node.subtree_profile_variable = btreemap! {};
    for MixedSite { pos, nuc } in &node.mixed {
      node.subtree_profile_variable.insert(*pos, prof_nuc[nuc].clone());
      let count = fixed_nuc_count.entry(seq[*pos]).or_insert(0);
      *count = count.saturating_sub(1);
    }
    for (pos, (anc, der)) in &node.mutations {
      node.subtree_profile_variable.insert(*pos, prof_nuc[der].clone());
      let count = fixed_nuc_count.entry(seq[*pos]).or_insert(0);
      *count = count.saturating_sub(1);
    }

    // This could be done more efficiently. We just need to look-up these positions, no need to save the flat vector.
    for rg in &node.undetermined {
      for pos in rg.0..rg.1 {
        node
          .subtree_profile_variable
          .insert(pos, prof_nuc[&gtr.alphabet.ambiguous].clone());
      }
    }

    for &n in gtr.alphabet.chars() {
      node.subtree_profile_fixed.insert(n, prof_nuc[&n].clone());
    }
  } else {
    // For internal nodes, we consider all positions that are variable in any of the children
    node.subtree_profile_variable = btreemap! {};
    let variable_positions: BTreeSet<usize> = children
      .iter()
      .flat_map(|c| {
        let c = c.read_arc().payload().read_arc();
        c.subtree_profile_variable.keys().copied().collect_vec()
      })
      .collect();

    for pos in variable_positions {
      let nuc = seq[pos];
      // if the parsimony sequence is ambiguous there is no information anywhere, can skip
      // TODO: if there is no gap in the alphabet, we can skip as well with gap.
      if nuc == gtr.alphabet.ambiguous {
        continue;
      }

      // Calculate the product of child messages
      let child_messages = children
        .iter()
        .map(|child| {
          let child = child.read_arc().payload().read_arc();

          let prof = child
            .subtree_profile_variable
            .get(&pos)
            .unwrap_or(&child.subtree_profile_fixed[&nuc]);

          child.expQt.dot(prof)
        })
        .collect_vec();
      let child_messages: Array2<f64> = stack_owned(Axis(0), &child_messages).unwrap();
      let vec = product_axis(&child_messages, Axis(0));
      let vec_norm = vec.sum();
      *logLH += vec_norm.ln();

      // Add position to variable states if the subleading states have a probability exceeding EPS
      if (vec.max().unwrap() < &((1.0 - EPS) * vec_norm)) || (nuc != gtr.alphabet.char(vec.argmax().unwrap())) {
        node.subtree_profile_variable.insert(pos, vec / vec_norm);
      }

      // This position is accounted for, hence we can subtract it from the count of fixed nucs
      let count = fixed_nuc_count.entry(nuc).or_insert(0);
      *count = count.saturating_sub(1);
    }

    // Collect contribution from the inert sites
    for &nuc in gtr.alphabet.chars() {
      if !fixed_nuc_count.contains_key(&nuc) {
        continue;
      }
      let child_messages = children
        .iter()
        .map(|child| {
          let child = child.read_arc().payload().read_arc();
          let prof = &child.subtree_profile_fixed[&nuc];
          child.expQt.dot(prof)
        })
        .collect_vec();
      let child_messages: Array2<f64> = stack_owned(Axis(0), &child_messages).unwrap();
      let vec = product_axis(&child_messages, Axis(0));
      let vec_norm = vec.sum();
      *logLH += (fixed_nuc_count[&nuc] as f64) * vec_norm.ln();
      node.subtree_profile_fixed.insert(nuc, vec / vec_norm);
    }

    // Add position that mutate towards the parent
    for (pos, (anc, der)) in &node.mutations {
      if !node.subtree_profile_variable.contains_key(pos) {
        let prof = node.subtree_profile_fixed[der].to_owned();
        node.subtree_profile_variable.insert(*pos, prof);
      }
    }
  }
  // TODO: we could save c.expQT.dot(xxx) on the edges. that would save some computation.
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::anc_graph::{Edge, Node};
  use crate::commands::ancestral::outgroup_profiles::outgroup_profiles;
  use crate::graph::create_graph_from_nwk::create_graph_from_nwk_str;
  use crate::gtr::gtr::GTRParams;
  use crate::io::json::{json_stringify, JsonPretty};
  use crate::seq::representation::{compress_sequences, post_order_intrusive, pre_order_intrusive};
  use crate::utils::random::get_random_number_generator;
  use crate::utils::string::vec_to_string;
  use crate::{o, pretty_assert_ulps_eq};
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

  #[rstest]
  fn test_subtree_profiles_traversal() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let mut rng = get_random_number_generator(Some(42));

    let inputs = BTreeMap::from([
      (o!("A"), o!("ACATCGCCNNA--G")),
      (o!("B"), o!("GCATCCCTGTA-NG")),
      (o!("C"), o!("CCGGCGATGTATTG")),
      (o!("D"), o!("TCGGCCGTGTRTTG")),
    ]);

    let L = inputs.first_key_value().unwrap().1.len();

    let graph = create_graph_from_nwk_str::<Node, Edge>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    compress_sequences(&inputs, &graph, &mut rng).unwrap();

    let root = graph.get_exactly_one_root()?;
    let mut root_seq = { root.read_arc().payload().read_arc().seq.clone() };

    let alphabet = Alphabet::new(AlphabetName::NucNogap)?;

    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.2, 0.3, 0.15, 0.45],
    })?;

    let mut logLH = 0.0;

    post_order_intrusive(&graph, &root, None, &mut root_seq, &mut |node, edge, seq| {
      let node = node.write_arc();

      let children: Vec<SafeNode<Node>> = graph
        .children_of(&node)
        .into_iter()
        .map(|(child, _)| child)
        .collect_vec();

      let mut node = node.payload().write_arc();
      let edge = edge.map(|e| e.payload().read_arc());
      subtree_profiles(&graph, &mut node, edge.as_deref(), &children, seq, &gtr, &mut logLH);
    });

    pretty_assert_ulps_eq!(-36.73309018328223, logLH, epsilon = 1e-5);

    pre_order_intrusive(&graph, &root, &mut root_seq, &mut |node, seq| {
      let node = node.write_arc();

      let parent = graph
        .one_parent_of(&node)
        .unwrap()
        .map(|parent| parent.read_arc().payload().read_arc());

      let mut node = node.payload().write_arc();
      outgroup_profiles(&mut node, parent.as_deref(), seq, &mut logLH, &gtr);
    });

    pretty_assert_ulps_eq!(-57.189205994979055, logLH, epsilon = 1e-5);



    Ok(())
  }
}
