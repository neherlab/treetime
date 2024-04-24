#![allow(clippy::default_trait_access)]
use crate::commands::ancestral::anc_graph::{Edge, Node};
use crate::graph::graph::{Graph, GraphNodeBackward, SafeEdge, SafeEdgeRef, SafeNode};
use crate::graph::node::{GraphNodeKey, Named};
use crate::seq::find_char_ranges::{find_ambiguous_ranges, find_gap_ranges};
use crate::seq::find_mixed_sites::{find_mixed_sites, MixedSite};
use crate::seq::range::range_contains;
use crate::seq::range_intersection::range_intersection_iter;
use crate::seq::range_union::range_union;
use crate::seq::sets::{sets_intersection, sets_union};
use crate::utils::manyzip::Manyzip;
use crate::utils::random::random_pop;
use crate::utils::string::vec_to_string;
use crate::{make_error, make_internal_report};
use eyre::{Report, WrapErr};
use itertools::{izip, Itertools};
use maplit::{btreemap, btreeset};
use parking_lot::RwLock;
use rand::Rng;
use std::collections::{BTreeMap, BTreeSet};
use std::ops::DerefMut;

pub fn get_common_length<K: AsRef<str>, V: AsRef<str>>(seqs: impl Iterator<Item = (K, V)>) -> Result<usize, Report> {
  let lengths = seqs
    .into_group_map_by(|(name, seq)| seq.as_ref().len())
    .into_iter()
    .collect_vec();

  match lengths[..] {
    [] => make_error!("No sequences found"),
    [(length, _)] => Ok(length),
    _ => {
      let message = lengths
        .into_iter()
        .sorted_by_key(|(length, _)| *length)
        .map(|(length, entries)| {
          let names = entries
            .iter()
            .map(|(name, _)| format!("    \"{}\"", name.as_ref()))
            .join("\n");
          format!("Length {length}:\n{names}")
        })
        .join("\n\n");

      make_error!("Sequences are expected to all have the same length, but found the following lengths:\n\n{message}")
    }
  }
  .wrap_err("When calculating length of sequences")
}

pub fn compress_sequences(
  seqs: &BTreeMap<String, String>,
  graph: &Graph<Node, Edge>,
  rng: &mut (impl Rng + Send + Sync + Clone),
) -> Result<(), Report> {
  let L = get_common_length(seqs.iter())?;

  graph.iter_depth_first_postorder_forward(
    |GraphNodeBackward {
       is_root,
       is_leaf,
       key,
       payload: mut n,
       children,
     }| {
      if n.is_leaf() {
        // At each terminal node, temporarily store the sequence and ranges of N, - and mixed sites
        n.seq = seqs[&n.name].chars().collect();

        n.gaps = find_gap_ranges(&n.seq);
        n.ambiguous = find_ambiguous_ranges(&n.seq);
        n.undetermined = range_union(&[n.gaps.clone(), n.ambiguous.clone()]); // TODO: avoid copy

        // n.mixed stores the exact character at each mixed positions, the non_consensus stores the possible states
        (n.mixed, n.non_consensus) = find_mixed_sites(&n.seq);
      } else {
        // Positions that are N or - in all children are still N or - in the parent
        let mut children = children.iter().map(|(node, _)| node.write()).collect_vec();
        let mut children = children.iter_mut().map(DerefMut::deref_mut).collect_vec();
        n.undetermined = range_intersection_iter(children.iter().map(|c| &c.undetermined)).collect();
        calculate_fitch_parsimony_in_place(&mut n, &mut children, L);

        // Deallocate full sequences from children
        children.iter_mut().for_each(|c| c.seq = vec![]);
      }
    },
  );

  root_seq_fill_non_consensus_inplace(graph, rng);

  gather_mutations_inplace(graph, rng);

  Ok(())
}

fn calculate_fitch_parsimony_in_place(n: &mut Node, children: &mut [&mut Node], L: usize) {
  // All sites that are not N or - but not fixed will need special treatment
  let non_consensus_positions: BTreeSet<usize> =
    children.iter().flat_map(|c| c.non_consensus.keys().copied()).collect();

  n.non_consensus = BTreeMap::new();
  n.seq = vec!['?'; L];

  // Introduce Ns to mark indeterminate positions
  for &(start, end) in &n.undetermined {
    for pos in start..end {
      n.seq[pos] = 'N';
    }
  }

  // Indeterminate in at least one child
  for pos in non_consensus_positions {
    let child_state_sets = gather_child_state_sets(children, pos);
    let intersection = sets_intersection(child_state_sets.clone().into_iter());
    match intersection.into_iter().collect_vec().as_slice() {
      [state] => {
        // One element in the intersection of child states. Write to sequence.
        n.seq[pos] = *state;
      }
      [] => {
        // Empty intersection of child states. Propagate union of child states.
        let states = sets_union(child_state_sets.into_iter()).into_iter().collect();
        n.non_consensus.insert(pos, states);

        n.seq[pos] = '~';
      }
      states => {
        // More than one element in the intersection of child states. Propagate intersection.
        let states = states.iter().copied().collect();
        n.non_consensus.insert(pos, states);

        n.seq[pos] = '~';
      }
    };
  }

  // Gather state sets for each position across child sequences
  let child_state_sets = Manyzip(children.iter().map(|c| c.seq.iter().copied()).collect_vec()).collect_vec();

  // Zip these states with node sequence
  let state_zip = n.seq.iter_mut().zip(child_state_sets);
  for (pos, (nuc, child_states)) in state_zip.enumerate() {
    if *nuc != '?' {
      continue;
    }

    let child_states_good = child_states
      .iter()
      .filter(|x| ['A', 'C', 'G', 'T'].contains(x))
      .unique()
      .collect_vec();

    match child_states_good.into_iter().copied().collect_vec().as_slice() {
      [state] => {
        // All children have the same state
        *nuc = *state;
      }
      [] => {
        // No child states. Impossible
        unreachable!("No child states. This is impossible");
      }
      states => {
        // Child states differ
        let states = states.iter().copied().collect::<BTreeSet<char>>(); // TODO: avoid copy
        n.non_consensus.insert(pos, states);
        *nuc = '~';
        // memorize child state to assign mutations
        for (c, cstate) in izip!(children.iter_mut(), child_states) {
          if !c.non_consensus.contains_key(&pos) && ['A', 'C', 'G', 'T'].contains(&cstate) {
            c.non_consensus.insert(pos, btreeset! {cstate});
          }
        }
      }
    }
  }
}

#[allow(clippy::map_entry)]
pub fn gather_child_state_sets(children: &mut [&mut Node], pos: usize) -> Vec<BTreeSet<char>> {
  let mut child_state_sets = vec![];
  for c in children {
    if c.non_consensus.contains_key(&pos) {
      child_state_sets.push(c.non_consensus[&pos].clone());
    } else if "ACGT".contains(c.seq[pos]) {
      // memorize child state to assign mutations
      let state = btreeset! {c.seq[pos]};
      child_state_sets.push(state.clone());
      c.non_consensus.insert(pos, state);
    }
  }
  child_state_sets
}

pub fn gather_consensus_child_states(children: &[&mut Node], pos: usize) -> Vec<char> {
  children
    .iter()
    .filter(|c| !range_contains(&c.undetermined, pos))
    .map(|c| c.seq[pos])
    .unique()
    .collect_vec()
}

fn root_seq_fill_non_consensus_inplace(graph: &Graph<Node, Edge>, rng: &mut impl Rng) {
  let roots = graph.get_root_payloads().collect_vec();
  if roots.len() > 1 {
    unimplemented!("Multiple roots are not supported yet");
  }

  let root = &mut *roots[0].write();
  root.non_consensus.iter_mut().for_each(|(pos, states)| {
    root.seq[*pos] = random_pop(states, rng);
  });

  root.nuc_composition = root
    .seq
    .iter()
    .counts()
    .into_iter()
    .map(|(letter, count)| (*letter, count))
    .collect();
}

fn gather_mutations_inplace(graph: &Graph<Node, Edge>, rng: &mut (impl Rng + Send + Sync + Clone)) {
  graph.iter_depth_first_preorder_forward_2(|node| {
    let node = node.write_arc();
    if node.is_leaf() {
      return;
    }

    let children = graph.children_of(&node);
    let mut node = node.payload().write_arc();

    children.into_iter().for_each(|(child, _)| {
      let mut child = child.write_arc().payload().write_arc();
      child.mutations = BTreeMap::new();
      child.nuc_composition = node.nuc_composition.clone();

      // we need this temporary sequence only for internal nodes
      if !child.is_leaf() {
        child.seq = node.seq.clone();
      }

      // all positions that potentially differ in the child c from the parent n are in `c.non_consensus`
      let child_states = child
        .non_consensus
        .iter_mut()
        .filter_map(|(&pos, states)| {
          let parent_state = node.seq[pos];
          (!states.contains(&parent_state)).then_some({
            // in this case we need a mutation to one state in states
            let child_state = random_pop(states, rng);
            (pos, parent_state, child_state)
          })
        })
        .collect_vec();

      child_states.into_iter().for_each(|(pos, parent_state, child_state)| {
        if !child.is_leaf() {
          child.seq[pos] = child_state;
        }
        child.mutations.insert(pos, (parent_state, child_state));
        *child.nuc_composition.entry(parent_state).or_insert(0) -= 1;
        *child.nuc_composition.entry(child_state).or_insert(0) += 1;
      });
    });

    if !node.is_root() {
      node.non_consensus = btreemap! {};
      node.seq = vec![];
    }
  });
}

pub fn reconstruct_leaf_sequences(graph: &Graph<Node, Edge>) -> Result<BTreeMap<String, String>, Report> {
  let root_seq = {
    let root = graph.get_exactly_one_root()?.read_arc().payload().read_arc();
    root.seq.clone()
  };

  graph
    .get_leaves()
    .into_iter()
    .map(|leaf| {
      let leaf = leaf.read_arc();
      let name = leaf.payload().read_arc().name().to_owned();
      let seq = decompress_leaf_sequence(graph, leaf.key(), &root_seq)?;
      Ok((name, seq))
    })
    .collect()
}

pub fn decompress_leaf_sequence(
  graph: &Graph<Node, Edge>,
  node_key: GraphNodeKey,
  root_seq: &[char],
) -> Result<String, Report> {
  let mut seq = root_seq.to_vec();
  let path = graph.path_from_root_to_node(node_key)?;

  for anc in path {
    let anc = anc.read_arc().payload().read_arc();
    // Apply mutations
    for (&pos, &(_, qry_seq)) in &anc.mutations {
      seq[pos] = qry_seq;
    }
  }

  let node = graph
    .get_node(node_key)
    .ok_or_else(|| make_internal_report!("No path from root to node '{node_key}'"))?
    .read_arc()
    .payload()
    .read_arc();

  // introduce N, gaps, and mixed sites
  for &(from, to) in &node.ambiguous {
    seq[from..to].iter_mut().for_each(|x| *x = 'N');
  }
  for &(from, to) in &node.gaps {
    seq[from..to].iter_mut().for_each(|x| *x = '-');
  }
  for &MixedSite { pos, nuc } in &node.mixed {
    seq[pos] = nuc;
  }

  Ok(vec_to_string(seq))
}

pub fn pre_order_intrusive<F>(graph: &Graph<Node, Edge>, node_arc: &SafeNode<Node>, seq: &mut [char], visitor: &mut F)
where
  F: FnMut(&SafeNode<Node>, &[char]),
{
  let node = node_arc.read_arc();

  for (pos, (_, der)) in &node.payload().read_arc().mutations {
    seq[*pos] = *der;
  }

  drop(node); // Prevents deadlock if visitors lock fof writing (likely)
  visitor(node_arc, seq);
  let node = node_arc.read_arc();

  let children = graph.children_of(&node).into_iter().map(|(child, _)| child);
  children.for_each(|child| pre_order_intrusive(graph, &child, seq, visitor));

  for (pos, (anc, _)) in &node.payload().read_arc().mutations {
    seq[*pos] = *anc;
  }
}

pub fn post_order_intrusive<F>(
  graph: &Graph<Node, Edge>,
  node_arc: &SafeNode<Node>,
  edge_arc: Option<&SafeEdge<Edge>>,
  seq: &mut [char],
  visitor: &mut F,
) where
  F: FnMut(&SafeNode<Node>, Option<&SafeEdgeRef<Edge>>, &[char]),
{
  let node = node_arc.read_arc();

  for (pos, (_, der)) in &node.payload().read_arc().mutations {
    seq[*pos] = *der;
  }

  graph
    .children_of(&node)
    .into_iter()
    .for_each(|(child, edge)| post_order_intrusive(graph, &child, Some(&edge), seq, visitor));

  drop(node); // Prevents deadlock if visitors lock fof writing (likely)
  visitor(node_arc, edge_arc.map(RwLock::read_arc).as_ref(), seq);
  let node = node_arc.read_arc();

  for (pos, (anc, _)) in &node.payload().read_arc().mutations {
    seq[*pos] = *anc;
  }
}

#[allow(dead_code)]
pub fn apply_muts_inplace(node: &Node, seq: &mut [char]) {
  for (&pos, &(_, der)) in &node.mutations {
    seq[pos] = der;
  }
}

pub fn apply_non_nuc_changes_inplace(node: &Node, seq: &mut [char]) {
  for &(from, to) in &node.ambiguous {
    seq[from..to].iter_mut().for_each(|x| *x = 'N');
  }
  for &(from, to) in &node.gaps {
    seq[from..to].iter_mut().for_each(|x| *x = '-');
  }
  for &MixedSite { pos, nuc } in &node.mixed {
    seq[pos] = nuc;
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::commands::ancestral::anc_reconstruction_fitch::ancestral_reconstruction_fitch;
  use crate::graph::create_graph_from_nwk::create_graph_from_nwk_str;
  use crate::graph::node::NodeType;
  use crate::io::json::{json_stringify, JsonPretty};
  use crate::o;
  use crate::utils::random::get_random_number_generator;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::Arc;

  #[rstest]
  fn test_seq_representation() -> Result<(), Report> {
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

    let mut actual: Vec<Node> = vec![];
    graph.get_nodes().iter().for_each(|node| {
      actual.push(node.read().payload().read().clone());
    });

    let expected = vec![
      Node {
        name: o!("root"),
        node_type: NodeType::Root(o!("root")),
        mutations: BTreeMap::from([]),
        gaps: vec![],
        ambiguous: vec![],
        undetermined: vec![],
        mixed: vec![],
        non_consensus: BTreeMap::from([
          (0, btreeset!['A', 'C', 'G']),
          (2, btreeset!['A']),
          (3, btreeset!['T']),
          (5, btreeset!['C']),
          (6, btreeset!['A', 'C']),
        ]),
        nuc_composition: btreemap! {
          'A' => 1,
          'C' => 2,
          'G' => 6,
          'T' => 5,
        },
        seq: vec!['T', 'C', 'G', 'G', 'C', 'G', 'G', 'T', 'G', 'T', 'A', 'T', 'T', 'G'],
        subtree_profile_variable: Default::default(),
        subtree_profile_fixed: Default::default(),
        profile_variable: Default::default(),
        profile_fixed: Default::default(),
        outgroup_profile_variable: Default::default(),
        outgroup_profile_fixed: Default::default(),
        expQt: Default::default(),
      },
      Node {
        name: o!("AB"),
        node_type: NodeType::Internal(o!("AB")),
        mutations: BTreeMap::from([
          (0, ('T', 'G')), //
          (2, ('G', 'A')), //
          (3, ('G', 'T')), //
          (6, ('G', 'C')), //
        ]),
        gaps: vec![],
        ambiguous: vec![],
        undetermined: vec![(11, 13)],
        mixed: vec![],
        non_consensus: BTreeMap::from([]),
        nuc_composition: btreemap! {
          'A' => 2,
          'C' => 3,
          'G' => 4,
          'T' => 5,
        },
        seq: vec![],
        subtree_profile_variable: Default::default(),
        subtree_profile_fixed: Default::default(),
        profile_variable: Default::default(),
        profile_fixed: Default::default(),
        outgroup_profile_variable: Default::default(),
        outgroup_profile_fixed: Default::default(),
        expQt: Default::default(),
      },
      Node {
        name: o!("A"),
        node_type: NodeType::Leaf(o!("A")),
        mutations: BTreeMap::from([
          (0, ('G', 'A')), //
          (7, ('T', 'C')), //
        ]),
        gaps: vec![(11, 13)],
        ambiguous: vec![(8, 10)],
        undetermined: vec![(8, 10), (11, 13)],
        mixed: vec![],
        non_consensus: BTreeMap::from([
          (0, btreeset![]), //
          (5, btreeset![]), //
          (7, btreeset![]), //
        ]),
        nuc_composition: btreemap! {
          'A' => 3,
          'C' => 4,
          'G' => 3,
          'T' => 4,
        },
        seq: vec![],
        subtree_profile_variable: Default::default(),
        subtree_profile_fixed: Default::default(),
        profile_variable: Default::default(),
        profile_fixed: Default::default(),
        outgroup_profile_variable: Default::default(),
        outgroup_profile_fixed: Default::default(),
        expQt: Default::default(),
      },
      Node {
        name: o!("B"),
        node_type: NodeType::Leaf(o!("B")),
        mutations: BTreeMap::from([
          (5, ('G', 'C')), //
        ]),
        gaps: vec![(11, 12)],
        ambiguous: vec![(12, 13)],
        undetermined: vec![(11, 13)],
        mixed: vec![],
        non_consensus: BTreeMap::from([
          (0, btreeset![]), //
          (5, btreeset![]), //
          (7, btreeset![]), //
        ]),
        nuc_composition: btreemap! {
          'A' => 2,
          'C' => 4,
          'G' => 3,
          'T' => 5,
        },
        seq: vec![],
        subtree_profile_variable: Default::default(),
        subtree_profile_fixed: Default::default(),
        profile_variable: Default::default(),
        profile_fixed: Default::default(),
        outgroup_profile_variable: Default::default(),
        outgroup_profile_fixed: Default::default(),
        expQt: Default::default(),
      },
      Node {
        name: o!("CD"),
        node_type: NodeType::Internal(o!("CD")),
        mutations: BTreeMap::from([]),
        gaps: vec![],
        ambiguous: vec![],
        undetermined: vec![],
        mixed: vec![],
        non_consensus: BTreeMap::from([]),
        nuc_composition: btreemap! {
          'A' => 1,
          'C' => 2,
          'G' => 6,
          'T' => 5,
        },
        seq: vec![],
        subtree_profile_variable: Default::default(),
        subtree_profile_fixed: Default::default(),
        profile_variable: Default::default(),
        profile_fixed: Default::default(),
        outgroup_profile_variable: Default::default(),
        outgroup_profile_fixed: Default::default(),
        expQt: Default::default(),
      },
      Node {
        name: o!("C"),
        node_type: NodeType::Leaf(o!("C")),
        mutations: BTreeMap::from([
          (0, ('T', 'C')), //
          (6, ('G', 'A')), //
        ]),
        gaps: vec![],
        ambiguous: vec![],
        undetermined: vec![],
        mixed: vec![],
        non_consensus: BTreeMap::from([
          (0, btreeset![]), //
          (5, btreeset![]), //
          (6, btreeset![]), //
        ]),
        nuc_composition: btreemap! {
          'A' => 2,
          'C' => 3,
          'G' => 5,
          'T' => 4,
        },
        seq: vec![],
        subtree_profile_variable: Default::default(),
        subtree_profile_fixed: Default::default(),
        profile_variable: Default::default(),
        profile_fixed: Default::default(),
        outgroup_profile_variable: Default::default(),
        outgroup_profile_fixed: Default::default(),
        expQt: Default::default(),
      },
      Node {
        name: o!("D"),
        node_type: NodeType::Leaf(o!("D")),
        mutations: BTreeMap::from([
          (5, ('G', 'C')), //
        ]),
        gaps: vec![],
        ambiguous: vec![],
        undetermined: vec![],
        mixed: vec![MixedSite::new(10, 'R')],
        non_consensus: BTreeMap::from([
          (0, btreeset![]),     //
          (5, btreeset![]),     //
          (6, btreeset![]),     //
          (10, btreeset!['G']), //
        ]),
        nuc_composition: btreemap! {
          'A' => 1,
          'C' => 3,
          'G' => 5,
          'T' => 5,
        },
        seq: vec![],
        subtree_profile_variable: Default::default(),
        subtree_profile_fixed: Default::default(),
        profile_variable: Default::default(),
        profile_fixed: Default::default(),
        outgroup_profile_variable: Default::default(),
        outgroup_profile_fixed: Default::default(),
        expQt: Default::default(),
      },
    ];

    assert_eq!(&expected, &actual);

    let nuc_counts: BTreeMap<String, usize> = actual
      .iter()
      .map(|node| (node.name.clone(), node.nuc_composition.values().sum()))
      .collect();

    let nuc_counts_expected = btreemap! {
      o!("A") => L,
      o!("AB") => L,
      o!("B") => L,
      o!("C") => L,
      o!("CD") => L,
      o!("D") => L,
      o!("root") => L,
    };

    assert_eq!(nuc_counts_expected, nuc_counts);

    Ok(())
  }

  #[rstest]
  fn test_seq_leaf_reconstruction() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let mut rng = get_random_number_generator(Some(42));

    let inputs = BTreeMap::from([
      (o!("A"), o!("ACATCGCCNNA--G")),
      (o!("B"), o!("GCATCCCTGTA-NG")),
      (o!("C"), o!("CCGGCGATGTATTG")),
      (o!("D"), o!("TCGGCCGTGTRTTG")),
    ]);

    let graph = create_graph_from_nwk_str::<Node, Edge>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    compress_sequences(&inputs, &graph, &mut rng).unwrap();

    let actual = reconstruct_leaf_sequences(&graph)?;

    assert_eq!(
      json_stringify(&inputs, JsonPretty(false))?,
      json_stringify(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  #[rstest]
  fn test_seq_ancestral_reconstruction() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let mut rng = get_random_number_generator(Some(42));

    let inputs = BTreeMap::from([
      (o!("A"), o!("ACATCGCCNNA--G")),
      (o!("B"), o!("GCATCCCTGTA-NG")),
      (o!("C"), o!("CCGGCGATGTATTG")),
      (o!("D"), o!("TCGGCCGTGTRTTG")),
    ]);

    #[rustfmt::skip]
    let expected = BTreeMap::from([
      (o!("AB"),   o!("GCATCGCTGTATTG")),
      (o!("CD"),   o!("TCGGCGGTGTATTG")),
      (o!("root"), o!("TCGGCGGTGTATTG")),
    ]);

    let graph = create_graph_from_nwk_str::<Node, Edge>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    compress_sequences(&inputs, &graph, &mut rng).unwrap();

    let actual = Arc::new(RwLock::new(BTreeMap::new()));
    ancestral_reconstruction_fitch(&graph, false, |node, seq| {
      actual
        .write_arc()
        .insert(node.name.clone(), vec_to_string(seq.to_owned()));
    })?;

    assert_eq!(
      json_stringify(&expected, JsonPretty(false))?,
      json_stringify(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  #[rstest]
  fn test_seq_ancestral_reconstruction_with_leaves() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let mut rng = get_random_number_generator(Some(42));

    let inputs = BTreeMap::from([
      (o!("A"), o!("ACATCGCCNNA--G")),
      (o!("B"), o!("GCATCCCTGTA-NG")),
      (o!("C"), o!("CCGGCGATGTATTG")),
      (o!("D"), o!("TCGGCCGTGTRTTG")),
    ]);

    #[rustfmt::skip]
    let expected = BTreeMap::from([
      (o!("A"),    o!("ACATCGCCNNA--G")),
      (o!("AB"),   o!("GCATCGCTGTATTG")),
      (o!("B"),    o!("GCATCCCTGTA-NG")),
      (o!("C"),    o!("CCGGCGATGTATTG")),
      (o!("CD"),   o!("TCGGCGGTGTATTG")),
      (o!("D"),    o!("TCGGCCGTGTRTTG")),
      (o!("root"), o!("TCGGCGGTGTATTG")),
    ]);

    let graph = create_graph_from_nwk_str::<Node, Edge>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    compress_sequences(&inputs, &graph, &mut rng).unwrap();

    let actual = Arc::new(RwLock::new(BTreeMap::new()));
    ancestral_reconstruction_fitch(&graph, true, |node, seq| {
      actual
        .write_arc()
        .insert(node.name.clone(), vec_to_string(seq.to_owned()));
    })?;

    assert_eq!(
      json_stringify(&expected, JsonPretty(false))?,
      json_stringify(&actual, JsonPretty(false))?
    );

    Ok(())
  }
}
