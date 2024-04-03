use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::{Graph, GraphNodeBackward};
use crate::graph::node::{GraphNode, GraphNodeKey, Named, NodeType, WithNwkComments};
use crate::seq::find_char_ranges::{find_ambiguous_ranges, find_gap_ranges};
use crate::seq::find_mixed_sites::{find_mixed_sites, MixedSite};
use crate::seq::range::range_contains;
use crate::seq::range_intersection::range_intersection_iter;
use crate::seq::range_union::range_union;
use crate::seq::sets::{sets_intersection, sets_union};
use crate::utils::random::random_pop;
use crate::utils::string::vec_to_string;
use crate::{make_error, make_internal_report};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use maplit::{btreemap, btreeset};
use rand::Rng;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::{Display, Formatter};
use std::ops::DerefMut;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Node {
  pub name: String,
  pub node_type: NodeType,

  pub mutations: BTreeMap<usize, (char, char)>,
  pub gaps: Vec<(usize, usize)>,
  pub ambiguous: Vec<(usize, usize)>,
  pub undetermined: Vec<(usize, usize)>,
  pub mixed: Vec<MixedSite>,
  pub non_consensus: BTreeMap<usize, BTreeSet<char>>,

  pub seq: Vec<char>,
}

impl Node {
  pub fn new(name: impl AsRef<str>, node_type: NodeType) -> Self {
    Self {
      name: name.as_ref().to_owned(),
      node_type,

      mutations: BTreeMap::from([]),
      gaps: vec![],
      ambiguous: vec![],
      undetermined: vec![],
      mixed: vec![],
      non_consensus: BTreeMap::new(),

      seq: vec![],
    }
  }

  /// Gather all states at a given position from all child nodes
  pub fn get_letter_disambiguated(&self, pos: usize) -> BTreeSet<char> {
    self
      // Get possible states if non-consensus
      .non_consensus.get(&pos).cloned()
      // Otherwise read the sequence at that position
      .unwrap_or_else(|| btreeset!{self.seq[pos]})
  }

  #[inline]
  pub const fn is_root(&self) -> bool {
    matches!(self.node_type, NodeType::Root(_))
  }

  #[inline]
  pub const fn is_internal(&self) -> bool {
    matches!(self.node_type, NodeType::Internal(_))
  }

  #[inline]
  pub const fn is_leaf(&self) -> bool {
    matches!(self.node_type, NodeType::Leaf(_))
  }
}

impl GraphNode for Node {
  fn root(name: &str) -> Self {
    Self::new(name, NodeType::Root(name.to_owned()))
  }

  fn internal(name: &str) -> Self {
    Self::new(name, NodeType::Internal(name.to_owned()))
  }

  fn leaf(name: &str) -> Self {
    Self::new(name, NodeType::Leaf(name.to_owned()))
  }

  fn set_node_type(&mut self, node_type: NodeType) {
    self.node_type = node_type;
  }
}

impl WithNwkComments for Node {}

impl Named for Node {
  fn name(&self) -> &str {
    &self.name
  }

  fn set_name(&mut self, name: &str) {
    self.name = name.to_owned();
  }
}

impl Display for Node {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    match &self.node_type {
      NodeType::Root(weight) => write!(f, "{weight:1.4}"),
      NodeType::Internal(weight) => write!(f, "{weight:1.4}"),
      NodeType::Leaf(name) => write!(f, "{name}"),
    }
  }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Edge {
  pub weight: f64,
}

impl GraphEdge for Edge {
  fn new(weight: f64) -> Self {
    Self { weight }
  }
}

impl Weighted for Edge {
  fn weight(&self) -> f64 {
    self.weight
  }
}

impl Display for Edge {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{:1.4}", &self.weight)
  }
}

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
  for pos in 0..L {
    // Skip ambiguous and gaps
    if range_contains(&n.undetermined, pos) {
      continue;
    }

    if non_consensus_positions.contains(&pos) {
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

          // Adjust children, now that we know they can disagree
          // TODO: avoid code duplication
          children.iter_mut().for_each(|c| {
            // TODO: double check that. Do we need to perform union here if position is already taken?
            c.non_consensus.entry(pos).or_insert_with(|| btreeset![c.seq[pos]]);
          });
        }
        states => {
          // More than one element in the intersection of child states. Propagate intersection.
          let states = states.iter().copied().collect();
          n.non_consensus.insert(pos, states);

          // Adjust children, now that we know they can disagree
          // TODO: avoid code duplication
          children.iter_mut().for_each(|c| {
            // TODO: double check that. Do we need to perform union here if position is already taken?
            c.non_consensus.entry(pos).or_insert_with(|| btreeset![c.seq[pos]]);
          });
        }
      };
    } else {
      let child_states = gather_consensus_child_states(children, pos);
      match child_states.as_slice() {
        [state] => {
          // Same character in all child sequences
          n.seq[pos] = *state;
        }
        [] => {
          // We are not in non-consensus position, but children have no consensus characters.
          // This is impossible, and the fact that we have to handle that might hint at bad code design.
          unreachable!("Child sequences contain no characters");
        }
        states => {
          // Child states are different, propagate union
          let states = states.iter().copied().collect();
          n.non_consensus.insert(pos, states);

          // Adjust children, now that we know they can disagree
          // TODO: avoid code duplication
          children.iter_mut().for_each(|c| {
            c.non_consensus.entry(pos).or_insert_with(|| btreeset! {c.seq[pos]});
          });
        }
      }
    }
  }
}

pub fn gather_child_state_sets(children: &[&mut Node], pos: usize) -> Vec<BTreeSet<char>> {
  children
    .iter()
    .filter(|c| !range_contains(&c.undetermined, pos))
    .map(|c| c.get_letter_disambiguated(pos).into_iter().collect::<BTreeSet<_>>())
    .unique()
    .collect_vec()
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
        child.mutations.insert(pos, (parent_state, child_state));
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
  let path = graph
    .path_from_leaf_to_root(node_key)?
    .ok_or_else(|| make_internal_report!("No path from root to node '{node_key}'"))?;

  for anc in path {
    let anc = anc.read_arc().payload().read_arc();
    // Apply mutations
    for (&pos, &(ref_char, _)) in &anc.mutations {
      seq[pos] = ref_char;
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

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::create_graph_from_nwk::create_graph_from_nwk_str;
  use crate::io::json::{json_stringify, JsonPretty};
  use crate::o;
  use crate::utils::random::get_random_number_generator;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

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
          (3, btreeset!['G']),
          (5, btreeset!['C']),
          (6, btreeset!['A', 'C']),
        ]),
        seq: vec!['T', 'C', 'G', 'T', 'C', 'G', 'G', 'T', 'G', 'T', 'A', 'T', 'T', 'G'],
      },
      Node {
        name: o!("AB"),
        node_type: NodeType::Internal(o!("AB")),
        mutations: BTreeMap::from([
          (0, ('T', 'G')), //
          (2, ('G', 'A')), //
          (6, ('G', 'C')), //
        ]),
        gaps: vec![],
        ambiguous: vec![],
        undetermined: vec![(11, 13)],
        mixed: vec![],
        non_consensus: BTreeMap::from([
          (0, btreeset!['A', 'G']),
          (5, btreeset!['G', 'C']),
          (7, btreeset!['C', 'T']),
          (2, btreeset!['A']),
          (3, btreeset!['T']),
          (6, btreeset!['C']),
        ]),
        seq: vec![],
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
          (0, btreeset!['A']), //
          (5, btreeset!['G']), //
          (7, btreeset!['C']), //
        ]),
        seq: vec![],
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
          (0, btreeset!['G']), //
          (5, btreeset!['C']), //
          (7, btreeset!['T']), //
        ]),
        seq: vec![],
      },
      Node {
        name: o!("CD"),
        node_type: NodeType::Internal(o!("CD")),
        mutations: BTreeMap::from([
          (3, ('T', 'G')), //
        ]),
        gaps: vec![],
        ambiguous: vec![],
        undetermined: vec![],
        mixed: vec![],
        non_consensus: BTreeMap::from([
          (0, btreeset!['C', 'T']),
          (5, btreeset!['G', 'C']),
          (6, btreeset!['A', 'G']),
          (2, btreeset!['G']),
          (3, btreeset!['G']),
        ]),
        seq: vec![],
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
          (0, btreeset!['C']), //
          (5, btreeset!['G']), //
          (6, btreeset!['A']), //
        ]),
        seq: vec![],
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
          (10, btreeset!['A', 'G']), //
          (0, btreeset!['T']),       //
          (5, btreeset!['C']),       //
          (6, btreeset!['G']),       //
        ]),
        seq: vec![],
      },
    ];

    assert_eq!(expected, actual);

    Ok(())
  }

  #[rstest]
  fn test_seq_reconstruction() -> Result<(), Report> {
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
}
