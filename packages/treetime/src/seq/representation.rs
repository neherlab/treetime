use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::{Graph, GraphNodeBackward};
use crate::graph::node::{GraphNode, Named, NodeType, WithNwkComments};
use crate::make_error;
use crate::seq::find_char_ranges::find_letter_ranges;
use crate::seq::find_mixed_sites::{find_mixed_sites, MixedSite};
use crate::seq::range_intersection::{range_intersection_iter, ranges_contain};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::{Display, Formatter};
use std::ops::Deref;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Node {
  pub name: String,
  pub node_type: NodeType,

  pub mutations: BTreeMap<usize, (char, char)>,
  pub gaps: Vec<(usize, usize)>,
  pub ambiguous: Vec<(usize, usize)>,
  pub mixed: Vec<MixedSite>,
  pub non_consensus: BTreeMap<usize, Vec<char>>,

  pub seq: Vec<char>,
}

impl Node {
  pub fn new(name: impl AsRef<str>, node_type: NodeType) -> Self {
    Self {
      name: name.as_ref().to_owned(),
      node_type,

      mutations: BTreeMap::new(),
      gaps: vec![],
      ambiguous: vec![],
      mixed: vec![],
      non_consensus: BTreeMap::new(),

      seq: vec![],
    }
  }

  /// Gather all states at a given position from all child nodes
  pub fn get_letter_disambiguated(&self, pos: usize) -> Vec<char> {
    self
      // Get possible states if non-consensus
      .non_consensus.get(&pos).cloned()
      // Otherwise read the sequence at that position
      .unwrap_or_else(|| vec![self.seq[pos]])
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
  fn root(name: &str, weight: f64) -> Self {
    Self::new(name, NodeType::Root(weight))
  }

  fn internal(name: &str, weight: f64) -> Self {
    Self::new(name, NodeType::Internal(weight))
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

/// Gather all states at a given position from all child nodes
#[allow(single_use_lifetimes)]
pub fn gather_child_states<'a>(pos: usize, children: impl Iterator<Item = &'a Node>) -> Vec<char> {
  children
    .flat_map(|c| c.get_letter_disambiguated(pos))
    .unique()
    .collect_vec()
}

#[allow(clippy::needless_range_loop)]
pub fn compress_sequences(seqs: &BTreeMap<String, String>, graph: &mut Graph<Node, Edge>) -> Result<(), Report> {
  let L = get_common_length(seqs.iter())?;

  graph.par_iter_breadth_first_backward(
    |GraphNodeBackward {
       is_root,
       is_leaf,
       key,
       payload: n,
       children,
     }| {
      if is_leaf {
        // At each terminal node, temporarily store the sequence and ranges of N, - and mixed sites
        n.seq = seqs[&n.name].chars().collect();
        n.ambiguous = find_letter_ranges(&n.seq, 'N');
        n.gaps = find_letter_ranges(&n.seq, '-');

        // n.mixed stores the exact character at each mixed positions, the non_consensus stores the possible states
        (n.mixed, n.non_consensus) = find_mixed_sites(&n.seq);
      } else {
        // Positions that are N or - in all children are still N or - in the parent
        let mut children = children.iter().map(|(node, _)| node.write()).collect_vec();
        n.ambiguous = range_intersection_iter(children.iter().map(|c| &c.ambiguous));
        n.gaps = range_intersection_iter(children.iter().map(|c| &c.gaps));

        // All sites that are not N or - but not fixed will need special treatment
        let non_consensus_positions: BTreeSet<usize> =
          children.iter().flat_map(|c| c.non_consensus.keys().copied()).collect();

        n.non_consensus = BTreeMap::new();
        let mut seq = vec![' '; L];
        for pos in 0..L {
          // Skip ambiguous and gaps
          if ranges_contain(&n.ambiguous, pos) || ranges_contain(&n.gaps, pos) {
            continue;
          }

          let states = gather_child_states(pos, children.iter().map(Deref::deref));
          match states[..] {
            [] => {}
            [state] => {
              seq[pos] = state;
            }
            _ => {
              n.non_consensus.insert(pos, states);
            }
          };
        }

        // Deallocate full sequences from children
        children.iter_mut().for_each(|c| c.seq = vec![]);

        n.seq = seq;
      }

      GraphTraversalContinuation::Continue
    },
  );

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::create_graph_from_nwk::create_graph_from_nwk_str;
  use crate::o;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

  #[rstest]
  fn test_seq_representation() -> Result<(), Report> {
    let seqs = BTreeMap::from([
      (o!("A"), o!("ACATCGCCNNA--G")),
      (o!("B"), o!("GCATCCCTGTA-TG")),
      (o!("C"), o!("CCGGCGATGTATTG")),
      (o!("D"), o!("TCGGCCGTGTRTTG")),
    ]);

    let mut graph = create_graph_from_nwk_str::<Node, Edge>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    compress_sequences(&seqs, &mut graph).unwrap();

    let mut actual: Vec<Node> = vec![];
    graph.get_nodes().iter().for_each(|node| {
      actual.push(node.read().payload().read().clone());
    });

    let expected = vec![
      Node {
        name: o!("root"),
        node_type: NodeType::Root(0.0),
        mutations: BTreeMap::new(),
        gaps: vec![],
        ambiguous: vec![],
        mixed: vec![],
        non_consensus: BTreeMap::from([
          (0, vec!['G', 'C', 'A', 'T']),
          (2, vec!['G', 'A']),
          (3, vec!['G', 'T']),
          (5, vec!['G', 'C']),
          (6, vec!['G', 'C', 'A']),
          (11, vec![' ', 'T']), // FIXME: empty char
        ]),
        seq: vec![],
      },
      Node {
        name: o!("AB"),
        node_type: NodeType::Leaf(o!("AB")),
        mutations: BTreeMap::new(),
        gaps: vec![(11, 12)],
        ambiguous: vec![],
        mixed: vec![],
        non_consensus: BTreeMap::from([
          (0, vec!['G', 'A']),
          (5, vec!['G', 'C']),
          (7, vec!['C', 'T']),
          (8, vec!['G', 'N']),
          (9, vec!['N', 'T']),
          (12, vec!['-', 'T']),
          (2, vec!['A']),
          (3, vec!['T']),
          (6, vec!['C']),
          (11, vec![]),
        ]),
        seq: vec![],
      },
      Node {
        name: o!("A"),
        node_type: NodeType::Leaf(o!("A")),
        mutations: BTreeMap::new(),
        gaps: vec![(11, 13)],
        ambiguous: vec![(8, 10)],
        mixed: vec![],
        non_consensus: BTreeMap::from([
          (0, vec!['A']),
          (5, vec!['G']),
          (7, vec!['C']),
          (8, vec!['N']),
          (9, vec!['N']),
          (12, vec!['-']),
        ]),
        seq: vec![],
      },
      Node {
        name: o!("B"),
        node_type: NodeType::Leaf(o!("B")),
        mutations: BTreeMap::new(),
        gaps: vec![(11, 12)],
        ambiguous: vec![],
        mixed: vec![],
        non_consensus: BTreeMap::from([
          (0, vec!['G']),
          (5, vec!['C']),
          (7, vec!['T']),
          (8, vec!['G']),
          (9, vec!['T']),
          (12, vec!['T']),
        ]),
        seq: vec![],
      },
      Node {
        name: o!("CD"),
        node_type: NodeType::Leaf(o!("CD")),
        mutations: BTreeMap::new(),
        gaps: vec![],
        ambiguous: vec![],
        mixed: vec![],
        non_consensus: BTreeMap::from([
          (0, vec!['C', 'T']),
          (5, vec!['G', 'C']),
          (6, vec!['G', 'A']),
          (2, vec!['G']),
          (3, vec!['G']),
          (11, vec!['T']),
        ]),
        seq: vec![],
      },
      Node {
        name: o!("C"),
        node_type: NodeType::Leaf(o!("C")),
        mutations: BTreeMap::new(),
        gaps: vec![],
        ambiguous: vec![],
        mixed: vec![],
        non_consensus: BTreeMap::from([(0, vec!['C']), (5, vec!['G']), (6, vec!['A'])]),
        seq: vec![],
      },
      Node {
        name: o!("D"),
        node_type: NodeType::Leaf(o!("D")),
        mutations: BTreeMap::new(),
        gaps: vec![],
        ambiguous: vec![],
        mixed: vec![MixedSite::new(10, 'R')],
        non_consensus: BTreeMap::from([(10, vec!['G', 'A']), (0, vec!['T']), (5, vec!['C']), (6, vec!['G'])]),
        seq: vec![],
      },
    ];

    assert_eq!(expected, actual);

    Ok(())
  }
}
