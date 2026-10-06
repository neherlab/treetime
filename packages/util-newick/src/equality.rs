use crate::types::{NewickEdgeData, NewickGraph, NewickHybrid, NewickLabel, NewickNodeData, NewickValue};
use crate::validate::{postorder, validate_graph};
use std::collections::BTreeMap;

pub(crate) fn graphs_equal(left: &NewickGraph, right: &NewickGraph, ordered: bool) -> bool {
  if left.rooted != right.rooted {
    return false;
  }
  match (validate_graph(left), validate_graph(right)) {
    (Ok(()), Ok(())) => {
      let mut interner = Interner::default();
      let left_id = interner.subgraph_id(left, ordered);
      let right_id = interner.subgraph_id(right, ordered);
      left_id == right_id
    },
    (Err(_), Err(_)) => fields_equal(left, right),
    (Ok(()), Err(_)) | (Err(_), Ok(())) => false,
  }
}

fn fields_equal(left: &NewickGraph, right: &NewickGraph) -> bool {
  left.root == right.root
    && left.nodes.len() == right.nodes.len()
    && left.edges.len() == right.edges.len()
    && left
      .nodes
      .iter()
      .zip(&right.nodes)
      .all(|(a, b)| NodeKey::of(a) == NodeKey::of(b) && a.children == b.children)
    && left
      .edges
      .iter()
      .zip(&right.edges)
      .all(|(a, b)| a.parent == b.parent && a.child == b.child && EdgeKey::of(&a.data) == EdgeKey::of(&b.data))
}

#[derive(Default)]
struct Interner<'g> {
  nodes: BTreeMap<NodeKey<'g>, usize>,
  edges: BTreeMap<EdgeKey<'g>, usize>,
  subgraphs: BTreeMap<(usize, Vec<(usize, usize)>), usize>,
}

impl<'g> Interner<'g> {
  fn subgraph_id(&mut self, graph: &'g NewickGraph, ordered: bool) -> Option<usize> {
    let order = postorder(graph)?;
    let mut ids = vec![0_usize; graph.nodes.len()];
    for node_idx in order {
      let node = &graph.nodes[node_idx];
      let mut children: Vec<(usize, usize)> = node
        .children
        .iter()
        .map(|&edge_idx| {
          let edge = &graph.edges[edge_idx];
          (intern(&mut self.edges, EdgeKey::of(&edge.data)), ids[edge.child])
        })
        .collect();
      if !ordered {
        children.sort_unstable();
      }
      let node_id = intern(&mut self.nodes, NodeKey::of(node));
      ids[node_idx] = intern(&mut self.subgraphs, (node_id, children));
    }
    Some(ids[graph.root])
  }
}

fn intern<K: Ord>(table: &mut BTreeMap<K, usize>, key: K) -> usize {
  let next = table.len();
  *table.entry(key).or_insert(next)
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct NodeKey<'g> {
  label: Option<&'g NewickLabel>,
  attrs: &'g BTreeMap<String, NewickValue>,
  raw_comments: &'g [String],
  hybrid: Option<&'g NewickHybrid>,
}

impl<'g> NodeKey<'g> {
  fn of(node: &'g NewickNodeData) -> Self {
    Self {
      label: node.label.as_ref(),
      attrs: &node.node_attrs,
      raw_comments: &node.raw_comments,
      hybrid: node.hybrid.as_ref(),
    }
  }
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct EdgeKey<'g> {
  length: Option<u64>,
  attrs: &'g BTreeMap<String, NewickValue>,
  raw_comments: &'g [String],
  is_acceptor: bool,
}

impl<'g> EdgeKey<'g> {
  fn of(edge: &'g NewickEdgeData) -> Self {
    Self {
      length: edge.branch_length.map(f64::to_bits),
      attrs: &edge.branch_attrs,
      raw_comments: &edge.raw_comments,
      is_acceptor: edge.is_acceptor,
    }
  }
}
