use crate::model::comment::{EdgeComment, NodeComment};
use crate::model::data::{NewickEdgeData, NewickHybrid, NewickNodeData, SupportSource};
use crate::model::graph::NewickGraph;
use crate::model::validate::validate_graph;
use std::collections::BTreeMap;

pub(crate) fn graphs_equal(left: &NewickGraph, right: &NewickGraph, ordered: bool) -> bool {
  if left.rooted() != right.rooted()
    || left.weight().map(f64::to_bits) != right.weight().map(f64::to_bits)
    || EdgeKey::of(left.root_edge()) != EdgeKey::of(right.root_edge())
  {
    return false;
  }
  match (validate_graph(left), validate_graph(right)) {
    (Ok(()), Ok(())) => {
      let mut interner = Interner::default();
      interner.subgraph_id(left, ordered) == interner.subgraph_id(right, ordered)
    },
    (Err(_), Err(_)) => fields_equal(left, right),
    (Ok(()), Err(_)) | (Err(_), Ok(())) => false,
  }
}

fn fields_equal(left: &NewickGraph, right: &NewickGraph) -> bool {
  left.root() == right.root()
    && left.node_count() == right.node_count()
    && left.edge_count() == right.edge_count()
    && left.nodes().zip(right.nodes()).all(|((idx, a), (_, b))| {
      NodeKey::of(a) == NodeKey::of(b)
        && left.child_edges(idx) == right.child_edges(idx)
        && left.parent_edges(idx) == right.parent_edges(idx)
    })
    && left.edges().zip(right.edges()).all(|((_, a), (_, b))| {
      a.parent() == b.parent() && a.child() == b.child() && EdgeKey::of(a.data()) == EdgeKey::of(b.data())
    })
}

#[derive(Default)]
struct Interner<'g> {
  nodes: BTreeMap<NodeKey<'g>, usize>,
  edges: BTreeMap<EdgeKey<'g>, usize>,
  subgraphs: BTreeMap<(usize, Vec<(usize, usize)>), usize>,
}

impl<'g> Interner<'g> {
  fn subgraph_id(&mut self, graph: &'g NewickGraph, ordered: bool) -> usize {
    let mut ids = vec![0_usize; graph.node_count()];
    for node in graph.postorder() {
      let mut children: Vec<(usize, usize)> = graph
        .child_edges(node)
        .iter()
        .map(|&edge| {
          let entry = graph.edge(edge);
          (intern(&mut self.edges, EdgeKey::of(entry.data())), ids[entry.child()])
        })
        .collect();
      if !ordered {
        children.sort_unstable();
      }
      let node_id = intern(&mut self.nodes, NodeKey::of(graph.node(node)));
      ids[node] = intern(&mut self.subgraphs, (node_id, children));
    }
    ids[graph.root()]
  }
}

fn intern<K: Ord>(table: &mut BTreeMap<K, usize>, key: K) -> usize {
  let next = table.len();
  *table.entry(key).or_insert(next)
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct NodeKey<'g> {
  name: Option<&'g str>,
  hybrid: Option<&'g NewickHybrid>,
  comments: &'g [NodeComment],
}

impl<'g> NodeKey<'g> {
  fn of(node: &'g NewickNodeData) -> Self {
    Self {
      name: node.name(),
      hybrid: node.hybrid(),
      comments: node.comments(),
    }
  }
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct EdgeKey<'g> {
  length: Option<u64>,
  support: Vec<u64>,
  support_source: Option<SupportSource>,
  probability: Option<u64>,
  is_acceptor: bool,
  occurrence_comments: &'g [NodeComment],
  comments: &'g [EdgeComment],
}

impl<'g> EdgeKey<'g> {
  fn of(edge: &'g NewickEdgeData) -> Self {
    Self {
      length: edge.branch_length().map(f64::to_bits),
      support: edge.support().iter().map(|value| value.to_bits()).collect(),
      support_source: (!edge.support().is_empty()).then(|| edge.support_source()),
      probability: edge.probability().map(f64::to_bits),
      is_acceptor: edge.is_acceptor(),
      occurrence_comments: edge.occurrence_comments(),
      comments: edge.comments(),
    }
  }
}
