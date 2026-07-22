use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{Described, GraphNode, Named};
use treetime_io::graphviz::{EdgeToGraphviz, NodeToGraphviz};
use treetime_io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk};

pub type GraphAncestral<D = ()> = Graph<NodeAncestral, EdgeAncestral, D>;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct NodeAncestral {
  pub name: Option<String>,
  pub desc: Option<String>,
  /// Input-tree branch support / bootstrap / posterior probability.
  pub confidence: Option<f64>,
}

impl NodeFromNwk for NodeAncestral {
  fn from_nwk(
    name: Option<impl AsRef<str>>,
    confidence: Option<f64>,
    _: &BTreeMap<String, String>,
  ) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      confidence,
      ..NodeAncestral::default()
    })
  }
}

impl NodeToNwk for NodeAncestral {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    BTreeMap::new()
  }
}

impl GraphNode for NodeAncestral {}

impl Named for NodeAncestral {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

impl Described for NodeAncestral {
  fn desc(&self) -> &Option<String> {
    &self.desc
  }

  fn set_desc(&mut self, desc: Option<String>) {
    self.desc = desc;
  }
}

impl NodeToGraphviz for NodeAncestral {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct EdgeAncestral {
  pub branch_length: Option<f64>,
}

impl GraphEdge for EdgeAncestral {}

impl HasBranchLength for EdgeAncestral {
  fn branch_length(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_branch_length(&mut self, weight: Option<f64>) {
    self.branch_length = weight;
  }
}

impl EdgeFromNwk for EdgeAncestral {
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self { branch_length })
  }
}

impl EdgeToNwk for EdgeAncestral {
  fn nwk_weight(&self) -> Option<f64> {
    self.branch_length()
  }
}

impl EdgeToGraphviz for EdgeAncestral {}

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::partition::sparse::{SparseEdgePartition, SparseNodePartition};
  use crate::partition::traits::MutationCommentProvider;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::indel::InDel;
  use crate::seq::mutation::Sub;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::{NodeCommentProvider, nwk_read_str};
  use treetime_primitives::{AsciiChar, Seq};

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn make_test_partition(
    graph: &GraphAncestral,
    length: usize,
    edge_subs: &[(usize, Vec<Sub>)],
  ) -> Result<PartitionMarginalSparse, Report> {
    let alphabet = Alphabet::default();
    let mut ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
    for (_, subs) in edge_subs {
      for s in subs {
        if s.pos() < length {
          ref_seq[s.pos()] = s.reff();
        }
      }
    }

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: alphabet.clone(),
      length,
      root_sequence: ref_seq.clone(),
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      let mut node_part = SparseNodePartition::empty(&alphabet);
      node_part.seq.sequence = ref_seq.clone();
      partition.nodes.insert(key, node_part);
    }

    let edges = graph.get_edges();
    for (idx, subs) in edge_subs {
      if let Some(edge) = edges.get(*idx) {
        let edge_key = edge.read_arc().key();
        let mut edge_part = SparseEdgePartition::with_fitch_subs(subs.clone());
        edge_part.set_ml_subs(subs.clone());
        partition.edges.insert(edge_key, edge_part);
      }
    }

    Ok(partition)
  }

  fn leaf_key(graph: &GraphAncestral) -> GraphNodeKey {
    graph.get_leaves()[0].read_arc().key()
  }

  #[test]
  fn test_mutation_comment_provider_formats_1_based_substitutions_and_indels() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let mut partition = make_test_partition(
      &graph,
      100,
      &[(
        0,
        vec![
          Sub::new(c(b'A'), 0_usize, c(b'T'))?,
          Sub::new(c(b'G'), 5_usize, c(b'C'))?,
        ],
      )],
    )?;
    let edge_key = graph.get_edges()[0].read_arc().key();
    partition
      .edges
      .get_mut(&edge_key)
      .expect("fixture edge partition must exist")
      .indels = vec![InDel::del((1, 3), Seq::try_from_str("CG")?)?];
    let provider = MutationCommentProvider::new(&partition, &graph);
    let comments = provider.node_comments(leaf_key(&graph))?;
    assert_eq!(comments.get("mutations").map(String::as_str), Some("A1T,C2-,G3-,G6C"));
    Ok(())
  }

  #[test]
  fn test_mutation_comment_provider_root_has_no_comments() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let partition = make_test_partition(&graph, 100, &[(0, vec![Sub::new(c(b'A'), 0_usize, c(b'T'))?])])?;
    let provider = MutationCommentProvider::new(&partition, &graph);
    let root_key = graph.get_roots()[0].read_arc().key();
    let comments = provider.node_comments(root_key)?;
    assert!(comments.is_empty());
    Ok(())
  }

  #[test]
  fn test_mutation_comment_provider_no_mutations_returns_empty() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let partition = make_test_partition(&graph, 100, &[(0, vec![])])?;
    let provider = MutationCommentProvider::new(&partition, &graph);
    let comments = provider.node_comments(leaf_key(&graph))?;
    assert!(comments.is_empty());
    Ok(())
  }

  #[test]
  fn test_mutation_comment_provider_sorts_by_position() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let partition = make_test_partition(
      &graph,
      100,
      &[(
        0,
        vec![
          Sub::new(c(b'C'), 50_usize, c(b'G'))?,
          Sub::new(c(b'A'), 10_usize, c(b'T'))?,
          Sub::new(c(b'G'), 30_usize, c(b'C'))?,
        ],
      )],
    )?;
    let provider = MutationCommentProvider::new(&partition, &graph);
    let comments = provider.node_comments(leaf_key(&graph))?;
    assert_eq!(comments.get("mutations").map(String::as_str), Some("A11T,G31C,C51G"));
    Ok(())
  }
}
