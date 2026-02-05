use crate::graph::edge::{GraphEdge, HasBranchLength};
use crate::graph::graph::Graph;
use crate::graph::node::{Described, GraphNode, Named};
use crate::io::graphviz::{EdgeToGraphviz, NodeToGraphviz};
use crate::io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions, format_weight};
use crate::o;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub type GraphAncestral = Graph<NodeAncestral, EdgeAncestral, ()>;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct NodeAncestral {
  pub name: Option<String>,
  pub desc: Option<String>,
}

impl NodeFromNwk for NodeAncestral {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      ..NodeAncestral::default()
    })
  }
}

impl NodeToNwk for NodeAncestral {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mut comments = BTreeMap::new();
    let mutations: String = "".to_owned(); // TODO: fill mutations
    comments.insert(o!("mutations"), mutations);
    comments
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

impl EdgeToGraphviz for EdgeAncestral {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .branch_length()
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.branch_length()
  }
}
