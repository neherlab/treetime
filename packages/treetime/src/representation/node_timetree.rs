use crate::distribution::distribution::Distribution;
use crate::graph::node::{GraphNode, Named};
use crate::io::nwk::NodeFromNwk;
use crate::representation::graph_ancestral::NodeAncestral;
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct NodeTimetree {
  pub name: Option<String>,
  pub desc: Option<String>,
  pub time: Option<f64>,
  pub time_before_present: Option<f64>,
  pub time_distribution: Option<Arc<Distribution>>,
  pub bad_branch: bool,
}

impl GraphNode for NodeTimetree {}

impl Named for NodeTimetree {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

impl NodeFromNwk for NodeTimetree {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      ..NodeTimetree::default()
    })
  }
}

impl From<&NodeAncestral> for NodeTimetree {
  fn from(node: &NodeAncestral) -> Self {
    Self {
      name: node.name.clone(),
      desc: node.desc.clone(),
      time: None,
      time_before_present: None,
      time_distribution: None,
      bad_branch: node.bad_branch,
    }
  }
}

impl crate::io::nwk::NodeToNwk for NodeTimetree {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    BTreeMap::new()
  }
}
