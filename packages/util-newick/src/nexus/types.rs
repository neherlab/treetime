use crate::model::comment::NewickComment;
use crate::model::graph::NewickGraph;
use crate::read::error::NewickWarning;
use crate::read::options::NewickTree;
use crate::write::options::NewickWriteOptions;
use deser::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct NexusFile {
  pub trees: Vec<NexusTree>,
  pub skipped: Vec<NexusCommand>,
  pub warnings: Vec<NewickWarning>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct NexusTree {
  pub name: String,
  pub tree: NewickTree,
  pub comments: Vec<NewickComment>,
}

impl NexusTree {
  pub fn as_ref(&self) -> NexusTreeRef<'_> {
    NexusTreeRef {
      name: &self.name,
      graph: &self.tree.graph,
      comments: &self.comments,
    }
  }
}

#[derive(Clone, Copy, Debug)]
pub struct NexusTreeRef<'t> {
  pub name: &'t str,
  pub graph: &'t NewickGraph,
  pub comments: &'t [NewickComment],
}

impl<'t> NexusTreeRef<'t> {
  pub fn new(name: &'t str, graph: &'t NewickGraph) -> Self {
    Self {
      name,
      graph,
      comments: &[],
    }
  }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct NexusCommand {
  pub block: Option<String>,
  pub command: String,
  pub line: usize,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct NexusWriteOptions {
  pub newick: NewickWriteOptions,
  pub translate: bool,
}
