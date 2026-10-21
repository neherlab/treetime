use crate::graph::Graph;
use crate::node::GraphNodeKey;
use crate::schedule::GraphSchedule;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use std::ops::Index;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NodeNames {
  by_node: BTreeMap<GraphNodeKey, String>,
  in_use: BTreeSet<String>,
  next_number: usize,
}

impl NodeNames {
  pub fn from_input(
    mut input: BTreeMap<GraphNodeKey, Option<String>>,
    graph: &Graph,
  ) -> Result<InputNodeNames, Report> {
    let mut names = Self {
      by_node: BTreeMap::new(),
      in_use: BTreeSet::new(),
      next_number: 0,
    };
    let mut duplicate_names: BTreeSet<String> = BTreeSet::new();
    for node in graph.get_nodes() {
      let key = node.key();
      let Some(name) = input.remove(&key).flatten().filter(|name| !name.is_empty()) else {
        continue;
      };
      if !names.in_use.insert(name.clone()) {
        duplicate_names.insert(name.clone());
      }
      names.by_node.insert(key, name);
    }

    let unnamed_nodes: BTreeSet<GraphNodeKey> = graph
      .get_nodes()
      .map(|node| node.key())
      .filter(|key| !names.by_node.contains_key(key))
      .collect();
    if !unnamed_nodes.is_empty() {
      for &key in GraphSchedule::new(graph)?.preorder() {
        if unnamed_nodes.contains(&key) {
          names.name_new_node(key);
        }
      }
    }

    Ok(InputNodeNames {
      names,
      unnamed_nodes,
      duplicate_names: duplicate_names.into_iter().collect(),
    })
  }

  pub fn name_new_node(&mut self, key: GraphNodeKey) {
    let mut name = generated_name(self.next_number);
    while self.in_use.contains(&name) {
      self.next_number += 1;
      name = generated_name(self.next_number);
    }
    self.next_number += 1;
    self.in_use.insert(name.clone());
    self.by_node.insert(key, name);
  }
}

impl Index<GraphNodeKey> for NodeNames {
  type Output = str;

  fn index(&self, key: GraphNodeKey) -> &str {
    &self.by_node[&key]
  }
}

pub struct InputNodeNames {
  pub names: NodeNames,
  pub unnamed_nodes: BTreeSet<GraphNodeKey>,
  pub duplicate_names: Vec<String>,
}

fn generated_name(number: usize) -> String {
  format!("NODE_{number:07}")
}
