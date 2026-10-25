use crate::graph::Graph;
use crate::key_map::NodeMap;
use crate::node::GraphNodeKey;
use crate::schedule::GraphSchedule;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use std::ops::Index;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NodeNames {
  by_node: NodeMap<String>,
  repeated_input_names: BTreeMap<GraphNodeKey, String>,
  in_use: BTreeSet<String>,
  next_number: usize,
}

impl NodeNames {
  pub fn from_input(mut input: NodeMap<Option<String>>, graph: &Graph) -> Result<InputNodeNames, Report> {
    let named: Vec<(GraphNodeKey, String)> = graph
      .get_nodes()
      .filter_map(|node| {
        let key = node.key();
        let name = input.remove(key).flatten().filter(|name| !name.is_empty())?;
        Some((key, name))
      })
      .collect();

    let mut names = Self {
      by_node: NodeMap::new(),
      repeated_input_names: BTreeMap::new(),
      in_use: named.iter().map(|(_, name)| name.clone()).collect(),
      next_number: 0,
    };
    let mut next_occurrence: BTreeMap<&str, usize> = BTreeMap::new();
    for (key, name) in &named {
      let Some(occurrence) = next_occurrence.get_mut(name.as_str()) else {
        next_occurrence.insert(name, 2);
        names.by_node.insert(*key, name.clone());
        continue;
      };
      let unique = names.unique_repeat_name(name, occurrence);
      names.by_node.insert(*key, unique);
      names.repeated_input_names.insert(*key, name.clone());
    }
    let duplicate_names = next_occurrence
      .into_iter()
      .filter(|(_, occurrence)| *occurrence > 2)
      .map(|(name, _)| name.to_owned())
      .collect();

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
      duplicate_names,
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

  pub fn pairing_name(&self, key: GraphNodeKey) -> &str {
    self
      .repeated_input_names
      .get(&key)
      .unwrap_or_else(|| &self.by_node[&key])
  }

  fn unique_repeat_name(&mut self, name: &str, occurrence: &mut usize) -> String {
    let mut unique = format!("{name}_{occurrence}");
    while self.in_use.contains(&unique) {
      *occurrence += 1;
      unique = format!("{name}_{occurrence}");
    }
    *occurrence += 1;
    self.in_use.insert(unique.clone());
    unique
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
