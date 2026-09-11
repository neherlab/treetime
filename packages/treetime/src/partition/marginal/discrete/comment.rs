use crate::partition::marginal::discrete::partition::PartitionMarginalDiscrete;
use eyre::Report;
use maplit::btreemap;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::NodeCommentProvider;

pub struct DiscreteCommentProvider<'a> {
  partition: &'a PartitionMarginalDiscrete,
  attribute: &'a str,
}

impl<'a> DiscreteCommentProvider<'a> {
  pub fn new(partition: &'a PartitionMarginalDiscrete, attribute: &'a str) -> Self {
    Self { partition, attribute }
  }
}

impl NodeCommentProvider for DiscreteCommentProvider<'_> {
  fn node_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report> {
    Ok(
      self
        .partition
        .get_reconstructed_trait(key)
        .map_or_else(BTreeMap::new, |trait_value| {
          btreemap! {
            self.attribute.to_owned() => trait_value,
          }
        }),
    )
  }
}

/// Newick/Nexus node-comment provider that reads a gathered per-node reconstructed-trait map.
///
/// Mirrors [`DiscreteCommentProvider`], but reads the trait from a value map instead of the partition,
/// so the tree writers no longer touch the partition during serialization.
pub struct DiscreteTraitCommentProvider<'a> {
  reconstructed_traits: &'a BTreeMap<GraphNodeKey, Option<String>>,
  attribute: &'a str,
}

impl<'a> DiscreteTraitCommentProvider<'a> {
  pub fn new(reconstructed_traits: &'a BTreeMap<GraphNodeKey, Option<String>>, attribute: &'a str) -> Self {
    Self {
      reconstructed_traits,
      attribute,
    }
  }
}

impl NodeCommentProvider for DiscreteTraitCommentProvider<'_> {
  fn node_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report> {
    Ok(
      self.reconstructed_traits[&key]
        .clone()
        .map_or_else(BTreeMap::new, |trait_value| {
          btreemap! {
            self.attribute.to_owned() => trait_value,
          }
        }),
    )
  }
}
