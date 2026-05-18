use crate::partition::discrete::PartitionDiscrete;
use maplit::btreemap;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::NodeCommentProvider;

/// Provide reconstructed discrete traits as node comments during tree serialization.
pub struct PartitionCommentProvider<'a> {
  partition: &'a PartitionDiscrete,
  attribute: &'a str,
}

impl<'a> PartitionCommentProvider<'a> {
  pub fn new(partition: &'a PartitionDiscrete, attribute: &'a str) -> Self {
    Self { partition, attribute }
  }
}

impl NodeCommentProvider for PartitionCommentProvider<'_> {
  fn node_comments(&self, key: GraphNodeKey) -> BTreeMap<String, String> {
    self
      .partition
      .get_reconstructed_trait(key)
      .map_or_else(BTreeMap::new, |trait_value| {
        btreemap! {
          self.attribute.to_owned() => trait_value,
        }
      })
  }
}
