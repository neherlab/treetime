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
