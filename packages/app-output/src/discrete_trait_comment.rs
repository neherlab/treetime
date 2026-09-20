use eyre::Report;
use maplit::btreemap;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::NodeCommentProvider;

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
