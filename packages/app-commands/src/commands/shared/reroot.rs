use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use treetime::clock::find_best_root::params::{RerootMethod, RerootSpec};
use treetime::make_report;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pair_by_name::first_key_by_name;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct RerootArgs {
  /// Reroot the tree by temporal-signal optimization.
  ///
  /// Defaults to least-squares when rerooting is enabled. Use --keep-root to keep the input root.
  #[cfg_attr(
    feature = "clap",
    clap(
      long = "reroot",
      value_enum,
      conflicts_with = "reroot_tips",
      help_heading = "Rooting"
    )
  )]
  reroot: Option<RerootMethodCli>,

  /// Reroot on the branch leading to a tip or the MRCA of a comma-separated tip list.
  #[cfg_attr(
    feature = "clap",
    clap(
      long = "reroot-tips",
      value_delimiter = ',',
      conflicts_with = "reroot",
      help_heading = "Rooting"
    )
  )]
  reroot_tips: Vec<String>,
}

impl RerootArgs {
  pub fn spec(&self, graph: &Graph, names: &BTreeMap<GraphNodeKey, Option<String>>) -> Result<RerootSpec, Report> {
    if self.reroot_tips.is_empty() {
      Ok(RerootSpec::Method(self.reroot.map(Into::into).unwrap_or_default()))
    } else {
      Ok(RerootSpec::Tips(resolve_reroot_tips(&self.reroot_tips, graph, names)?))
    }
  }
}

pub fn resolve_reroot_tips(
  tips: &[String],
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<Vec<GraphNodeKey>, Report> {
  let first_leaf_of = first_key_by_name(graph.get_leaves().map(|leaf| leaf.key()), names);
  tips
    .iter()
    .map(|tip| {
      first_leaf_of
        .get(tip.as_str())
        .copied()
        .ok_or_else(|| make_report!("Reroot tip not found: {tip}"))
    })
    .collect()
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "RerootMethod")]
pub enum RerootMethodCli {
  #[default]
  #[cfg_attr(feature = "clap", value(alias = "best"))]
  LeastSquares,
  MinDev,
  Oldest,
  ClockFilter,
}

impl From<RerootMethodCli> for RerootMethod {
  fn from(method: RerootMethodCli) -> Self {
    match method {
      RerootMethodCli::LeastSquares => RerootMethod::LeastSquares,
      RerootMethodCli::MinDev => RerootMethod::MinDev,
      RerootMethodCli::Oldest => RerootMethod::Oldest,
      RerootMethodCli::ClockFilter => RerootMethod::ClockFilter,
    }
  }
}

#[cfg(test)]
mod tests {
  use crate::commands::shared::reroot::{RerootArgs, RerootMethodCli, resolve_reroot_tips};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::clock::find_best_root::params::{RerootMethod, RerootSpec};
  use treetime::o;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::assert_error;

  const TREE: &str = "((A:0.1,B:0.1)AB:0.1,(A:0.1,C:0.1)AC:0.1)root;";

  #[test]
  fn test_reroot_args_default_spec_is_least_squares() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let args = RerootArgs::default();

    let actual = args.spec(&tree.graph, &tree.names()).unwrap();

    assert_eq!(RerootSpec::Method(RerootMethod::LeastSquares), actual);
  }

  #[test]
  fn test_reroot_args_method_spec() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let args = RerootArgs {
      reroot: Some(RerootMethodCli::MinDev),
      ..RerootArgs::default()
    };

    let actual = args.spec(&tree.graph, &tree.names()).unwrap();

    assert_eq!(RerootSpec::Method(RerootMethod::MinDev), actual);
  }

  #[test]
  fn test_reroot_args_tips_spec_resolves_leaf_keys() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let args = RerootArgs {
      reroot_tips: vec![o!("B"), o!("C")],
      ..RerootArgs::default()
    };

    let actual = args.spec(&tree.graph, &tree.names()).unwrap();

    let expected = vec![helpers::leaf_keys(&tree, "B")[0], helpers::leaf_keys(&tree, "C")[0]];
    assert_eq!(RerootSpec::Tips(expected), actual);
  }

  #[test]
  fn test_reroot_resolve_tips_takes_the_first_leaf_of_a_repeated_name_in_input_order() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let leaves_named_a = helpers::leaf_keys(&tree, "A");

    let actual = resolve_reroot_tips(&[o!("C"), o!("A")], &tree.graph, &tree.names()).unwrap();

    assert_eq!(2, leaves_named_a.len());
    assert_eq!(vec![helpers::leaf_keys(&tree, "C")[0], leaves_named_a[0]], actual);
  }

  #[test]
  fn test_reroot_resolve_tips_of_an_empty_list_is_empty() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();

    let actual = resolve_reroot_tips(&[], &tree.graph, &tree.names()).unwrap();

    assert_eq!(Vec::<GraphNodeKey>::new(), actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::missing_tip(        "Z",  "Reroot tip not found: Z")]
  #[case::internal_node_name( "AB", "Reroot tip not found: AB")]
  #[trace]
  fn test_reroot_resolve_tips_rejects_names_of_no_leaf(#[case] tip: &str, #[case] expected: &str) {
    let tree = nwk_read(TREE.as_bytes()).unwrap();

    let result = resolve_reroot_tips(&[tip.to_owned()], &tree.graph, &tree.names());

    assert_error!(result, expected);
  }

  mod helpers {
    use itertools::Itertools;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::NwkParse;

    pub(super) fn leaf_keys(tree: &NwkParse, name: &str) -> Vec<GraphNodeKey> {
      let names = tree.names();
      tree
        .graph
        .get_leaves()
        .map(|leaf| leaf.key())
        .filter(|key| names[key].as_deref() == Some(name))
        .sorted()
        .collect()
    }
  }
}
