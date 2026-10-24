use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::mirror_enum::impl_from_mirror;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use treetime::reroot::params::{RerootMethod, RerootSpec};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::pair_by_name::key_by_name;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct RerootArgs {
  /// Method that chooses the root of the tree
  ///
  /// `least-squares` and `oldest` use the sampling dates; `min-dev` and `midpoint` use only the branch lengths, root
  /// once on the input branch lengths, and count every tip. Without a rooting setting the tree is rerooted by
  /// least-squares. Use --outgroup to root on named tips, or --keep-root to keep the input root.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_enum, conflicts_with = "outgroup", help_heading = "Rooting")
  )]
  pub(crate) reroot: Option<RerootMethodCli>,
}

#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct RootingArgs {
  /// Root the tree on the branch above the most recent common ancestor of the named tips
  ///
  /// With one tip, the root is placed on the branch leading to that tip.
  #[schemars(example = ["outgroup_tip"])]
  #[cfg_attr(
    feature = "clap",
    clap(long, num_args = 1.., value_delimiter = ',', conflicts_with = "reroot", help_heading = "Rooting")
  )]
  pub outgroup: Vec<String>,

  /// Keep the root of the input tree instead of rerooting it.
  #[cfg_attr(feature = "clap", clap(long, conflicts_with_all = ["reroot", "outgroup"], help_heading = "Rooting"))]
  pub keep_root: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RerootTarget {
  Method(RerootMethod),
  Outgroup(Vec<String>),
}

impl RerootTarget {
  pub fn unless_keep_root(keep_root: bool, reroot: Option<RerootMethod>, outgroup: Vec<String>) -> Option<Self> {
    (!keep_root).then(|| Self::new(reroot, outgroup))
  }

  pub fn new(reroot: Option<RerootMethod>, outgroup: Vec<String>) -> Self {
    if outgroup.is_empty() {
      Self::Method(reroot.unwrap_or_default())
    } else {
      Self::Outgroup(outgroup)
    }
  }

  pub fn spec(&self, graph: &Graph, names: &NodeNames) -> Result<RerootSpec, Report> {
    match self {
      Self::Method(method) => Ok(RerootSpec::Method(*method)),
      Self::Outgroup(tips) => Ok(RerootSpec::Tips(resolve_outgroup(tips, graph, names)?)),
    }
  }
}

pub fn resolve_outgroup(tips: &[String], graph: &Graph, names: &NodeNames) -> Result<Vec<GraphNodeKey>, Report> {
  let leaf_of = key_by_name(graph.get_leaves().map(|leaf| leaf.key()), names);
  tips
    .iter()
    .map(|tip| {
      leaf_of
        .get(tip.as_str())
        .copied()
        .ok_or_else(|| argument_error(&["outgroup"], format!("Outgroup tip not found: {tip}")))
    })
    .collect()
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[schemars(rename = "RerootMethod")]
pub enum RerootMethodCli {
  /// Minimizes the residuals of the root-to-tip regression on the sampling dates
  #[default]
  LeastSquares,
  /// Minimizes the variance of root-to-tip distances over all tips; more robust than midpoint
  MinDev,
  /// Roots at the middle of the branch above the oldest dated tip
  Oldest,
  /// Roots at the middle of the longest tip-to-tip path; uses only the two farthest tips
  Midpoint,
}

impl_from_mirror!(RerootMethodCli => RerootMethod { LeastSquares, MinDev, Oldest, Midpoint });

#[cfg(test)]
mod tests {
  use crate::commands::shared::reroot::{RerootTarget, resolve_outgroup};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::o;
  use treetime::reroot::params::{RerootMethod, RerootSpec};
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::assert_error;

  const TREE: &str = "((A:0.1,B:0.1)AB:0.1,(A:0.1,C:0.1)AC:0.1)root;";

  #[test]
  fn test_reroot_target_default_spec_is_least_squares() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();

    let actual = RerootTarget::new(None, vec![]).spec(&tree.graph, &tree.names).unwrap();

    assert_eq!(RerootSpec::Method(RerootMethod::LeastSquares), actual);
  }

  #[test]
  fn test_reroot_target_method_spec() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();

    let actual = RerootTarget::new(Some(RerootMethod::MinDev), vec![])
      .spec(&tree.graph, &tree.names)
      .unwrap();

    assert_eq!(RerootSpec::Method(RerootMethod::MinDev), actual);
  }

  #[test]
  fn test_reroot_target_outgroup_spec_resolves_leaf_keys() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();

    let actual = RerootTarget::new(None, vec![o!("B"), o!("C")])
      .spec(&tree.graph, &tree.names)
      .unwrap();

    let expected = vec![helpers::leaf_keys(&tree, "B")[0], helpers::leaf_keys(&tree, "C")[0]];
    assert_eq!(RerootSpec::Tips(expected), actual);
  }

  #[test]
  fn test_reroot_resolve_outgroup_takes_the_leaf_of_a_repeated_name_in_input_order() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let leaves_named_a = helpers::leaf_keys(&tree, "A");

    let actual = resolve_outgroup(&[o!("C"), o!("A")], &tree.graph, &tree.names).unwrap();

    assert_eq!(2, leaves_named_a.len());
    assert_eq!(vec![helpers::leaf_keys(&tree, "C")[0], leaves_named_a[0]], actual);
  }

  #[test]
  fn test_reroot_resolve_outgroup_of_an_empty_list_is_empty() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();

    let actual = resolve_outgroup(&[], &tree.graph, &tree.names).unwrap();

    assert_eq!(Vec::<GraphNodeKey>::new(), actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::missing_tip(        "Z",  "Outgroup tip not found: Z")]
  #[case::internal_node_name( "AB", "Outgroup tip not found: AB")]
  #[trace]
  fn test_reroot_resolve_outgroup_rejects_names_of_no_leaf(#[case] tip: &str, #[case] expected: &str) {
    let tree = nwk_read(TREE.as_bytes()).unwrap();

    let result = resolve_outgroup(&[tip.to_owned()], &tree.graph, &tree.names);

    assert_error!(result, expected);
  }

  mod helpers {
    use itertools::Itertools;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::NwkParse;

    pub(super) fn leaf_keys(tree: &NwkParse, name: &str) -> Vec<GraphNodeKey> {
      tree
        .graph
        .get_leaves()
        .map(|leaf| leaf.key())
        .filter(|key| tree.names.pairing_name(*key) == name)
        .sorted()
        .collect()
    }
  }
}
