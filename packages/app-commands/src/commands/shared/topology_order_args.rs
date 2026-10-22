use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::leaf_order::leaf_order;
#[cfg(feature = "clap")]
use clap::ValueHint;
use deser::{Deserialize, Serialize};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::pair_by_name::pair_by_name;
use treetime_graph::topology_order::{TopologyOrderPreset, TopologyOrderSpec, TopologyOrderTargetAggregate};
use treetime_io::name_list::name_list_read_file;
use treetime_io::nwk::NewickDialect;
use treetime_io::tree::tree_read_file;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

#[derive(Debug, Clone, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct TopologyOrderArgs {
  /// Canonical topology ordering preset.
  #[cfg_attr(feature = "clap", clap(long, value_enum, help_heading = "Tree ordering"))]
  pub topology_order: Option<TopologyOrderArg>,

  /// Tree whose leaf order sets the target order
  ///
  /// Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-list`,
  /// the leaf order of the input tree is the target order.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_hint = ValueHint::FilePath, value_name = "FILE", conflicts_with = "topology_order_list", help_heading = "Tree ordering")
  )]
  #[schemars(extend("x-path" = "input"))]
  pub topology_order_tree: Option<PathBuf>,

  /// File with leaf names in the target order, one name per line
  ///
  /// Used by `--topology-order target-order` and `target-order-reverse`. Without this flag and `--topology-order-tree`,
  /// the leaf order of the input tree is the target order.
  #[cfg_attr(
    feature = "clap",
    clap(long, value_hint = ValueHint::FilePath, value_name = "FILE", help_heading = "Tree ordering")
  )]
  #[schemars(extend("x-path" = "input"))]
  pub topology_order_list: Option<PathBuf>,

  /// Aggregate used to map a subtree to a target-order position.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = TopologyOrderTargetAggregateArg::default(), help_heading = "Tree ordering"))]
  #[default(TopologyOrderTargetAggregateArg::Mean)]
  pub topology_order_aggregate: TopologyOrderTargetAggregateArg,
}

#[derive(Debug, Clone)]
pub(crate) struct TopologyOrder {
  preset: TopologyOrderPreset,
  target: TargetOrderSource,
  aggregate: TopologyOrderTargetAggregate,
}

impl TopologyOrder {
  pub(crate) fn new(
    topology_order: Option<TopologyOrderArg>,
    topology_order_tree: Option<PathBuf>,
    topology_order_list: Option<PathBuf>,
    topology_order_aggregate: TopologyOrderTargetAggregateArg,
  ) -> Result<Self, Report> {
    let preset = topology_order.map_or(TopologyOrderPreset::DescendantCount, Into::into);
    let target_fields_present = topology_order_tree.is_some()
      || topology_order_list.is_some()
      || topology_order_aggregate != TopologyOrderTargetAggregateArg::Mean;
    if target_fields_present && !preset.is_target_order() {
      return Err(argument_error(
        &[
          "topology_order",
          "topology_order_tree",
          "topology_order_list",
          "topology_order_aggregate",
        ],
        "--topology-order-tree, --topology-order-list and --topology-order-aggregate require --topology-order=target-order or target-order-reverse",
      ));
    }
    let target = match (topology_order_tree, topology_order_list) {
      (Some(path), _) => TargetOrderSource::Tree(path),
      (None, Some(path)) => TargetOrderSource::List(path),
      (None, None) => TargetOrderSource::InputTree,
    };
    Ok(Self {
      preset,
      target,
      aggregate: topology_order_aggregate.into(),
    })
  }

  pub(crate) fn resolve(
    &self,
    graph: &Graph,
    names: &NodeNames,
    input_order: Option<Vec<GraphNodeKey>>,
    tree_dialect: NewickDialect,
  ) -> Result<TopologyOrderSpec, Report> {
    let target_order = if self.preset.is_target_order() {
      self.target_order(graph, names, input_order, tree_dialect)?
    } else {
      BTreeMap::new()
    };
    Ok(TopologyOrderSpec {
      preset: self.preset,
      target_order,
      target_aggregate: self.aggregate,
    })
  }

  fn target_order(
    &self,
    graph: &Graph,
    names: &NodeNames,
    input_order: Option<Vec<GraphNodeKey>>,
    tree_dialect: NewickDialect,
  ) -> Result<BTreeMap<GraphNodeKey, usize>, Report> {
    let target_names = match &self.target {
      TargetOrderSource::Tree(path) => {
        let nwk_parsed = tree_read_file(path, tree_dialect).wrap_err("When reading the topology order tree")?;
        leaf_order(&nwk_parsed.graph)
          .into_iter()
          .map(|key| nwk_parsed.names.pairing_name(key).to_owned())
          .collect_vec()
      },
      TargetOrderSource::List(path) => name_list_read_file(path).wrap_err("When reading the topology order list")?,
      TargetOrderSource::InputTree => {
        let order = input_order.unwrap_or_else(|| leaf_order(graph));
        return Ok(
          order
            .into_iter()
            .enumerate()
            .map(|(position, key)| (key, position))
            .collect(),
        );
      },
    };
    Ok(target_positions(target_names, graph, names))
  }
}

#[derive(Debug, Clone)]
enum TargetOrderSource {
  InputTree,
  Tree(PathBuf),
  List(PathBuf),
}

pub(crate) fn target_positions(
  target_names: Vec<String>,
  graph: &Graph,
  names: &NodeNames,
) -> BTreeMap<GraphNodeKey, usize> {
  let positions = target_names
    .into_iter()
    .enumerate()
    .map(|(position, name)| (name, position));
  pair_by_name(graph.get_leaves().map(|leaf| leaf.key()), names, positions).by_node
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum TopologyOrderArg {
  Keep,
  DescendantCount,
  DescendantCountReverse,
  Height,
  HeightReverse,
  Divergence,
  DivergenceReverse,
  Label,
  LabelReverse,
  TargetOrder,
  TargetOrderReverse,
}

impl From<TopologyOrderArg> for TopologyOrderPreset {
  fn from(value: TopologyOrderArg) -> Self {
    match value {
      TopologyOrderArg::Keep => Self::Keep,
      TopologyOrderArg::DescendantCount => Self::DescendantCount,
      TopologyOrderArg::DescendantCountReverse => Self::DescendantCountReverse,
      TopologyOrderArg::Height => Self::Height,
      TopologyOrderArg::HeightReverse => Self::HeightReverse,
      TopologyOrderArg::Divergence => Self::Divergence,
      TopologyOrderArg::DivergenceReverse => Self::DivergenceReverse,
      TopologyOrderArg::Label => Self::Label,
      TopologyOrderArg::LabelReverse => Self::LabelReverse,
      TopologyOrderArg::TargetOrder => Self::TargetOrder,
      TopologyOrderArg::TargetOrderReverse => Self::TargetOrderReverse,
    }
  }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum TopologyOrderTargetAggregateArg {
  #[default]
  Mean,
  Median,
}

impl From<TopologyOrderTargetAggregateArg> for TopologyOrderTargetAggregate {
  fn from(value: TopologyOrderTargetAggregateArg) -> Self {
    match value {
      TopologyOrderTargetAggregateArg::Mean => Self::Mean,
      TopologyOrderTargetAggregateArg::Median => Self::Median,
    }
  }
}
