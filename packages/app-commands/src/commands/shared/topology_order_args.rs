use crate::commands::shared::leaf_order::leaf_order;
#[cfg(feature = "clap")]
use clap::ValueHint;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_with::skip_serializing_none;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pair_by_name::pair_by_name;
use treetime_graph::topology_order::{TopologyOrderPreset, TopologyOrderSpec, TopologyOrderTargetAggregate};
use treetime_io::name_list::name_list_read_file;
use treetime_io::tree::tree_read_file;
use treetime_utils::{make_error, make_report};

#[skip_serializing_none]
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
#[deser(skip_serializing_optionals)]
#[serde(default, deny_unknown_fields)]
#[deser(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct TopologyOrderArgs {
  /// Order tree topology before writing output files.
  #[cfg_attr(feature = "clap", clap(long, value_enum, help_heading = "Tree ordering"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order_target_source"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order_target_file"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order_target_aggregate"))]
  pub ladderize: Option<LadderizeArg>,

  /// Canonical topology ordering preset.
  #[cfg_attr(feature = "clap", clap(long, value_enum, help_heading = "Tree ordering"))]
  pub topology_order: Option<TopologyOrderArg>,

  /// Source for target-order topology sorting.
  #[cfg_attr(feature = "clap", clap(long, value_enum, help_heading = "Tree ordering"))]
  pub topology_order_target_source: Option<TopologyOrderTargetSourceArg>,

  /// File used by list or reference-topology target-order sources.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath, help_heading = "Tree ordering"))]
  #[schemars(extend("x-path" = "input"))]
  pub topology_order_target_file: Option<PathBuf>,

  /// Aggregate used to map a subtree to a target-order position.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = TopologyOrderTargetAggregateArg::default(), help_heading = "Tree ordering"))]
  #[default(TopologyOrderTargetAggregateArg::Mean)]
  pub topology_order_target_aggregate: TopologyOrderTargetAggregateArg,
}

impl TopologyOrderArgs {
  pub fn resolve_topology_order(
    &self,
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    input_order: Option<Vec<GraphNodeKey>>,
  ) -> Result<TopologyOrderSpec, Report> {
    self.validate()?;

    let preset = match (self.ladderize, self.topology_order) {
      (Some(LadderizeArg::None), None) => TopologyOrderPreset::Keep,
      (Some(LadderizeArg::Ascending), None) => TopologyOrderPreset::DescendantCount,
      (Some(LadderizeArg::Descending), None) => TopologyOrderPreset::DescendantCountReverse,
      (None, Some(topology_order)) => topology_order.into(),
      (None, None) => TopologyOrderPreset::DescendantCount,
      (Some(_), Some(_)) => {
        return make_error!("--ladderize cannot be combined with --topology-order");
      },
    };

    let target_order = if preset.is_target_order() {
      self.target_order(graph, names, input_order)?
    } else {
      BTreeMap::new()
    };

    Ok(TopologyOrderSpec {
      preset,
      target_order,
      target_aggregate: self.topology_order_target_aggregate.into(),
    })
  }

  fn validate(&self) -> Result<(), Report> {
    if self.ladderize.is_some()
      && (self.topology_order.is_some()
        || self.topology_order_target_source.is_some()
        || self.topology_order_target_file.is_some()
        || self.topology_order_target_aggregate != TopologyOrderTargetAggregateArg::Mean)
    {
      return make_error!("--ladderize cannot be combined with --topology-order* options");
    }

    let target_fields_present = self.topology_order_target_source.is_some()
      || self.topology_order_target_file.is_some()
      || self.topology_order_target_aggregate != TopologyOrderTargetAggregateArg::Mean;
    let target_mode = self
      .topology_order
      .is_some_and(|order| TopologyOrderPreset::from(order).is_target_order());

    if target_fields_present && !target_mode {
      return make_error!(
        "--topology-order-target-* options require --topology-order=target-order or target-order-reverse"
      );
    }

    if target_mode
      && matches!(
        self.topology_order_target_source,
        Some(TopologyOrderTargetSourceArg::ReferenceTopology | TopologyOrderTargetSourceArg::List)
      )
      && self.topology_order_target_file.is_none()
    {
      return make_error!("--topology-order-target-file is required for reference-topology and list target sources");
    }

    Ok(())
  }

  fn target_order(
    &self,
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    input_order: Option<Vec<GraphNodeKey>>,
  ) -> Result<BTreeMap<GraphNodeKey, usize>, Report> {
    let target_names = match self
      .topology_order_target_source
      .unwrap_or(TopologyOrderTargetSourceArg::Input)
    {
      TopologyOrderTargetSourceArg::Input => {
        let order = input_order.unwrap_or_else(|| leaf_order(graph));
        return Ok(
          order
            .into_iter()
            .enumerate()
            .map(|(position, key)| (key, position))
            .collect(),
        );
      },
      TopologyOrderTargetSourceArg::ReferenceTopology => {
        let path = self
          .topology_order_target_file
          .as_ref()
          .ok_or_else(|| make_report!("--topology-order-target-file is required for reference-topology"))?;
        let nwk_parsed = tree_read_file(path).wrap_err("When reading target reference topology")?;
        let ref_names = nwk_parsed.names();
        leaf_order(&nwk_parsed.graph)
          .into_iter()
          .filter_map(|key| ref_names[&key].clone())
          .collect_vec()
      },
      TopologyOrderTargetSourceArg::List => {
        let path = self
          .topology_order_target_file
          .as_ref()
          .ok_or_else(|| make_report!("--topology-order-target-file is required for list"))?;
        name_list_read_file(path, b'\n').wrap_err("When reading the topology order target list")?
      },
    };
    Ok(target_positions(target_names, graph, names))
  }
}

pub(crate) fn target_positions(
  target_names: Vec<String>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> BTreeMap<GraphNodeKey, usize> {
  let positions = target_names
    .into_iter()
    .enumerate()
    .map(|(position, name)| (name, position));
  pair_by_name(graph.get_leaves().map(|leaf| leaf.key()), names, positions).by_node
}

#[derive(
  Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize,
)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum LadderizeArg {
  None,
  Ascending,
  Descending,
}

#[derive(
  Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize,
)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum TopologyOrderArg {
  Keep,
  #[cfg_attr(feature = "clap", value(alias = "ladderize"))]
  DescendantCount,
  #[cfg_attr(feature = "clap", value(alias = "ladderize-reverse"))]
  DescendantCountReverse,
  Height,
  HeightReverse,
  Divergence,
  DivergenceReverse,
  Label,
  #[cfg_attr(feature = "clap", value(alias = "alphabetical-reverse"))]
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

#[derive(
  Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize,
)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum TopologyOrderTargetSourceArg {
  Input,
  ReferenceTopology,
  List,
}

#[derive(
  Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize,
)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
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
