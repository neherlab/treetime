//! Command-line tree-ordering arguments.
//!
//! Owns the `--ladderize`/`--topology-order*` flags and the reads they imply: a target-order source
//! can name a reference-topology Newick file or a plain list file, and this module reads both. It
//! parses the flags and file contents into a `treetime_graph::topology_order::TopologyOrderSpec`,
//! which the command applies to the graph before the output writers serialize it. Ordering is a
//! graph transform, not an output-encoding concern, so it lives with the adapter, not in `app-output`.

#[cfg(feature = "clap")]
use clap::ValueHint;
use eyre::{Report, WrapErr};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::topology_order::{TopologyOrderPreset, TopologyOrderSpec, TopologyOrderTargetAggregate};
use treetime_io::nwk::nwk_read_file;
use treetime_utils::io::fs::read_file_to_string;
use treetime_utils::{make_error, make_report};

#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
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
    input_order: Option<Vec<String>>,
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
      vec![]
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
      .is_some_and(|order| order.into_preset().is_target_order());

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
    input_order: Option<Vec<String>>,
  ) -> Result<Vec<String>, Report> {
    match self
      .topology_order_target_source
      .unwrap_or(TopologyOrderTargetSourceArg::Input)
    {
      TopologyOrderTargetSourceArg::Input => input_order.map_or_else(|| leaf_order(graph, names), Ok),
      TopologyOrderTargetSourceArg::ReferenceTopology => {
        let path = self
          .topology_order_target_file
          .as_ref()
          .ok_or_else(|| make_report!("--topology-order-target-file is required for reference-topology"))?;
        let nwk_parsed = nwk_read_file(path).wrap_err("When reading target reference topology")?;
        let ref_names = nwk_parsed.names();
        let ref_graph = nwk_parsed.graph;
        leaf_order(&ref_graph, &ref_names)
      },
      TopologyOrderTargetSourceArg::List => {
        let path = self
          .topology_order_target_file
          .as_ref()
          .ok_or_else(|| make_report!("--topology-order-target-file is required for list"))?;
        let contents = read_file_to_string(path)
          .wrap_err_with(|| format!("When reading topology order target list '{}'", path.display()))?;
        Ok(
          contents
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(str::to_owned)
            .collect(),
        )
      },
    }
  }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
pub enum LadderizeArg {
  None,
  Ascending,
  Descending,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
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

impl TopologyOrderArg {
  fn into_preset(self) -> TopologyOrderPreset {
    self.into()
  }
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

#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
pub enum TopologyOrderTargetSourceArg {
  Input,
  ReferenceTopology,
  List,
}

impl std::fmt::Display for TopologyOrderTargetSourceArg {
  fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    match self {
      Self::Input => write!(f, "input"),
      Self::ReferenceTopology => write!(f, "reference-topology"),
      Self::List => write!(f, "list"),
    }
  }
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
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

fn leaf_order(graph: &Graph, names: &BTreeMap<GraphNodeKey, Option<String>>) -> Result<Vec<String>, Report> {
  graph
    .get_leaves()
    .map(|leaf| {
      let key = leaf.key();
      names[&key]
        .clone()
        .ok_or_else(|| make_report!("Leaf node {key} has no name"))
    })
    .collect()
}
