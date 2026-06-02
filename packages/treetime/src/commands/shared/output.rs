#[cfg(feature = "clap")]
use clap::ValueHint;
use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::path::{Path, PathBuf};
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_graph::topology_order::{TopologyOrderPreset, TopologyOrderSpec, TopologyOrderTargetAggregate};
use treetime_io::graph::GraphWriteOptions;
use treetime_io::nwk::{EdgeFromNwk, NodeFromNwk, nwk_read_file};
use treetime_utils::{make_error, make_report};

/// Output destination shared by every command.
///
/// `--outdir` (short `-O`) names a directory into which the standard set of output files is written
/// under fixed names. Each command also exposes per-type `--output-*` flags that override the path of
/// a single output. An override path of `-` writes that output to standard output;
/// compression is chosen from the path extension by the underlying writer. This mirrors the Nextclade
/// scheme of explicit per-type output flags layered over a directory default: type is selected by the
/// flag, never inferred from the file extension.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct OutputArgs {
  /// Directory to write the standard set of output files to
  #[cfg_attr(feature = "clap", clap(long, short = 'O', value_hint = ValueHint::AnyPath))]
  pub outdir: PathBuf,

  #[cfg_attr(feature = "clap", clap(flatten))]
  pub topology_order_args: TopologyOrderArgs,
}

impl OutputArgs {
  /// Resolve a single output to its destination path.
  ///
  /// Returns the explicit per-type override when given, otherwise `outdir/<default_name>`.
  pub fn resolve(&self, override_path: Option<&Path>, default_name: &str) -> PathBuf {
    match override_path {
      Some(path) => path.to_path_buf(),
      None => self.outdir.join(default_name),
    }
  }

  pub fn graph_write_options<N, E, D>(&self, graph: &Graph<N, E, D>) -> Result<GraphWriteOptions, Report>
  where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Sync + Send,
  {
    self.topology_order_args.graph_write_options(graph, None)
  }

  pub fn graph_write_options_with_input_order<N, E, D>(
    &self,
    graph: &Graph<N, E, D>,
    input_order: Vec<String>,
  ) -> Result<GraphWriteOptions, Report>
  where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Sync + Send,
  {
    self.topology_order_args.graph_write_options(graph, Some(input_order))
  }
}

#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct TopologyOrderArgs {
  /// Order tree topology before writing output files.
  #[cfg_attr(feature = "clap", clap(long, value_enum))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order_target_source"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order_target_file"))]
  #[cfg_attr(feature = "clap", clap(conflicts_with = "topology_order_target_aggregate"))]
  pub ladderize: Option<LadderizeArg>,

  /// Canonical topology ordering preset.
  #[cfg_attr(feature = "clap", clap(long, value_enum))]
  pub topology_order: Option<TopologyOrderArg>,

  /// Source for target-order topology sorting.
  #[cfg_attr(feature = "clap", clap(long, value_enum))]
  pub topology_order_target_source: Option<TopologyOrderTargetSourceArg>,

  /// File used by list or reference-topology target-order sources.
  #[cfg_attr(feature = "clap", clap(long, value_hint = ValueHint::FilePath))]
  pub topology_order_target_file: Option<PathBuf>,

  /// Aggregate used to map a subtree to a target-order position.
  #[cfg_attr(feature = "clap", clap(long, value_enum, default_value_t = TopologyOrderTargetAggregateArg::default()))]
  #[default(TopologyOrderTargetAggregateArg::Mean)]
  pub topology_order_target_aggregate: TopologyOrderTargetAggregateArg,
}

impl TopologyOrderArgs {
  fn graph_write_options<N, E, D>(
    &self,
    graph: &Graph<N, E, D>,
    input_order: Option<Vec<String>>,
  ) -> Result<GraphWriteOptions, Report>
  where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Sync + Send,
  {
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
      self.target_order(graph, input_order)?
    } else {
      vec![]
    };

    Ok(GraphWriteOptions {
      topology_order: TopologyOrderSpec {
        preset,
        target_order,
        target_aggregate: self.topology_order_target_aggregate.into(),
      },
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

  fn target_order<N, E, D>(
    &self,
    graph: &Graph<N, E, D>,
    input_order: Option<Vec<String>>,
  ) -> Result<Vec<String>, Report>
  where
    N: GraphNode + Named,
    E: GraphEdge,
    D: Sync + Send,
  {
    match self
      .topology_order_target_source
      .unwrap_or(TopologyOrderTargetSourceArg::Input)
    {
      TopologyOrderTargetSourceArg::Input => input_order.map_or_else(|| leaf_order(graph), Ok),
      TopologyOrderTargetSourceArg::ReferenceTopology => {
        let path = self
          .topology_order_target_file
          .as_ref()
          .ok_or_else(|| make_report!("--topology-order-target-file is required for reference-topology"))?;
        let graph =
          nwk_read_file::<OrderNode, OrderEdge, ()>(path).wrap_err("When reading target reference topology")?;
        leaf_order(&graph)
      },
      TopologyOrderTargetSourceArg::List => {
        let path = self
          .topology_order_target_file
          .as_ref()
          .ok_or_else(|| make_report!("--topology-order-target-file is required for list"))?;
        let contents = std::fs::read_to_string(path)
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

#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
pub enum LadderizeArg {
  None,
  Ascending,
  Descending,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
pub enum TopologyOrderArg {
  Keep,
  #[cfg_attr(feature = "clap", value(alias = "ladderize"))]
  DescendantCount,
  #[cfg_attr(feature = "clap", value(alias = "ladderize-reverse"))]
  DescendantCountReverse,
  Height,
  HeightReverse,
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
      TopologyOrderArg::Label => Self::Label,
      TopologyOrderArg::LabelReverse => Self::LabelReverse,
      TopologyOrderArg::TargetOrder => Self::TargetOrder,
      TopologyOrderArg::TargetOrderReverse => Self::TargetOrderReverse,
    }
  }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
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

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
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

fn leaf_order<N, E, D>(graph: &Graph<N, E, D>) -> Result<Vec<String>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Sync + Send,
{
  graph
    .get_leaves()
    .into_iter()
    .map(|leaf| {
      let leaf = leaf.read_arc();
      leaf
        .payload()
        .read_arc()
        .name()
        .map(|name| name.as_ref().to_owned())
        .ok_or_else(|| make_report!("Leaf node {} has no name", leaf.key()))
    })
    .collect()
}

#[derive(Clone, Debug, Default)]
struct OrderNode {
  name: Option<String>,
}

impl GraphNode for OrderNode {}

impl Named for OrderNode {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|name| name.as_ref().to_owned());
  }
}

impl NodeFromNwk for OrderNode {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|name| name.as_ref().to_owned()),
    })
  }
}

#[derive(Clone, Debug, Default)]
struct OrderEdge {
  branch_length: Option<f64>,
}

impl GraphEdge for OrderEdge {}

impl HasBranchLength for OrderEdge {
  fn branch_length(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_branch_length(&mut self, branch_length: Option<f64>) {
    self.branch_length = branch_length;
  }
}

impl EdgeFromNwk for OrderEdge {
  fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
    Ok(Self { branch_length: weight })
  }
}
