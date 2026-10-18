use crate::optimize::indel::IndelPolicy;
use crate::optimize::partition_rates::PartitionRates;
use smart_default::SmartDefault;
use strum_macros::Display;

#[derive(Copy, Clone, Debug, PartialEq, Eq, SmartDefault)]
pub struct TopologyOps {
  #[default = true]
  pub collapse_short_branches: bool,
  #[default = true]
  pub merge_siblings: bool,
  #[default = true]
  pub flip_parent_child: bool,
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault)]
pub enum BranchLengthMode {
  Input,
  #[default]
  Marginal,
}

#[derive(Copy, Clone, Debug, Display, PartialEq, Eq, SmartDefault)]
#[strum(serialize_all = "kebab-case")]
pub enum BranchOptMethod {
  #[default]
  Brent,
  Newton,
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum InitialGuessMode {
  #[default]
  Auto,
  Always,
  Never,
}

#[derive(Copy, Debug, Clone, PartialEq, Eq)]
pub enum ExistingBranchLengths {
  Keep,
  Overwrite,
}

#[derive(Clone, Copy, Debug)]
pub struct OptimizeLoopSettings {
  pub max_iter: usize,
  pub dp: f64,
  pub damping: f64,
  pub opt_method: BranchOptMethod,
  pub indel_policy: IndelPolicy,
  pub topology_ops: TopologyOps,
  pub partition_rates: PartitionRates,
}
