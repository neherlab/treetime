use smart_default::SmartDefault;

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

#[derive(Copy, Clone, Debug, PartialEq, Eq, SmartDefault)]
pub enum BranchOptMethod {
  Brent,
  #[default]
  BrentSqrt,
  BrentLog,
  Newton,
  NewtonSqrt,
  NewtonLog,
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
