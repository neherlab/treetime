use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use treetime::clock::find_best_root::params::{RerootMethod, RerootSpec};

#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub(crate) struct RerootArgs {
  /// Reroot the tree by temporal-signal optimization.
  ///
  /// Defaults to least-squares when rerooting is enabled. Use --keep-root to keep the input root.
  #[cfg_attr(feature = "clap", clap(long = "reroot", value_enum, conflicts_with = "reroot_tips"))]
  reroot: Option<RerootMethodCli>,

  /// Reroot on the branch leading to a tip or the MRCA of a comma-separated tip list.
  #[cfg_attr(
    feature = "clap",
    clap(long = "reroot-tips", value_delimiter = ',', conflicts_with = "reroot")
  )]
  reroot_tips: Vec<String>,
}

impl RerootArgs {
  pub(crate) fn spec(&self) -> RerootSpec {
    if self.reroot_tips.is_empty() {
      RerootSpec::Method(self.reroot.map(Into::into).unwrap_or_default())
    } else {
      RerootSpec::Tips(self.reroot_tips.clone())
    }
  }
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, SmartDefault, Serialize, Deserialize, JsonSchema)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[serde(rename_all = "kebab-case")]
#[schemars(rename = "RerootMethod")]
pub(crate) enum RerootMethodCli {
  #[default]
  #[cfg_attr(feature = "clap", value(alias = "best"))]
  LeastSquares,
  MinDev,
  Oldest,
  #[cfg_attr(feature = "clap", value(alias = "clock-filter"))]
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
  use crate::commands::shared::reroot::{RerootArgs, RerootMethodCli};
  use pretty_assertions::assert_eq;
  use treetime::clock::find_best_root::params::{RerootMethod, RerootSpec};
  use treetime::o;

  #[test]
  fn test_reroot_args_default_spec_is_least_squares() {
    let args = RerootArgs::default();

    let actual = args.spec();

    assert_eq!(RerootSpec::Method(RerootMethod::LeastSquares), actual);
  }

  #[test]
  fn test_reroot_args_method_spec() {
    let args = RerootArgs {
      reroot: Some(RerootMethodCli::MinDev),
      ..RerootArgs::default()
    };

    let actual = args.spec();

    assert_eq!(RerootSpec::Method(RerootMethod::MinDev), actual);
  }

  #[test]
  fn test_reroot_args_tips_spec() {
    let args = RerootArgs {
      reroot_tips: vec![o!("A"), o!("B")],
      ..RerootArgs::default()
    };

    let actual = args.spec();

    assert_eq!(RerootSpec::Tips(vec![o!("A"), o!("B")]), actual);
  }
}
