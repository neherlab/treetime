use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;

#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct RerootArgs {
  /// Reroot the tree by temporal-signal optimization.
  ///
  /// Defaults to least-squares when rerooting is enabled. Use --keep-root to keep the input root.
  #[cfg_attr(feature = "clap", clap(long = "reroot", value_enum, conflicts_with = "reroot_tips"))]
  pub reroot: Option<RerootMethod>,

  /// Reroot on the branch leading to a tip or the MRCA of a comma-separated tip list.
  #[cfg_attr(
    feature = "clap",
    clap(long = "reroot-tips", value_delimiter = ',', conflicts_with = "reroot")
  )]
  pub reroot_tips: Vec<String>,
}

impl RerootArgs {
  pub fn spec(&self) -> RerootSpec {
    if self.reroot_tips.is_empty() {
      RerootSpec::Method(self.reroot.unwrap_or_default())
    } else {
      RerootSpec::Tips(self.reroot_tips.clone())
    }
  }
}

#[cfg(test)]
mod tests {
  use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
  use crate::commands::shared::reroot::RerootArgs;
  use crate::o;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_reroot_args_default_spec_is_least_squares() {
    let args = RerootArgs::default();

    let actual = args.spec();

    assert_eq!(RerootSpec::Method(RerootMethod::LeastSquares), actual);
  }

  #[test]
  fn test_reroot_args_method_spec() {
    let args = RerootArgs {
      reroot: Some(RerootMethod::MinDev),
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
