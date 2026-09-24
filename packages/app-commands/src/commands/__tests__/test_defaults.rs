#[cfg(test)]
mod tests {
  use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
  use crate::commands::clock::args::TreetimeClockArgsRaw;
  use crate::commands::homoplasy::args::TreetimeHomoplasyArgsRaw;
  use crate::commands::mugration::args::TreetimeMugrationArgsRaw;
  use crate::commands::optimize::args::TreetimeOptimizeArgsRaw;
  use crate::commands::prune::args::TreetimePruneArgsRaw;
  use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_defaults_ancestral_clap_matches_config() {
    assert_eq!(
      helpers::config_defaults::<TreetimeAncestralArgsRaw>(),
      helpers::clap_defaults::<TreetimeAncestralArgsRaw>()
    );
  }

  #[test]
  fn test_defaults_clock_clap_matches_config() {
    assert_eq!(
      helpers::config_defaults::<TreetimeClockArgsRaw>(),
      helpers::clap_defaults::<TreetimeClockArgsRaw>()
    );
  }

  #[test]
  fn test_defaults_homoplasy_clap_matches_config() {
    assert_eq!(
      helpers::config_defaults::<TreetimeHomoplasyArgsRaw>(),
      helpers::clap_defaults::<TreetimeHomoplasyArgsRaw>()
    );
  }

  #[test]
  fn test_defaults_mugration_clap_matches_config() {
    assert_eq!(
      helpers::config_defaults::<TreetimeMugrationArgsRaw>(),
      helpers::clap_defaults::<TreetimeMugrationArgsRaw>()
    );
  }

  #[test]
  fn test_defaults_optimize_clap_matches_config() {
    assert_eq!(
      helpers::config_defaults::<TreetimeOptimizeArgsRaw>(),
      helpers::clap_defaults::<TreetimeOptimizeArgsRaw>()
    );
  }

  #[test]
  fn test_defaults_prune_clap_matches_config() {
    assert_eq!(
      helpers::config_defaults::<TreetimePruneArgsRaw>(),
      helpers::clap_defaults::<TreetimePruneArgsRaw>()
    );
  }

  #[test]
  fn test_defaults_timetree_clap_matches_config() {
    assert_eq!(
      helpers::config_defaults::<TreetimeTimetreeArgsRaw>(),
      helpers::clap_defaults::<TreetimeTimetreeArgsRaw>()
    );
  }

  mod helpers {
    use clap::Parser;
    use serde::Serialize;
    use serde_json::Value;

    pub(super) fn config_defaults<T: Default + Serialize>() -> Value {
      serde_json::to_value(T::default()).unwrap()
    }

    pub(super) fn clap_defaults<T: Parser + Serialize>() -> Value {
      serde_json::to_value(T::try_parse_from(["treetime"]).unwrap()).unwrap()
    }
  }
}
