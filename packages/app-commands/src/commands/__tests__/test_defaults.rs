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
  use rstest::rstest;
  use serde_json::Value;

  #[rustfmt::skip]
  #[rstest]
  #[case::ancestral( helpers::config_and_clap_defaults::<TreetimeAncestralArgsRaw>())]
  #[case::clock(     helpers::config_and_clap_defaults::<TreetimeClockArgsRaw>())]
  #[case::homoplasy( helpers::config_and_clap_defaults::<TreetimeHomoplasyArgsRaw>())]
  #[case::mugration( helpers::config_and_clap_defaults::<TreetimeMugrationArgsRaw>())]
  #[case::optimize(  helpers::config_and_clap_defaults::<TreetimeOptimizeArgsRaw>())]
  #[case::prune(     helpers::config_and_clap_defaults::<TreetimePruneArgsRaw>())]
  #[case::timetree(  helpers::config_and_clap_defaults::<TreetimeTimetreeArgsRaw>())]
  #[trace]
  fn test_defaults_clap_matches_config(#[case] (config, clap): (Value, Value)) {
    assert_eq!(config, clap);
  }

  mod helpers {
    use clap::Parser;
    use deser::Serialize;
    use serde_json::Value;
    use treetime_utils::io::json::to_json_value;

    pub(super) fn config_and_clap_defaults<T: Default + Parser + Serialize>() -> (Value, Value) {
      let config = to_json_value(&T::default()).unwrap();
      let clap = to_json_value(&T::try_parse_from(["treetime"]).unwrap()).unwrap();
      (config, clap)
    }
  }
}
