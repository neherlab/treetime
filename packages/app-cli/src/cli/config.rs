use crate::cli::diagnostics::entry::check_command_config;
use crate::cli::diagnostics::source::{ConfigSource, parse_config_document};
use crate::cli::pipeline::types::SCHEMA_KEY;
use crate::cli::schema::command_schema;
use clap::ArgMatches;
use clap::parser::ValueSource;
use eyre::Report;
use schemars::JsonSchema;
use serde::Serialize;
use serde::de::DeserializeOwned;
use serde_json::Value;
use std::collections::BTreeSet;
use std::path::PathBuf;
use treetime_utils::io::fs::read_file_to_string;

pub(crate) fn overlay_config<T>(args: &mut T, matches: &ArgMatches) -> Result<(), Report>
where
  T: Serialize + DeserializeOwned + Default + JsonSchema,
{
  let Some(config_path) = matches.get_one::<PathBuf>("config") else {
    return Ok(());
  };

  let explicit: BTreeSet<String> = matches
    .ids()
    .filter(|id| matches.value_source(id.as_str()) == Some(ValueSource::CommandLine))
    .map(|id| id.to_string())
    .collect();

  let text = read_file_to_string(config_path)?;
  let source = ConfigSource::new(config_path.display().to_string(), text.clone());
  let mut file_value = parse_config_document(&source, &text)?;

  if let Value::Object(map) = &mut file_value {
    map.remove(SCHEMA_KEY);
  }

  let mut merged = serde_json::to_value(T::default())?;
  merge_value(&mut merged, &file_value);
  let cli = serde_json::to_value(&*args)?;
  apply_cli_overrides(&mut merged, &cli, &explicit);

  check_command_config(&source, &merged, &command_schema::<T>())?;

  *args = serde_json::from_value(merged)?;
  Ok(())
}

fn merge_value(base: &mut Value, overlay: &Value) {
  match (base, overlay) {
    (Value::Object(base), Value::Object(overlay)) => {
      for (key, value) in overlay {
        merge_value(base.entry(key.clone()).or_insert(Value::Null), value);
      }
    },
    (base, overlay) => *base = overlay.clone(),
  }
}

fn apply_cli_overrides(merged: &mut Value, cli: &Value, explicit: &BTreeSet<String>) {
  let (Value::Object(merged), Value::Object(cli)) = (merged, cli) else {
    return;
  };
  for (key, cli_child) in cli {
    if explicit.contains(key) {
      merged.insert(key.clone(), cli_child.clone());
    } else if let (Some(merged_child @ Value::Object(_)), Value::Object(_)) = (merged.get_mut(key), cli_child) {
      apply_cli_overrides(merged_child, cli_child, explicit);
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use serde_json::json;

  fn ids(names: &[&str]) -> BTreeSet<String> {
    names.iter().map(|s| (*s).to_owned()).collect()
  }

  #[test]
  fn test_config_apply_cli_overrides_explicit_leaf_wins() {
    let mut merged = json!({ "tree": "from-config.nwk", "seed": 1 });
    let cli = json!({ "tree": "from-cli.nwk", "seed": 999 });
    apply_cli_overrides(&mut merged, &cli, &ids(&["tree"]));
    assert_eq!(json!({ "tree": "from-cli.nwk", "seed": 1 }), merged);
  }

  #[test]
  fn test_config_apply_cli_overrides_recurses_into_flatten_container() {
    let mut merged = json!({ "alignment": { "aln": "from-config.fasta", "vcf": null } });
    let cli = json!({ "alignment": { "aln": "from-cli.fasta", "vcf": null } });
    apply_cli_overrides(&mut merged, &cli, &ids(&["aln"]));
    assert_eq!(json!({ "alignment": { "aln": "from-cli.fasta", "vcf": null } }), merged);
  }

  #[test]
  fn test_config_apply_cli_overrides_no_explicit_keeps_config() {
    let config = json!({ "alignment": { "aln": "keep.fasta" }, "seed": 7 });
    let mut merged = config.clone();
    let cli = json!({ "alignment": { "aln": "ignored.fasta" }, "seed": 42 });
    apply_cli_overrides(&mut merged, &cli, &ids(&[]));
    assert_eq!(config, merged);
  }

  #[test]
  fn test_config_apply_cli_overrides_replaces_array_leaf_wholesale() {
    let mut merged = json!({ "cdses": ["a", "b", "c"] });
    let cli = json!({ "cdses": ["x"] });
    apply_cli_overrides(&mut merged, &cli, &ids(&["cdses"]));
    assert_eq!(json!({ "cdses": ["x"] }), merged);
  }

  mod end_to_end {
    use crate::cli::config::overlay_config;
    use crate::cli::treetime_cli::TreetimeArgs;
    use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
    use crate::commands::clock::args::TreetimeClockArgsRaw;
    use crate::commands::homoplasy::args::TreetimeHomoplasyArgsRaw;
    use crate::commands::mugration::args::TreetimeMugrationArgsRaw;
    use crate::commands::optimize::args::TreetimeOptimizeArgsRaw;
    use crate::commands::prune::args::TreetimePruneArgsRaw;
    use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
    use clap::{CommandFactory, FromArgMatches};
    use eyre::Report;
    use indoc::indoc;
    use pretty_assertions::assert_eq;
    use std::fs;
    use std::path::{Path, PathBuf};
    use tempfile::tempdir;
    use treetime_utils::{assert_error, pretty_assert_ulps_eq};

    fn parse_timetree(argv: &[&str]) -> TreetimeTimetreeArgsRaw {
      let matches = TreetimeArgs::command().get_matches_from(argv);
      let sub = matches.subcommand_matches("timetree").unwrap();
      let mut args = TreetimeTimetreeArgsRaw::from_arg_matches(sub).unwrap();
      overlay_config(&mut args, sub).unwrap();
      args
    }

    fn write_config(dir: &Path) -> PathBuf {
      let path = dir.join("timetree.yaml");
      fs::write(
        &path,
        indoc! {r"
          skyline_n_points: 4
          coalescent_confidence: 3.5
        "},
      )
      .unwrap();
      path
    }

    #[test]
    fn test_config_timetree_file_overrides_defaults() {
      let dir = tempdir().unwrap();
      let path = write_config(dir.path());
      let args = parse_timetree(&["treetime", "timetree", "--config", path.to_str().unwrap()]);
      assert_eq!(4, args.skyline_n_points);
      pretty_assert_ulps_eq!(3.5, args.coalescent_confidence);
      pretty_assert_ulps_eq!(50.0, args.gen_per_year);
    }

    #[test]
    fn test_config_timetree_cli_flag_overrides_file() {
      let dir = tempdir().unwrap();
      let path = write_config(dir.path());
      let args = parse_timetree(&[
        "treetime",
        "timetree",
        "--config",
        path.to_str().unwrap(),
        "--skyline-n-points",
        "8",
      ]);
      assert_eq!(8, args.skyline_n_points);
      pretty_assert_ulps_eq!(3.5, args.coalescent_confidence);
    }

    #[test]
    fn test_config_timetree_absent_keeps_defaults() {
      let args = parse_timetree(&["treetime", "timetree"]);
      assert_eq!(20, args.skyline_n_points);
      pretty_assert_ulps_eq!(2.0, args.coalescent_confidence);
    }

    fn resolve_ancestral(argv: &[&str]) -> Result<TreetimeAncestralArgs, Report> {
      let matches = TreetimeArgs::command().get_matches_from(argv);
      let sub = matches.subcommand_matches("ancestral").unwrap();
      let mut args = TreetimeAncestralArgsRaw::from_arg_matches(sub).unwrap();
      overlay_config(&mut args, sub)?;
      TreetimeAncestralArgs::try_from(args)
    }

    #[test]
    fn test_config_ancestral_config_satisfies_required_tree() {
      let dir = tempdir().unwrap();
      let path = dir.path().join("ancestral.yaml");
      fs::write(
        &path,
        indoc! {r"
        tree: from-config.nwk
      "},
      )
      .unwrap();
      let args = resolve_ancestral(&["treetime", "ancestral", "--config", path.to_str().unwrap()]).unwrap();
      assert_eq!(Path::new("from-config.nwk"), args.tree());
    }

    #[test]
    fn test_config_ancestral_accepts_schema_key() {
      let dir = tempdir().unwrap();
      let path = dir.path().join("ancestral.yaml");
      fs::write(
        &path,
        indoc! {r#"
          "$schema": "https://raw.githubusercontent.com/neherlab/treetime/rust/packages/schemas/input-config-ancestral.schema.json"
          tree: from-config.nwk
        "#},
      )
      .unwrap();
      let args = resolve_ancestral(&["treetime", "ancestral", "--config", path.to_str().unwrap()]).unwrap();
      assert_eq!(Path::new("from-config.nwk"), args.tree());
    }

    #[test]
    fn test_config_ancestral_missing_tree_everywhere_errors() {
      let result = resolve_ancestral(&["treetime", "ancestral"]);
      assert_error!(
        result,
        "the following required arguments were not provided:\n  --tree <TREE>"
      );
    }

    mod reject_unknown_key {
      use super::*;

      #[test]
      fn test_config_ancestral_rejects_unknown_key() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("config.yaml");
        fs::write(&path, "definitely_not_a_real_field: 1\n").unwrap();
        let matches =
          TreetimeArgs::command().get_matches_from(["treetime", "ancestral", "--config", path.to_str().unwrap()]);
        let sub = matches.subcommand_matches("ancestral").unwrap();
        let mut args = TreetimeAncestralArgsRaw::from_arg_matches(sub).unwrap();
        let result = overlay_config(&mut args, sub);
        assert_error!(
          result,
          "invalid configuration: unknown field `definitely_not_a_real_field`"
        );
      }

      #[test]
      fn test_config_clock_rejects_unknown_key() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("config.yaml");
        fs::write(&path, "definitely_not_a_real_field: 1\n").unwrap();
        let matches =
          TreetimeArgs::command().get_matches_from(["treetime", "clock", "--config", path.to_str().unwrap()]);
        let sub = matches.subcommand_matches("clock").unwrap();
        let mut args = TreetimeClockArgsRaw::from_arg_matches(sub).unwrap();
        let result = overlay_config(&mut args, sub);
        assert_error!(
          result,
          "invalid configuration: unknown field `definitely_not_a_real_field`"
        );
      }

      #[test]
      fn test_config_timetree_rejects_unknown_key() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("config.yaml");
        fs::write(&path, "definitely_not_a_real_field: 1\n").unwrap();
        let matches =
          TreetimeArgs::command().get_matches_from(["treetime", "timetree", "--config", path.to_str().unwrap()]);
        let sub = matches.subcommand_matches("timetree").unwrap();
        let mut args = TreetimeTimetreeArgsRaw::from_arg_matches(sub).unwrap();
        let result = overlay_config(&mut args, sub);
        assert_error!(
          result,
          "invalid configuration: unknown field `definitely_not_a_real_field`"
        );
      }

      #[test]
      fn test_config_optimize_rejects_unknown_key() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("config.yaml");
        fs::write(&path, "definitely_not_a_real_field: 1\n").unwrap();
        let matches =
          TreetimeArgs::command().get_matches_from(["treetime", "optimize", "--config", path.to_str().unwrap()]);
        let sub = matches.subcommand_matches("optimize").unwrap();
        let mut args = TreetimeOptimizeArgsRaw::from_arg_matches(sub).unwrap();
        let result = overlay_config(&mut args, sub);
        assert_error!(
          result,
          "invalid configuration: unknown field `definitely_not_a_real_field`"
        );
      }

      #[test]
      fn test_config_prune_rejects_unknown_key() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("config.yaml");
        fs::write(&path, "definitely_not_a_real_field: 1\n").unwrap();
        let matches =
          TreetimeArgs::command().get_matches_from(["treetime", "prune", "--config", path.to_str().unwrap()]);
        let sub = matches.subcommand_matches("prune").unwrap();
        let mut args = TreetimePruneArgsRaw::from_arg_matches(sub).unwrap();
        let result = overlay_config(&mut args, sub);
        assert_error!(
          result,
          "invalid configuration: unknown field `definitely_not_a_real_field`"
        );
      }

      #[test]
      fn test_config_mugration_rejects_unknown_key() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("config.yaml");
        fs::write(&path, "definitely_not_a_real_field: 1\n").unwrap();
        let matches =
          TreetimeArgs::command().get_matches_from(["treetime", "mugration", "--config", path.to_str().unwrap()]);
        let sub = matches.subcommand_matches("mugration").unwrap();
        let mut args = TreetimeMugrationArgsRaw::from_arg_matches(sub).unwrap();
        let result = overlay_config(&mut args, sub);
        assert_error!(
          result,
          "invalid configuration: unknown field `definitely_not_a_real_field`"
        );
      }

      #[test]
      fn test_config_homoplasy_rejects_unknown_key() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("config.yaml");
        fs::write(&path, "definitely_not_a_real_field: 1\n").unwrap();
        let matches =
          TreetimeArgs::command().get_matches_from(["treetime", "homoplasy", "--config", path.to_str().unwrap()]);
        let sub = matches.subcommand_matches("homoplasy").unwrap();
        let mut args = TreetimeHomoplasyArgsRaw::from_arg_matches(sub).unwrap();
        let result = overlay_config(&mut args, sub);
        assert_error!(
          result,
          "invalid configuration: unknown field `definitely_not_a_real_field`"
        );
      }

      #[test]
      fn test_config_ancestral_rejects_unknown_nested_key() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("config.yaml");
        fs::write(
          &path,
          indoc! {r"
          model_args:
            not_a_model_field: 1
        "},
        )
        .unwrap();
        let matches =
          TreetimeArgs::command().get_matches_from(["treetime", "ancestral", "--config", path.to_str().unwrap()]);
        let sub = matches.subcommand_matches("ancestral").unwrap();
        let mut args = TreetimeAncestralArgsRaw::from_arg_matches(sub).unwrap();
        let result = overlay_config(&mut args, sub);
        assert_error!(result, "invalid configuration: unknown field `model_args`");
      }
    }

    #[test]
    fn test_config_timetree_reads_json_document() {
      let dir = tempdir().unwrap();
      let path = dir.path().join("config.json");
      fs::write(&path, r#"{ "skyline_n_points": 4 }"#).unwrap();
      let args = parse_timetree(&["treetime", "timetree", "--config", path.to_str().unwrap()]);
      assert_eq!(4, args.skyline_n_points);
    }

    #[test]
    fn test_config_timetree_preserves_scientific_notation_precision() {
      let dir = tempdir().unwrap();
      let path = dir.path().join("config.yaml");
      fs::write(&path, "clock_rate: 5.7e-05\n").unwrap();
      let args = parse_timetree(&["treetime", "timetree", "--config", path.to_str().unwrap()]);
      pretty_assert_ulps_eq!(5.7e-05, args.clock_rate.unwrap());
    }
  }

  mod required_args {
    use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
    use crate::commands::clock::args::{TreetimeClockArgs, TreetimeClockArgsRaw};
    use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, TreetimeHomoplasyArgsRaw};
    use crate::commands::mugration::args::{TreetimeMugrationArgs, TreetimeMugrationArgsRaw};
    use crate::commands::optimize::args::{TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
    use crate::commands::prune::args::{TreetimePruneArgs, TreetimePruneArgsRaw};
    use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
    use pretty_assertions::assert_eq;
    use treetime_utils::assert_error;

    #[test]
    fn test_config_required_ancestral_missing_tree_errors() {
      assert_error!(
        TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw::default()),
        "the following required arguments were not provided:\n  --tree <TREE>"
      );
    }

    #[test]
    fn test_config_required_ancestral_with_tree_ok() {
      let raw = TreetimeAncestralArgsRaw {
        tree: Some("tree.nwk".into()),
        ..Default::default()
      };
      let args = TreetimeAncestralArgs::try_from(raw).unwrap();
      assert_eq!(std::path::Path::new("tree.nwk"), args.tree());
    }

    #[test]
    fn test_config_required_optimize_missing_tree_errors() {
      assert_error!(
        TreetimeOptimizeArgs::try_from(TreetimeOptimizeArgsRaw::default()),
        "the following required arguments were not provided:\n  --tree <TREE>"
      );
    }

    #[test]
    fn test_config_required_prune_missing_tree_errors() {
      assert_error!(
        TreetimePruneArgs::try_from(TreetimePruneArgsRaw::default()),
        "the following required arguments were not provided:\n  --tree <TREE>"
      );
    }

    #[test]
    fn test_config_required_clock_missing_metadata_errors() {
      assert_error!(
        TreetimeClockArgs::try_from(TreetimeClockArgsRaw::default()),
        "the following required arguments were not provided:\n  --metadata <METADATA>"
      );
    }

    #[test]
    fn test_config_required_mugration_missing_both_lists_both() {
      assert_error!(
        TreetimeMugrationArgs::try_from(TreetimeMugrationArgsRaw::default()),
        "the following required arguments were not provided:\n  --metadata <METADATA>\n  --attribute <ATTRIBUTE>"
      );
    }

    #[test]
    fn test_config_required_mugration_only_attribute_missing_lists_attribute() {
      let raw = TreetimeMugrationArgsRaw {
        metadata: Some("metadata.tsv".into()),
        ..Default::default()
      };
      assert_error!(
        TreetimeMugrationArgs::try_from(raw),
        "the following required arguments were not provided:\n  --attribute <ATTRIBUTE>"
      );
    }

    #[test]
    fn test_config_required_mugration_only_metadata_missing_lists_metadata() {
      let raw = TreetimeMugrationArgsRaw {
        attribute: Some("country".to_owned()),
        ..Default::default()
      };
      assert_error!(
        TreetimeMugrationArgs::try_from(raw),
        "the following required arguments were not provided:\n  --metadata <METADATA>"
      );
    }

    #[test]
    fn test_config_required_homoplasy_missing_tree_errors() {
      assert_error!(
        TreetimeHomoplasyArgs::try_from(TreetimeHomoplasyArgsRaw::default()),
        "the following required arguments were not provided:\n  --tree <TREE>"
      );
    }

    #[test]
    fn test_config_required_timetree_defaults_ok() {
      TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw::default()).unwrap();
    }
  }

  mod round_trip {
    use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
    use crate::commands::clock::args::TreetimeClockArgsRaw;
    use crate::commands::homoplasy::args::TreetimeHomoplasyArgsRaw;
    use crate::commands::mugration::args::TreetimeMugrationArgsRaw;
    use crate::commands::optimize::args::TreetimeOptimizeArgsRaw;
    use crate::commands::prune::args::TreetimePruneArgsRaw;
    use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
    use pretty_assertions::assert_eq;
    use serde_json::{Value, to_value};

    #[test]
    fn test_config_round_trip_ancestral() {
      let value: Value = to_value(TreetimeAncestralArgsRaw::default()).unwrap();
      let back: TreetimeAncestralArgsRaw = serde::Deserialize::deserialize(&value).unwrap();
      assert_eq!(value, to_value(back).unwrap());
    }

    #[test]
    fn test_config_round_trip_clock() {
      let value: Value = to_value(TreetimeClockArgsRaw::default()).unwrap();
      let back: TreetimeClockArgsRaw = serde::Deserialize::deserialize(&value).unwrap();
      assert_eq!(value, to_value(back).unwrap());
    }

    #[test]
    fn test_config_round_trip_homoplasy() {
      let value: Value = to_value(TreetimeHomoplasyArgsRaw::default()).unwrap();
      let back: TreetimeHomoplasyArgsRaw = serde::Deserialize::deserialize(&value).unwrap();
      assert_eq!(value, to_value(back).unwrap());
    }

    #[test]
    fn test_config_round_trip_mugration() {
      let value: Value = to_value(TreetimeMugrationArgsRaw::default()).unwrap();
      let back: TreetimeMugrationArgsRaw = serde::Deserialize::deserialize(&value).unwrap();
      assert_eq!(value, to_value(back).unwrap());
    }

    #[test]
    fn test_config_round_trip_optimize() {
      let value: Value = to_value(TreetimeOptimizeArgsRaw::default()).unwrap();
      let back: TreetimeOptimizeArgsRaw = serde::Deserialize::deserialize(&value).unwrap();
      assert_eq!(value, to_value(back).unwrap());
    }

    #[test]
    fn test_config_round_trip_prune() {
      let value: Value = to_value(TreetimePruneArgsRaw::default()).unwrap();
      let back: TreetimePruneArgsRaw = serde::Deserialize::deserialize(&value).unwrap();
      assert_eq!(value, to_value(back).unwrap());
    }

    #[test]
    fn test_config_round_trip_timetree() {
      let value: Value = to_value(TreetimeTimetreeArgsRaw::default()).unwrap();
      let back: TreetimeTimetreeArgsRaw = serde::Deserialize::deserialize(&value).unwrap();
      assert_eq!(value, to_value(back).unwrap());
    }
  }
}
