use crate::cli::diagnostics::entry::check_command_config;
use crate::cli::diagnostics::source::{ConfigSource, RawDiagnostic, render_and_bail};
use crate::cli::schema::command_schema;
use clap::ArgMatches;
use clap::parser::ValueSource;
use eyre::Report;
use schemars::JsonSchema;
use serde::Serialize;
use serde::de::DeserializeOwned;
use serde_json::Value;
use std::collections::HashSet;
use std::path::PathBuf;
use treetime_utils::io::fs::read_file_to_string;

/// Overlay a `--config` file onto CLI-parsed command args when `--config` was given.
///
/// Precedence, highest to lowest: an explicit command-line flag, then the config file, then the
/// built-in default. The file holds the command's configuration object (JSON or YAML). It is layered
/// over the serialized defaults, then explicit command-line leaves are layered on top, producing the
/// effective configuration. That effective object is validated against the command's schema through
/// the shared diagnostics layer, so a bad value renders a caret-annotated error pointing into the
/// file (a leaf supplied only on the command line has no file span, so it falls back to a plain
/// message). The typed value is materialized only after validation passes.
pub fn overlay_config<T>(args: &mut T, matches: &ArgMatches) -> Result<(), Report>
where
  T: Serialize + DeserializeOwned + Default + JsonSchema,
{
  let Some(config_path) = matches.get_one::<PathBuf>("config") else {
    return Ok(());
  };

  // Ids the user supplied on the command line, as opposed to defaults or the config file. Only these
  // override the file.
  let explicit: HashSet<String> = matches
    .ids()
    .filter(|id| matches.value_source(id.as_str()) == Some(ValueSource::CommandLine))
    .map(|id| id.to_string())
    .collect();

  let text = read_file_to_string(config_path)?;
  let source = ConfigSource::new(config_path.display().to_string(), text.clone());
  let file_value: Value = match serde_yaml::from_str(&text) {
    Ok(file_value) => file_value,
    Err(err) => {
      render_and_bail(
        &source,
        "invalid configuration",
        vec![RawDiagnostic::new(
          "config::syntax",
          format!("could not parse config: {err}"),
        )],
      )?;
      unreachable!("render_and_bail returns an error whenever diagnostics are present");
    },
  };

  let mut merged = serde_json::to_value(T::default())?;
  merge_value(&mut merged, &file_value);
  let cli = serde_json::to_value(&*args)?;
  apply_cli_overrides(&mut merged, &cli, &explicit);

  check_command_config(&source, &merged, &command_schema::<T>())?;

  *args = serde_json::from_value(merged)?;
  Ok(())
}

/// Layer `overlay` onto `base`, recursing through nested objects so a partial config overrides only
/// the leaves it sets, leaving sibling defaults intact.
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

/// Copy every explicitly-set command-line leaf from `cli` into `merged`, recursing through the
/// nested objects that `clap(flatten)` produces so flattened arg groups are handled uniformly.
///
/// Arg ids are unique within a command, so a key in `explicit` always denotes a leaf the user set on
/// the command line; such leaves replace the config value wholesale. Objects that are not themselves
/// args are flatten containers, descended to reach the explicit leaves inside them.
fn apply_cli_overrides(merged: &mut Value, cli: &Value, explicit: &HashSet<String>) {
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

  fn ids(names: &[&str]) -> HashSet<String> {
    names.iter().map(|s| (*s).to_owned()).collect()
  }

  // An explicitly-set leaf overrides the config value; unset leaves keep the config value.
  #[test]
  fn test_config_apply_cli_overrides_explicit_leaf_wins() {
    let mut merged = json!({ "tree": "from-config.nwk", "seed": 1 });
    let cli = json!({ "tree": "from-cli.nwk", "seed": 999 });
    apply_cli_overrides(&mut merged, &cli, &ids(&["tree"]));
    assert_eq!(json!({ "tree": "from-cli.nwk", "seed": 1 }), merged);
  }

  // Explicit leaves inside a `clap(flatten)` container are reached by recursion; the container name
  // itself is not an arg id.
  #[test]
  fn test_config_apply_cli_overrides_recurses_into_flatten_container() {
    let mut merged = json!({ "alignment": { "aln": "from-config.fasta", "vcf": null } });
    let cli = json!({ "alignment": { "aln": "from-cli.fasta", "vcf": null } });
    apply_cli_overrides(&mut merged, &cli, &ids(&["aln"]));
    assert_eq!(json!({ "alignment": { "aln": "from-cli.fasta", "vcf": null } }), merged);
  }

  // With nothing explicit, the config object is left untouched.
  #[test]
  fn test_config_apply_cli_overrides_no_explicit_keeps_config() {
    let config = json!({ "alignment": { "aln": "keep.fasta" }, "seed": 7 });
    let mut merged = config.clone();
    let cli = json!({ "alignment": { "aln": "ignored.fasta" }, "seed": 42 });
    apply_cli_overrides(&mut merged, &cli, &ids(&[]));
    assert_eq!(config, merged);
  }

  // A whole-value leaf (array from `value_delimiter`) is replaced wholesale, not merged elementwise.
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
    use crate::cli::validate::Validate;
    use clap::{CommandFactory, FromArgMatches};
    use eyre::Report;
    use indoc::indoc;
    use pretty_assertions::assert_eq;
    use std::fs;
    use std::path::{Path, PathBuf};
    use tempfile::tempdir;
    use treetime::commands::ancestral::args::TreetimeAncestralArgs;
    use treetime::commands::timetree::args::TreetimeTimetreeArgs;
    use treetime_utils::{assert_error, pretty_assert_ulps_eq};

    // Drive the real parse path: full clap parse (which records value sources), then the `--config`
    // overlay, exactly as `treetime_parse_cli_args` does for a subcommand.
    fn parse_timetree(argv: &[&str]) -> TreetimeTimetreeArgs {
      let matches = TreetimeArgs::command().get_matches_from(argv);
      let sub = matches.subcommand_matches("timetree").unwrap();
      let mut args = TreetimeTimetreeArgs::from_arg_matches(sub).unwrap();
      overlay_config(&mut args, sub).unwrap();
      args
    }

    fn write_config(dir: &Path) -> PathBuf {
      // A partial YAML config: only two fields set, everything else must fall back to defaults.
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

    // Without overriding flags, config values win over defaults and unspecified fields keep defaults.
    #[test]
    fn test_config_timetree_file_overrides_defaults() {
      let dir = tempdir().unwrap();
      let path = write_config(dir.path());
      let args = parse_timetree(&["treetime", "timetree", "--config", path.to_str().unwrap()]);
      assert_eq!(4, args.skyline_n_points); // from config
      pretty_assert_ulps_eq!(3.5, args.coalescent_confidence); // from config
      pretty_assert_ulps_eq!(50.0, args.gen_per_year); // default, untouched by config
    }

    // An explicit command-line flag overrides the config file; other config values remain in effect.
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
      assert_eq!(8, args.skyline_n_points); // explicit CLI beats config
      pretty_assert_ulps_eq!(3.5, args.coalescent_confidence); // still from config
    }

    // Without `--config`, parsing is unaffected: defaults stand.
    #[test]
    fn test_config_timetree_absent_keeps_defaults() {
      let args = parse_timetree(&["treetime", "timetree"]);
      assert_eq!(20, args.skyline_n_points);
      pretty_assert_ulps_eq!(2.0, args.coalescent_confidence);
    }

    // Full resolve path for a command with a required argument: parse, overlay `--config`, validate.
    fn resolve_ancestral(argv: &[&str]) -> Result<TreetimeAncestralArgs, Report> {
      let matches = TreetimeArgs::command().get_matches_from(argv);
      let sub = matches.subcommand_matches("ancestral").unwrap();
      let mut args = TreetimeAncestralArgs::from_arg_matches(sub).unwrap();
      overlay_config(&mut args, sub)?;
      args.validate()?;
      Ok(args)
    }

    // The headline of Part 1: a required argument may be supplied entirely by the config file, with no
    // corresponding command-line flag, and validation passes.
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

    // When a required argument is present in neither the command line nor the config file, validation
    // errors with the clap-style message.
    #[test]
    fn test_config_ancestral_missing_tree_everywhere_errors() {
      let result = resolve_ancestral(&["treetime", "ancestral"]);
      assert_error!(
        result,
        "the following required arguments were not provided:\n  --tree <TREE>"
      );
    }
  }
}
