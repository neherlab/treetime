#[cfg(test)]
mod tests {
  use crate::confine::PathPolicy;
  use app_commands::command::AppCommand;
  use helpers::{Dirs, dirs};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::fs;
  use std::os::unix::fs::symlink;
  use treetime_utils::assert_error;

  #[test]
  fn test_confine_relative_input_resolves_inside_data_dir() {
    let Dirs { data, outside, out } = dirs();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({ "tree": "sub/tree.nwk", "alignment": ["sub/aln.fasta"] });
    policy.confine(AppCommand::Ancestral, &mut config, out.path()).unwrap();
    let data = data.path().canonicalize().unwrap();
    assert_eq!(
      json!({
        "tree": data.join("sub/tree.nwk"),
        "alignment": [data.join("sub/aln.fasta")],
        "output_all": out.path(),
      }),
      config
    );
    drop(outside);
  }

  #[test]
  fn test_confine_absolute_input_inside_data_dir_is_accepted() {
    let Dirs { data, out, .. } = dirs();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let tree = data.path().canonicalize().unwrap().join("sub/tree.nwk");
    let mut config = json!({ "tree": tree });
    policy.confine(AppCommand::Optimize, &mut config, out.path()).unwrap();
    assert_eq!(json!(tree), config["tree"]);
  }

  #[test]
  fn test_confine_rejects_parent_directory_escape() {
    let Dirs { data, out, outside } = dirs();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let escape = format!(
      "../{}/secret.nwk",
      outside.path().file_name().unwrap().to_string_lossy()
    );
    let mut config = json!({ "tree": escape });
    assert_error!(
      policy.confine(AppCommand::Optimize, &mut config, out.path()),
      format!("input `{escape}` of setting `tree` is outside the directories the server reads inputs from")
    );
  }

  #[test]
  fn test_confine_rejects_absolute_path_outside_roots() {
    let Dirs { data, out, outside } = dirs();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let secret = outside.path().join("secret.nwk");
    let mut config = json!({ "metadata": secret, "tree": "sub/tree.nwk" });
    assert_error!(
      policy.confine(AppCommand::Clock, &mut config, out.path()),
      format!(
        "input `{}` of setting `metadata` is outside the directories the server reads inputs from",
        secret.display()
      )
    );
  }

  #[test]
  fn test_confine_rejects_symlink_that_resolves_outside() {
    let Dirs { data, out, outside } = dirs();
    symlink(outside.path().join("secret.nwk"), data.path().join("link.nwk")).unwrap();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({ "tree": "link.nwk" });
    assert_error!(
      policy.confine(AppCommand::Prune, &mut config, out.path()),
      "input `link.nwk` of setting `tree` is outside the directories the server reads inputs from"
    );
  }

  #[test]
  fn test_confine_rejects_symlinked_directory_that_resolves_outside() {
    let Dirs { data, out, outside } = dirs();
    symlink(outside.path(), data.path().join("linked")).unwrap();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({ "alignment": ["sub/aln.fasta", "linked/secret.nwk"] });
    assert_error!(
      policy.confine(AppCommand::Ancestral, &mut config, out.path()),
      "input `linked/secret.nwk` of setting `alignment` is outside the directories the server reads inputs from"
    );
  }

  #[test]
  fn test_confine_accepts_symlink_that_stays_inside() {
    let Dirs { data, out, .. } = dirs();
    symlink(data.path().join("sub/tree.nwk"), data.path().join("alias.nwk")).unwrap();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({ "tree": "alias.nwk" });
    policy.confine(AppCommand::Prune, &mut config, out.path()).unwrap();
    assert_eq!(
      json!(data.path().canonicalize().unwrap().join("sub/tree.nwk")),
      config["tree"]
    );
  }

  #[test]
  fn test_confine_accepts_extra_input_root() {
    let Dirs { data, out, outside } = dirs();
    let policy = PathPolicy::new(data.path(), &[outside.path().to_path_buf()]).unwrap();
    let secret = outside.path().join("secret.nwk");
    let mut config = json!({ "tree": secret });
    policy.confine(AppCommand::Prune, &mut config, out.path()).unwrap();
    assert_eq!(json!(secret.canonicalize().unwrap()), config["tree"]);
  }

  #[test]
  fn test_confine_missing_input_is_an_error() {
    let Dirs { data, out, .. } = dirs();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({ "tree": "missing.nwk" });
    assert_error!(
      policy.confine(AppCommand::Prune, &mut config, out.path()),
      "input `missing.nwk` of setting `tree` cannot be read: No such file or directory (os error 2)"
    );
  }

  #[test]
  fn test_confine_discards_client_output_paths() {
    let Dirs { data, out, .. } = dirs();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({
      "tree": "sub/tree.nwk",
      "output_all": "/etc",
      "output_tree_nwk": "/etc/passwd",
      "output_clock_model": "../x.json",
      "plot_rtt": "/tmp/plot.png",
      "output_selection": ["nwk"],
      "output_nwk_style": ["beast"],
    });
    policy.confine(AppCommand::Clock, &mut config, out.path()).unwrap();
    assert_eq!(
      json!({
        "tree": data.path().canonicalize().unwrap().join("sub/tree.nwk"),
        "output_all": out.path(),
        "output_selection": ["nwk"],
        "output_nwk_style": ["beast"],
      }),
      config
    );
  }

  #[test]
  fn test_confine_translation_template_expands_listed_cdses() {
    let Dirs { data, out, outside } = dirs();
    fs::write(data.path().join("sub/translation_E.fasta"), ">a\nM\n").unwrap();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({
      "tree": "sub/tree.nwk",
      "translations": "sub/translation_{cds}.fasta",
      "cdses": ["E"],
    });
    policy.confine(AppCommand::Ancestral, &mut config, out.path()).unwrap();

    let escape = format!("../../{}/secret", outside.path().file_name().unwrap().to_string_lossy());
    let mut config = json!({
      "tree": "sub/tree.nwk",
      "translations": "sub/{cds}.nwk",
      "cdses": [escape],
    });
    assert_error!(
      policy.confine(AppCommand::Ancestral, &mut config, out.path()),
      format!(
        "input `{}/sub/{escape}.nwk` of setting `translations` is outside the directories the server reads inputs from",
        data.path().canonicalize().unwrap().display()
      )
    );
  }

  #[test]
  fn test_confine_rejects_missing_data_dir() {
    let Dirs { data, .. } = dirs();
    let missing = data.path().join("missing");
    assert_error!(
      PathPolicy::new(&missing, &[]),
      format!(
        "When resolving the data directory: When resolving directory '{}': No such file or directory (os error 2)",
        missing.display()
      )
    );
  }

  mod helpers {
    use std::fs;
    use tempfile::{TempDir, tempdir};

    pub(super) struct Dirs {
      pub data: TempDir,
      pub out: TempDir,
      pub outside: TempDir,
    }

    pub(super) fn dirs() -> Dirs {
      let data = tempdir().unwrap();
      fs::create_dir_all(data.path().join("sub")).unwrap();
      fs::write(data.path().join("sub/tree.nwk"), "(A:1,B:1);\n").unwrap();
      fs::write(data.path().join("sub/aln.fasta"), ">A\nACGT\n>B\nACGA\n").unwrap();
      let outside = tempdir().unwrap();
      fs::write(outside.path().join("secret.nwk"), "(S:1,T:1);\n").unwrap();
      Dirs {
        data,
        out: tempdir().unwrap(),
        outside,
      }
    }
  }
}
