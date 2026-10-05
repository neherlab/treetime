#[cfg(test)]
mod tests {
  use crate::confine::PathPolicy;
  use app_commands::command::AppCommand;
  use helpers::{Dirs, dirs};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::env::current_dir;
  use std::fs;
  use std::os::unix::fs::symlink;
  use std::path::Path;
  use treetime_utils::assert_error;

  #[test]
  fn test_confine_relative_input_resolves_inside_examples_dir() {
    let Dirs { data, outside } = dirs();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({ "tree": "sub/tree.nwk", "alignment": ["sub/aln.fasta"] });
    policy.confine(AppCommand::Ancestral, &mut config).unwrap();
    let data = data.path().canonicalize().unwrap();
    assert_eq!(
      json!({
        "tree": data.join("sub/tree.nwk"),
        "alignment": [data.join("sub/aln.fasta")],
      }),
      config
    );
    drop(outside);
  }

  #[test]
  fn test_confine_absolute_input_inside_examples_dir_is_accepted() {
    let Dirs { data, .. } = dirs();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let tree = data.path().canonicalize().unwrap().join("sub/tree.nwk");
    let mut config = json!({ "tree": tree });
    policy.confine(AppCommand::Optimize, &mut config).unwrap();
    assert_eq!(json!(tree), config["tree"]);
  }

  #[test]
  fn test_confine_rejects_parent_directory_escape() {
    let Dirs { data, outside } = dirs();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let escape = format!(
      "../{}/secret.nwk",
      outside.path().file_name().unwrap().to_string_lossy()
    );
    let mut config = json!({ "tree": escape });
    assert_error!(
      policy.confine(AppCommand::Optimize, &mut config),
      format!("input `{escape}` of setting `tree` is outside the directories the server reads inputs from")
    );
  }

  #[test]
  fn test_confine_rejects_absolute_path_outside_roots() {
    let Dirs { data, outside } = dirs();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let secret = outside.path().join("secret.nwk");
    let mut config = json!({ "metadata": secret, "tree": "sub/tree.nwk" });
    assert_error!(
      policy.confine(AppCommand::Clock, &mut config),
      format!(
        "input `{}` of setting `metadata` is outside the directories the server reads inputs from",
        secret.display()
      )
    );
  }

  #[test]
  fn test_confine_rejects_symlink_that_resolves_outside() {
    let Dirs { data, outside } = dirs();
    symlink(outside.path().join("secret.nwk"), data.path().join("link.nwk")).unwrap();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({ "tree": "link.nwk" });
    assert_error!(
      policy.confine(AppCommand::Prune, &mut config),
      "input `link.nwk` of setting `tree` is outside the directories the server reads inputs from"
    );
  }

  #[test]
  fn test_confine_rejects_symlinked_directory_that_resolves_outside() {
    let Dirs { data, outside } = dirs();
    symlink(outside.path(), data.path().join("linked")).unwrap();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({ "alignment": ["sub/aln.fasta", "linked/secret.nwk"] });
    assert_error!(
      policy.confine(AppCommand::Ancestral, &mut config),
      "input `linked/secret.nwk` of setting `alignment` is outside the directories the server reads inputs from"
    );
  }

  #[test]
  fn test_confine_accepts_symlink_that_stays_inside() {
    let Dirs { data, .. } = dirs();
    symlink(data.path().join("sub/tree.nwk"), data.path().join("alias.nwk")).unwrap();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({ "tree": "alias.nwk" });
    policy.confine(AppCommand::Prune, &mut config).unwrap();
    assert_eq!(
      json!(data.path().canonicalize().unwrap().join("sub/tree.nwk")),
      config["tree"]
    );
  }

  #[test]
  fn test_confine_accepts_extra_input_root() {
    let Dirs { data, outside } = dirs();
    let policy = PathPolicy::new(data.path(), &[outside.path().to_path_buf()]).unwrap();
    let secret = outside.path().join("secret.nwk");
    let mut config = json!({ "tree": secret });
    policy.confine(AppCommand::Prune, &mut config).unwrap();
    assert_eq!(json!(secret.canonicalize().unwrap()), config["tree"]);
  }

  #[test]
  fn test_confine_relative_input_missing_from_examples_dir_resolves_from_working_directory() {
    let policy = PathPolicy::new(Path::new("../../data"), &[]).unwrap();
    let mut config = json!({ "tree": "../../data/zika/20/tree.nwk" });
    policy.confine(AppCommand::Prune, &mut config).unwrap();
    assert_eq!(
      json!(Path::new("../../data/zika/20/tree.nwk").canonicalize().unwrap()),
      config["tree"]
    );
  }

  #[test]
  fn test_confine_working_directory_input_outside_examples_dir_is_rejected() {
    let policy = PathPolicy::new(Path::new("../../data"), &[]).unwrap();
    let mut config = json!({ "tree": "Cargo.toml" });
    assert_error!(
      policy.confine(AppCommand::Prune, &mut config),
      "input `Cargo.toml` of setting `tree` is outside the directories the server reads inputs from"
    );
  }

  #[test]
  fn test_confine_missing_input_is_an_error() {
    let Dirs { data, .. } = dirs();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({ "tree": "missing.nwk" });
    assert_error!(
      policy.confine(AppCommand::Prune, &mut config),
      "input `missing.nwk` of setting `tree` cannot be read: No such file or directory (os error 2)"
    );
  }

  #[test]
  fn test_confine_translation_template_expands_listed_cdses() {
    let Dirs { data, outside } = dirs();
    fs::write(data.path().join("sub/translation_E.fasta"), ">a\nM\n").unwrap();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({
      "tree": "sub/tree.nwk",
      "translations": "sub/translation_{cds}.fasta",
      "cdses": ["E"],
    });
    policy.confine(AppCommand::Ancestral, &mut config).unwrap();

    let escape = format!("../../{}/secret", outside.path().file_name().unwrap().to_string_lossy());
    let mut config = json!({
      "tree": "sub/tree.nwk",
      "translations": "sub/{cds}.nwk",
      "cdses": [escape],
    });
    assert_error!(
      policy.confine(AppCommand::Ancestral, &mut config),
      format!(
        "input `{}/sub/{escape}.nwk` of setting `translations` is outside the directories the server reads inputs from",
        data.path().canonicalize().unwrap().display()
      )
    );
  }

  #[test]
  fn test_confine_translation_template_is_rewritten_to_the_checked_examples_dir_location() {
    let Dirs { data, .. } = dirs();
    fs::create_dir_all(data.path().join("src")).unwrap();
    fs::write(data.path().join("src/lib.rs"), ">a\nM\n").unwrap();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({
      "tree": "sub/tree.nwk",
      "translations": "src/{cds}.rs",
      "cdses": ["lib"],
    });
    policy.confine(AppCommand::Ancestral, &mut config).unwrap();
    assert_eq!(
      json!(data.path().canonicalize().unwrap().join("src/{cds}.rs")),
      config["translations"]
    );
  }

  #[test]
  fn test_confine_translation_template_missing_from_examples_dir_is_rewritten_to_the_working_directory() {
    let Dirs { data, .. } = dirs();
    let policy = PathPolicy::new(data.path(), &[Path::new("src").to_path_buf()]).unwrap();
    let mut config = json!({
      "tree": "sub/tree.nwk",
      "translations": "src/{cds}.rs",
      "cdses": ["lib"],
    });
    policy.confine(AppCommand::Ancestral, &mut config).unwrap();
    assert_eq!(
      json!(current_dir().unwrap().join("src/{cds}.rs")),
      config["translations"]
    );
  }

  #[test]
  fn test_confine_translation_template_in_working_directory_outside_roots_is_rejected() {
    let Dirs { data, .. } = dirs();
    let policy = PathPolicy::new(data.path(), &[]).unwrap();
    let mut config = json!({
      "tree": "sub/tree.nwk",
      "translations": "src/{cds}.rs",
      "cdses": ["lib"],
    });
    assert_error!(
      policy.confine(AppCommand::Ancestral, &mut config),
      format!(
        "input `{}` of setting `translations` is outside the directories the server reads inputs from",
        current_dir().unwrap().join("src/lib.rs").display()
      )
    );
  }

  #[test]
  fn test_confine_rejects_missing_examples_dir() {
    let Dirs { data, .. } = dirs();
    let missing = data.path().join("missing");
    assert_error!(
      PathPolicy::new(&missing, &[]),
      format!(
        "When resolving the examples folder: When resolving directory '{}': No such file or directory (os error 2)",
        missing.display()
      )
    );
  }

  mod helpers {
    use std::fs;
    use tempfile::{TempDir, tempdir};

    pub(super) struct Dirs {
      pub data: TempDir,
      pub outside: TempDir,
    }

    pub(super) fn dirs() -> Dirs {
      let data = tempdir().unwrap();
      fs::create_dir_all(data.path().join("sub")).unwrap();
      fs::write(data.path().join("sub/tree.nwk"), "(A:1,B:1);\n").unwrap();
      fs::write(data.path().join("sub/aln.fasta"), ">A\nACGT\n>B\nACGA\n").unwrap();
      let outside = tempdir().unwrap();
      fs::write(outside.path().join("secret.nwk"), "(S:1,T:1);\n").unwrap();
      Dirs { data, outside }
    }
  }
}
