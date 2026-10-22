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
  use std::slice;
  use treetime_utils::assert_error;

  #[test]
  fn test_confine_absolute_inputs_inside_a_root_are_accepted() {
    let Dirs { data, .. } = dirs();
    let root = data.path().canonicalize().unwrap();
    let policy = PathPolicy::new(slice::from_ref(&root)).unwrap();
    let mut config = json!({ "tree": root.join("sub/tree.nwk"), "alignment": [root.join("sub/aln.fasta")] });
    policy.confine(AppCommand::Ancestral, &mut config).unwrap();
    assert_eq!(
      json!({
        "tree": root.join("sub/tree.nwk"),
        "alignment": [root.join("sub/aln.fasta")],
      }),
      config
    );
  }

  #[test]
  fn test_confine_relative_input_resolves_from_the_working_directory() {
    let policy = PathPolicy::new(&[Path::new("../../data").to_path_buf()]).unwrap();
    let mut config = json!({ "tree": "../../data/zika/20/tree.nwk" });
    policy.confine(AppCommand::Prune, &mut config).unwrap();
    assert_eq!(
      json!(Path::new("../../data/zika/20/tree.nwk").canonicalize().unwrap()),
      config["tree"]
    );
  }

  #[test]
  fn test_confine_working_directory_input_outside_roots_is_rejected() {
    let Dirs { data, .. } = dirs();
    let policy = PathPolicy::new(&[data.path().to_path_buf()]).unwrap();
    let mut config = json!({ "tree": "Cargo.toml" });
    assert_error!(
      policy.confine(AppCommand::Prune, &mut config),
      "input `Cargo.toml` of setting `tree` is outside the directories the server reads inputs from"
    );
  }

  #[test]
  fn test_confine_without_roots_rejects_every_input() {
    let Dirs { data, .. } = dirs();
    let tree = data.path().join("sub/tree.nwk");
    let policy = PathPolicy::new(&[]).unwrap();
    let mut config = json!({ "tree": tree });
    assert_error!(
      policy.confine(AppCommand::Prune, &mut config),
      format!(
        "input `{}` of setting `tree` is outside the directories the server reads inputs from",
        tree.display()
      )
    );
  }

  #[test]
  fn test_confine_rejects_parent_directory_escape() {
    let Dirs { data, outside } = dirs();
    let policy = PathPolicy::new(&[data.path().to_path_buf()]).unwrap();
    let escape = data
      .path()
      .join("sub/../..")
      .join(outside.path().file_name().unwrap())
      .join("secret.nwk");
    let mut config = json!({ "tree": escape });
    assert_error!(
      policy.confine(AppCommand::Optimize, &mut config),
      format!(
        "input `{}` of setting `tree` is outside the directories the server reads inputs from",
        escape.display()
      )
    );
  }

  #[test]
  fn test_confine_rejects_absolute_path_outside_roots() {
    let Dirs { data, outside } = dirs();
    let policy = PathPolicy::new(&[data.path().to_path_buf()]).unwrap();
    let secret = outside.path().join("secret.nwk");
    let mut config = json!({ "metadata": secret, "tree": data.path().join("sub/tree.nwk") });
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
    let link = data.path().join("link.nwk");
    symlink(outside.path().join("secret.nwk"), &link).unwrap();
    let policy = PathPolicy::new(&[data.path().to_path_buf()]).unwrap();
    let mut config = json!({ "tree": link });
    assert_error!(
      policy.confine(AppCommand::Prune, &mut config),
      format!(
        "input `{}` of setting `tree` is outside the directories the server reads inputs from",
        link.display()
      )
    );
  }

  #[test]
  fn test_confine_rejects_symlinked_directory_that_resolves_outside() {
    let Dirs { data, outside } = dirs();
    symlink(outside.path(), data.path().join("linked")).unwrap();
    let policy = PathPolicy::new(&[data.path().to_path_buf()]).unwrap();
    let secret = data.path().join("linked/secret.nwk");
    let mut config = json!({ "alignment": [data.path().join("sub/aln.fasta"), secret] });
    assert_error!(
      policy.confine(AppCommand::Ancestral, &mut config),
      format!(
        "input `{}` of setting `alignment` is outside the directories the server reads inputs from",
        secret.display()
      )
    );
  }

  #[test]
  fn test_confine_accepts_symlink_that_stays_inside() {
    let Dirs { data, .. } = dirs();
    let alias = data.path().join("alias.nwk");
    symlink(data.path().join("sub/tree.nwk"), &alias).unwrap();
    let policy = PathPolicy::new(&[data.path().to_path_buf()]).unwrap();
    let mut config = json!({ "tree": alias });
    policy.confine(AppCommand::Prune, &mut config).unwrap();
    assert_eq!(
      json!(data.path().canonicalize().unwrap().join("sub/tree.nwk")),
      config["tree"]
    );
  }

  #[test]
  fn test_confine_accepts_inputs_of_any_root() {
    let Dirs { data, outside } = dirs();
    let policy = PathPolicy::new(&[data.path().to_path_buf(), outside.path().to_path_buf()]).unwrap();
    let secret = outside.path().join("secret.nwk");
    let mut config = json!({ "tree": secret });
    policy.confine(AppCommand::Prune, &mut config).unwrap();
    assert_eq!(json!(secret.canonicalize().unwrap()), config["tree"]);
  }

  #[test]
  fn test_confine_missing_input_is_an_error() {
    let Dirs { data, .. } = dirs();
    let missing = data.path().join("missing.nwk");
    let policy = PathPolicy::new(&[data.path().to_path_buf()]).unwrap();
    let mut config = json!({ "tree": missing });
    assert_error!(
      policy.confine(AppCommand::Prune, &mut config),
      format!(
        "input `{}` of setting `tree` cannot be read: No such file or directory (os error 2)",
        missing.display()
      )
    );
  }

  #[test]
  fn test_confine_translation_template_expands_listed_cdses() {
    let Dirs { data, outside } = dirs();
    let root = data.path().canonicalize().unwrap();
    fs::write(root.join("sub/translation_E.fasta"), ">a\nM\n").unwrap();
    let policy = PathPolicy::new(slice::from_ref(&root)).unwrap();
    let mut config = json!({
      "tree": root.join("sub/tree.nwk"),
      "translations": root.join("sub/translation_{cds}.fasta"),
      "cdses": ["E"],
    });
    policy.confine(AppCommand::Ancestral, &mut config).unwrap();

    let escape = format!("../../{}/secret", outside.path().file_name().unwrap().to_string_lossy());
    let mut config = json!({
      "tree": root.join("sub/tree.nwk"),
      "translations": root.join("sub/{cds}.nwk"),
      "cdses": [escape],
    });
    assert_error!(
      policy.confine(AppCommand::Ancestral, &mut config),
      format!(
        "input `{}/sub/{escape}.nwk` of setting `translations` is outside the directories the server reads inputs from",
        root.display()
      )
    );
  }

  #[test]
  fn test_confine_relative_translation_template_is_rewritten_to_the_working_directory() {
    let Dirs { data, .. } = dirs();
    let policy = PathPolicy::new(&[data.path().to_path_buf(), Path::new("src").to_path_buf()]).unwrap();
    let mut config = json!({
      "tree": data.path().join("sub/tree.nwk"),
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
    let policy = PathPolicy::new(&[data.path().to_path_buf()]).unwrap();
    let mut config = json!({
      "tree": data.path().join("sub/tree.nwk"),
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
  fn test_confine_rejects_missing_input_root() {
    let Dirs { data, .. } = dirs();
    let missing = data.path().join("missing");
    assert_error!(
      PathPolicy::new(slice::from_ref(&missing)),
      format!(
        "When resolving an input directory: When resolving directory '{}': No such file or directory (os error 2)",
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
