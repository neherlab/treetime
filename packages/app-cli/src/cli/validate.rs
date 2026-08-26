use eyre::Report;
use treetime::commands::ancestral::args::TreetimeAncestralArgs;
use treetime::commands::clock::args::TreetimeClockArgs;
use treetime::commands::homoplasy::args::TreetimeHomoplasyArgs;
use treetime::commands::mugration::args::TreetimeMugrationArgs;
use treetime::commands::optimize::args::TreetimeOptimizeArgs;
use treetime::commands::prune::args::TreetimePruneArgs;
use treetime::commands::timetree::args::TreetimeTimetreeArgs;
use treetime_utils::make_error;

/// Validate a command's arguments after CLI parsing and the `--config` merge.
///
/// This is the single place where required-argument presence is checked. The check lives here rather
/// than in clap because `--config` can legitimately supply a value that clap parsing never saw:
/// required fields are `Option<T>` (so clap does not reject them at parse time) and are enforced here,
/// once the merged precedence (explicit CLI > config file > default) has produced the final value.
///
/// Cross-argument interaction rules belong here too, so validation does not scatter across command
/// run code.
pub trait Validate {
  fn validate(&self) -> Result<(), Report>;
}

/// Emit a clap-style "required arguments were not provided" error listing every missing flag.
fn require(missing: &[&str]) -> Result<(), Report> {
  if missing.is_empty() {
    Ok(())
  } else {
    let list = missing.join("\n  ");
    make_error!("the following required arguments were not provided:\n  {list}")
  }
}

impl Validate for TreetimeAncestralArgs {
  fn validate(&self) -> Result<(), Report> {
    let mut missing = Vec::new();
    if self.tree.is_none() {
      missing.push("--tree <TREE>");
    }
    require(&missing)
  }
}

impl Validate for TreetimePruneArgs {
  fn validate(&self) -> Result<(), Report> {
    let mut missing = Vec::new();
    if self.tree.is_none() {
      missing.push("--tree <TREE>");
    }
    require(&missing)
  }
}

impl Validate for TreetimeOptimizeArgs {
  fn validate(&self) -> Result<(), Report> {
    let mut missing = Vec::new();
    if self.tree.is_none() {
      missing.push("--tree <TREE>");
    }
    require(&missing)
  }
}

impl Validate for TreetimeClockArgs {
  fn validate(&self) -> Result<(), Report> {
    let mut missing = Vec::new();
    if self.metadata.is_none() {
      missing.push("--metadata <METADATA>");
    }
    require(&missing)
  }
}

impl Validate for TreetimeMugrationArgs {
  fn validate(&self) -> Result<(), Report> {
    let mut missing = Vec::new();
    if self.metadata.is_none() {
      missing.push("--metadata <METADATA>");
    }
    if self.attribute.is_none() {
      missing.push("--attribute <ATTRIBUTE>");
    }
    require(&missing)
  }
}

impl Validate for TreetimeHomoplasyArgs {
  fn validate(&self) -> Result<(), Report> {
    // Homoplasy flattens the ancestral args, so its required `--tree` lives there.
    self.ancestral_args.validate()
  }
}

impl Validate for TreetimeTimetreeArgs {
  fn validate(&self) -> Result<(), Report> {
    // Timetree has no required-only argument: tree, dates, and alignment are each optional at parse
    // time and their combinations are resolved downstream.
    Ok(())
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use treetime_utils::assert_error;

  #[test]
  fn test_validate_ancestral_missing_tree_errors() {
    let args = TreetimeAncestralArgs::default();
    assert_error!(
      args.validate(),
      "the following required arguments were not provided:\n  --tree <TREE>"
    );
  }

  #[test]
  fn test_validate_ancestral_with_tree_ok() {
    let args = TreetimeAncestralArgs {
      tree: Some("tree.nwk".into()),
      ..Default::default()
    };
    args.validate().unwrap();
  }

  #[test]
  fn test_validate_mugration_missing_both_lists_both() {
    let args = TreetimeMugrationArgs::default();
    assert_error!(
      args.validate(),
      "the following required arguments were not provided:\n  --metadata <METADATA>\n  --attribute <ATTRIBUTE>"
    );
  }

  #[test]
  fn test_validate_mugration_only_attribute_missing_lists_attribute() {
    let args = TreetimeMugrationArgs {
      metadata: Some("metadata.tsv".into()),
      ..Default::default()
    };
    assert_error!(
      args.validate(),
      "the following required arguments were not provided:\n  --attribute <ATTRIBUTE>"
    );
  }

  #[test]
  fn test_validate_mugration_only_metadata_missing_lists_metadata() {
    let args = TreetimeMugrationArgs {
      attribute: Some("country".to_owned()),
      ..Default::default()
    };
    assert_error!(
      args.validate(),
      "the following required arguments were not provided:\n  --metadata <METADATA>"
    );
  }

  #[test]
  fn test_validate_clock_missing_metadata_errors() {
    let args = TreetimeClockArgs::default();
    assert_error!(
      args.validate(),
      "the following required arguments were not provided:\n  --metadata <METADATA>"
    );
  }

  #[test]
  fn test_validate_clock_with_metadata_ok() {
    let args = TreetimeClockArgs {
      metadata: Some("metadata.tsv".into()),
      ..Default::default()
    };
    args.validate().unwrap();
  }

  #[test]
  fn test_validate_optimize_missing_tree_errors() {
    let args = TreetimeOptimizeArgs::default();
    assert_error!(
      args.validate(),
      "the following required arguments were not provided:\n  --tree <TREE>"
    );
  }

  #[test]
  fn test_validate_optimize_with_tree_ok() {
    let args = TreetimeOptimizeArgs {
      tree: Some("tree.nwk".into()),
      ..Default::default()
    };
    args.validate().unwrap();
  }

  #[test]
  fn test_validate_prune_missing_tree_errors() {
    let args = TreetimePruneArgs::default();
    assert_error!(
      args.validate(),
      "the following required arguments were not provided:\n  --tree <TREE>"
    );
  }

  #[test]
  fn test_validate_prune_with_tree_ok() {
    let args = TreetimePruneArgs {
      tree: Some("tree.nwk".into()),
      ..Default::default()
    };
    args.validate().unwrap();
  }

  #[test]
  fn test_validate_homoplasy_missing_tree_errors() {
    let args = TreetimeHomoplasyArgs::default();
    assert_error!(
      args.validate(),
      "the following required arguments were not provided:\n  --tree <TREE>"
    );
  }

  #[test]
  fn test_validate_homoplasy_with_tree_ok() {
    let args = TreetimeHomoplasyArgs {
      ancestral_args: TreetimeAncestralArgs {
        tree: Some("tree.nwk".into()),
        ..Default::default()
      },
      ..Default::default()
    };
    args.validate().unwrap();
  }

  #[test]
  fn test_validate_timetree_defaults_ok() {
    TreetimeTimetreeArgs::default().validate().unwrap();
  }
}
