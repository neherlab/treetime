#[cfg(test)]
mod tests {
  use crate::commands::shared::output::{CommandKind, NwkStyleArg, OutputCoreArgs, OutputSelection};
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use tempfile::TempDir;
  use treetime_io::graph::TreeWriteKind;
  use treetime_io::nwk::NwkStyle;
  use treetime_utils::assert_error;

  fn nwk(style: NwkStyle) -> TreeWriteKind {
    TreeWriteKind::nwk(style)
  }

  fn nexus(style: NwkStyle) -> TreeWriteKind {
    TreeWriteKind::nexus(style)
  }

  // --- Tier resolution behavior ---

  #[test]
  fn test_resolve_output_all_default_selection() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Ancestral, &[], &[]).unwrap();

    assert!(resolved.tree_outputs.contains_key(&nwk(NwkStyle::Plain)));
    assert!(resolved.tree_outputs.contains_key(&nexus(NwkStyle::Plain)));
    assert!(
      !resolved.tree_outputs.contains_key(&TreeWriteKind::GraphJson),
      "GraphJson is not a default output"
    );
    assert!(
      !resolved.tree_outputs.contains_key(&TreeWriteKind::Dot),
      "Dot is not a default output"
    );
  }

  #[test]
  fn test_resolve_output_all_with_selection_filter() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = args
      .resolve(CommandKind::Ancestral, &[OutputSelection::Nwk], &[])
      .unwrap();

    assert_eq!(resolved.tree_outputs.len(), 1);
    assert!(resolved.tree_outputs.contains_key(&nwk(NwkStyle::Plain)));
  }

  #[test]
  fn test_resolve_per_file_only_no_output_all() {
    let args = OutputCoreArgs {
      output_tree_nwk: Some(PathBuf::from("/tmp/test.nwk")),
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Ancestral, &[], &[]).unwrap();

    assert_eq!(
      resolved.tree_outputs,
      btreemap! { nwk(NwkStyle::Plain) => PathBuf::from("/tmp/test.nwk") }
    );
  }

  #[test]
  fn test_resolve_no_outputs_errors() {
    let args = OutputCoreArgs::default();
    let result = args.resolve(CommandKind::Ancestral, &[], &[]);
    assert_error!(
      result,
      "No output flags provided. At least one is required: --output-all or one of the --output-tree-* / --output-* flags"
    );
  }

  #[test]
  fn test_resolve_rejects_duplicate_output_destinations() {
    let path = PathBuf::from("same-output");
    let args = OutputCoreArgs {
      output_tree_nwk: Some(path.clone()),
      ..OutputCoreArgs::default()
    };

    let result = args.resolve(
      CommandKind::Ancestral,
      &[],
      &[(OutputSelection::Gtr, Some(path.as_path()))],
    );

    assert_error!(
      result,
      "Output destination 'same-output' is selected more than once (Nwk(NwkWriteSpec { style: Plain }) and --output-gtr)"
    );
  }

  #[test]
  fn test_resolve_non_tree_explicit_path() {
    let dir = TempDir::new().unwrap();
    let gtr_path = dir.path().join("my.gtr.json");
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = args
      .resolve(
        CommandKind::Ancestral,
        &[],
        &[(OutputSelection::Gtr, Some(gtr_path.as_path()))],
      )
      .unwrap();

    assert_eq!(&resolved.non_tree_outputs[&OutputSelection::Gtr], &gtr_path);
  }

  #[test]
  fn test_resolve_non_tree_default_path_uses_command_stem() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = args
      .resolve(CommandKind::Ancestral, &[], &[(OutputSelection::Gtr, None)])
      .unwrap();

    assert_eq!(
      &resolved.non_tree_outputs[&OutputSelection::Gtr],
      &dir.path().join("ancestral.gtr.json")
    );
  }

  #[test]
  fn test_resolve_tree_default_paths_use_command_stem() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Clock, &[OutputSelection::Nwk], &[]).unwrap();
    assert_eq!(
      &resolved.tree_outputs[&nwk(NwkStyle::Plain)],
      &dir.path().join("clock.nwk")
    );
  }

  #[test]
  fn test_resolve_timetree_defaults_include_auspice() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Timetree, &[], &[]).unwrap();
    assert!(resolved.tree_outputs.contains_key(&TreeWriteKind::Auspice));
  }

  // --- All-selection expansion ---

  #[rustfmt::skip]
  #[rstest]
  #[case::ancestral(CommandKind::Ancestral)]
  #[case::timetree (CommandKind::Timetree)]
  #[case::optimize (CommandKind::Optimize)]
  #[case::mugration(CommandKind::Mugration)]
  #[case::clock    (CommandKind::Clock)]
  #[case::prune    (CommandKind::Prune)]
  #[trace]
  fn test_resolve_all_expands_to_complete_tree_set(#[case] command: CommandKind) {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = args.resolve(command, &[OutputSelection::All], &[]).unwrap();
    let actual = resolved.tree_outputs.keys().cloned().collect();
    let expected = btreeset! {
      nwk(NwkStyle::Plain),
      nexus(NwkStyle::Plain),
      TreeWriteKind::Auspice,
      TreeWriteKind::Phyloxml,
      TreeWriteKind::PhyloxmlJson,
      TreeWriteKind::MatPb,
      TreeWriteKind::MatJson,
      TreeWriteKind::GraphJson,
      TreeWriteKind::Dot,
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_resolve_explicit_mat_available_on_clock() {
    let args = OutputCoreArgs {
      output_tree_mat_pb: Some(PathBuf::from("/tmp/out.pb")),
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Clock, &[], &[]).unwrap();
    assert_eq!(
      &PathBuf::from("/tmp/out.pb"),
      &resolved.tree_outputs[&TreeWriteKind::MatPb]
    );
  }

  #[test]
  fn test_resolve_explicit_auspice_available_on_clock() {
    let args = OutputCoreArgs {
      output_tree_auspice: Some(PathBuf::from("/tmp/out.auspice.json")),
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Clock, &[], &[]).unwrap();
    assert_eq!(
      &PathBuf::from("/tmp/out.auspice.json"),
      &resolved.tree_outputs[&TreeWriteKind::Auspice]
    );
  }

  // --- S2 interaction table (acceptance criteria) ---

  #[test]
  fn test_s2_row1_output_all_default_style() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = args
      .resolve(
        CommandKind::Ancestral,
        &[OutputSelection::Nwk, OutputSelection::Nexus],
        &[],
      )
      .unwrap();
    let expected = btreemap! {
      nwk(NwkStyle::Plain) => dir.path().join("ancestral.nwk"),
      nexus(NwkStyle::Plain) => dir.path().join("ancestral.nexus"),
    };
    assert_eq!(resolved.tree_outputs, expected);
  }

  #[test]
  fn test_s2_row2_output_all_two_styles() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      output_nwk_style: vec![NwkStyleArg::Plain, NwkStyleArg::Beast],
      ..Default::default()
    };
    let resolved = args
      .resolve(
        CommandKind::Ancestral,
        &[OutputSelection::Nwk, OutputSelection::Nexus],
        &[],
      )
      .unwrap();
    let expected = btreemap! {
      nwk(NwkStyle::Plain) => dir.path().join("ancestral.nwk"),
      nwk(NwkStyle::Beast) => dir.path().join("ancestral.annotated.nwk"),
      nexus(NwkStyle::Plain) => dir.path().join("ancestral.nexus"),
      nexus(NwkStyle::Beast) => dir.path().join("ancestral.annotated.nexus"),
    };
    assert_eq!(resolved.tree_outputs, expected);
  }

  #[test]
  fn test_s2_row3_per_file_default_style() {
    let args = OutputCoreArgs {
      output_tree_nwk: Some(PathBuf::from("my.nwk")),
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Ancestral, &[], &[]).unwrap();
    assert_eq!(
      resolved.tree_outputs,
      btreemap! { nwk(NwkStyle::Plain) => PathBuf::from("my.nwk") }
    );
  }

  #[test]
  fn test_s2_row4_per_file_single_non_plain_style_used_as_is() {
    let args = OutputCoreArgs {
      output_tree_nwk: Some(PathBuf::from("my.nwk")),
      output_nwk_style: vec![NwkStyleArg::Beast],
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Ancestral, &[], &[]).unwrap();
    assert_eq!(
      resolved.tree_outputs,
      btreemap! { nwk(NwkStyle::Beast) => PathBuf::from("my.nwk") }
    );
  }

  #[test]
  fn test_s2_row5_per_file_two_styles_insert_secondary() {
    let args = OutputCoreArgs {
      output_tree_nwk: Some(PathBuf::from("my.nwk")),
      output_nwk_style: vec![NwkStyleArg::Plain, NwkStyleArg::Beast],
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Ancestral, &[], &[]).unwrap();
    let expected = btreemap! {
      nwk(NwkStyle::Plain) => PathBuf::from("my.nwk"),
      nwk(NwkStyle::Beast) => PathBuf::from("my.annotated.nwk"),
    };
    assert_eq!(resolved.tree_outputs, expected);
  }

  #[test]
  fn test_s2_row6_per_file_overrides_output_all_nwk_only() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      output_tree_nwk: Some(PathBuf::from("my.nwk")),
      output_nwk_style: vec![NwkStyleArg::Plain, NwkStyleArg::Beast],
      ..Default::default()
    };
    let resolved = args
      .resolve(
        CommandKind::Ancestral,
        &[OutputSelection::Nwk, OutputSelection::Nexus],
        &[],
      )
      .unwrap();
    let expected = btreemap! {
      nwk(NwkStyle::Plain) => PathBuf::from("my.nwk"),
      nwk(NwkStyle::Beast) => PathBuf::from("my.annotated.nwk"),
      nexus(NwkStyle::Plain) => dir.path().join("ancestral.nexus"),
      nexus(NwkStyle::Beast) => dir.path().join("ancestral.annotated.nexus"),
    };
    assert_eq!(resolved.tree_outputs, expected);
  }

  #[test]
  fn test_s2_row7_per_file_supplements_selection() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      output_tree_nwk: Some(PathBuf::from("my.nwk")),
      ..Default::default()
    };
    let resolved = args
      .resolve(CommandKind::Ancestral, &[OutputSelection::Nexus], &[])
      .unwrap();
    let expected = btreemap! {
      nwk(NwkStyle::Plain) => PathBuf::from("my.nwk"),
      nexus(NwkStyle::Plain) => dir.path().join("ancestral.nexus"),
    };
    assert_eq!(resolved.tree_outputs, expected);
  }

  #[test]
  fn test_s2_row8_two_per_file_two_non_plain_styles() {
    let args = OutputCoreArgs {
      output_tree_nwk: Some(PathBuf::from("my.nwk")),
      output_tree_nexus: Some(PathBuf::from("my.nexus")),
      output_nwk_style: vec![NwkStyleArg::Beast, NwkStyleArg::Nhx],
      ..Default::default()
    };
    let resolved = args.resolve(CommandKind::Ancestral, &[], &[]).unwrap();
    let expected: BTreeMap<TreeWriteKind, PathBuf> = btreemap! {
      nwk(NwkStyle::Beast) => PathBuf::from("my.annotated.nwk"),
      nwk(NwkStyle::Nhx) => PathBuf::from("my.nhx.nwk"),
      nexus(NwkStyle::Beast) => PathBuf::from("my.annotated.nexus"),
      nexus(NwkStyle::Nhx) => PathBuf::from("my.nhx.nexus"),
    };
    assert_eq!(resolved.tree_outputs, expected);
  }
}
