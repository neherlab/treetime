#[cfg(test)]
mod tests {
  use crate::commands::shared::output::{CommandKind, OutputSelection};
  use pretty_assertions::assert_eq;
  use strum::IntoEnumIterator;

  #[test]
  fn test_output_selection_is_tree_covers_all_tree_variants() {
    let tree_variants: Vec<OutputSelection> = OutputSelection::iter().filter(|s| s.is_tree()).collect();
    assert_eq!(
      tree_variants,
      vec![
        OutputSelection::Nwk,
        OutputSelection::NwkAnnotated,
        OutputSelection::NwkNhx,
        OutputSelection::Nexus,
        OutputSelection::NexusAnnotated,
        OutputSelection::NexusNhx,
        OutputSelection::Auspice,
        OutputSelection::Phyloxml,
        OutputSelection::PhyloxmlJson,
        OutputSelection::MatPb,
        OutputSelection::MatJson,
        OutputSelection::GraphJson,
        OutputSelection::Dot,
      ]
    );
  }

  #[test]
  fn test_output_selection_non_tree_variants() {
    let non_tree: Vec<OutputSelection> = OutputSelection::iter()
      .filter(|s| !s.is_tree() && !s.is_meta())
      .collect();
    assert_eq!(
      non_tree,
      vec![
        OutputSelection::AugurNodeData,
        OutputSelection::Gtr,
        OutputSelection::ClockModel,
        OutputSelection::Confidence,
        OutputSelection::ReconstructedNucFasta,
        OutputSelection::ReconstructedAaFasta,
        OutputSelection::TraitsCsv,
        OutputSelection::ClockCsv,
      ]
    );
  }

  #[test]
  fn test_output_selection_all_is_meta() {
    assert!(OutputSelection::All.is_meta());
    assert!(!OutputSelection::All.is_tree());
  }

  #[test]
  fn test_output_selection_every_tree_variant_has_tree_write_kind() {
    for variant in OutputSelection::iter() {
      if variant.is_tree() {
        assert!(
          variant.to_tree_write_kind().is_some(),
          "{variant:?} is_tree() but to_tree_write_kind() is None"
        );
      }
    }
  }

  #[test]
  fn test_output_selection_non_tree_has_no_tree_write_kind() {
    for variant in OutputSelection::iter() {
      if !variant.is_tree() {
        assert!(
          variant.to_tree_write_kind().is_none(),
          "{variant:?} is not a tree but to_tree_write_kind() returned Some"
        );
      }
    }
  }

  #[test]
  fn test_output_selection_every_variant_has_extension() {
    for variant in OutputSelection::iter() {
      if !variant.is_meta() {
        assert!(!variant.extension().is_empty(), "{variant:?} has empty extension");
      }
    }
  }

  #[test]
  fn test_output_selection_every_variant_has_flag_name() {
    for variant in OutputSelection::iter() {
      let flag = variant.flag_name();
      assert!(!flag.is_empty(), "{variant:?} has empty flag_name");
      assert!(flag.starts_with("--"), "{variant:?} flag_name doesn't start with --");
    }
  }

  #[test]
  fn test_command_kind_default_outputs_subset_of_available() {
    for command in [
      CommandKind::Ancestral,
      CommandKind::Timetree,
      CommandKind::Optimize,
      CommandKind::Mugration,
      CommandKind::Clock,
      CommandKind::Prune,
    ] {
      let available = command.available_outputs();
      let defaults = command.default_outputs();
      assert!(
        defaults.is_subset(&available),
        "{command:?}: default outputs {:?} not subset of available {:?}",
        defaults.difference(&available).collect::<Vec<_>>(),
        available
      );
    }
  }

  #[test]
  fn test_command_kind_default_tree_outputs_always_include_nwk_nexus() {
    for command in [
      CommandKind::Ancestral,
      CommandKind::Timetree,
      CommandKind::Optimize,
      CommandKind::Mugration,
      CommandKind::Clock,
      CommandKind::Prune,
    ] {
      let defaults = command.default_outputs();
      assert!(
        defaults.contains(&OutputSelection::Nwk),
        "{command:?} defaults missing Nwk"
      );
      assert!(
        defaults.contains(&OutputSelection::Nexus),
        "{command:?} defaults missing Nexus"
      );
    }
  }

  #[test]
  fn test_command_kind_stems_non_empty() {
    for command in [
      CommandKind::Ancestral,
      CommandKind::Timetree,
      CommandKind::Optimize,
      CommandKind::Mugration,
      CommandKind::Clock,
      CommandKind::Prune,
    ] {
      assert!(!command.stem().is_empty(), "{command:?} has empty stem");
    }
  }

  #[test]
  fn test_command_kind_timetree_defaults_include_auspice() {
    let defaults = CommandKind::Timetree.default_outputs();
    assert!(defaults.contains(&OutputSelection::Auspice));
  }

  #[test]
  fn test_command_kind_non_timetree_defaults_exclude_auspice() {
    for command in [
      CommandKind::Ancestral,
      CommandKind::Optimize,
      CommandKind::Mugration,
      CommandKind::Clock,
      CommandKind::Prune,
    ] {
      let defaults = command.default_outputs();
      assert!(
        !defaults.contains(&OutputSelection::Auspice),
        "{command:?} defaults should not include Auspice"
      );
    }
  }

  #[test]
  fn test_command_kind_graph_json_dot_not_in_defaults() {
    for command in [
      CommandKind::Ancestral,
      CommandKind::Timetree,
      CommandKind::Optimize,
      CommandKind::Mugration,
      CommandKind::Clock,
      CommandKind::Prune,
    ] {
      let defaults = command.default_outputs();
      assert!(
        !defaults.contains(&OutputSelection::GraphJson),
        "{command:?} defaults should not include GraphJson"
      );
      assert!(
        !defaults.contains(&OutputSelection::Dot),
        "{command:?} defaults should not include Dot"
      );
    }
  }

  #[test]
  fn test_command_kind_aa_fasta_available_but_not_default() {
    let available = CommandKind::Ancestral.available_outputs();
    let defaults = CommandKind::Ancestral.default_outputs();
    assert!(available.contains(&OutputSelection::ReconstructedAaFasta));
    assert!(!defaults.contains(&OutputSelection::ReconstructedAaFasta));
  }

  #[test]
  fn test_command_kind_auspice_not_available_on_non_timetree() {
    for command in [
      CommandKind::Ancestral,
      CommandKind::Optimize,
      CommandKind::Mugration,
      CommandKind::Clock,
      CommandKind::Prune,
    ] {
      let available = command.available_outputs();
      assert!(
        !available.contains(&OutputSelection::Auspice),
        "{command:?} should not have Auspice in available outputs"
      );
    }
  }

  #[test]
  fn test_command_kind_phyloxml_mat_not_available() {
    for command in [
      CommandKind::Ancestral,
      CommandKind::Timetree,
      CommandKind::Optimize,
      CommandKind::Mugration,
      CommandKind::Clock,
      CommandKind::Prune,
    ] {
      let available = command.available_outputs();
      assert!(!available.contains(&OutputSelection::Phyloxml));
      assert!(!available.contains(&OutputSelection::PhyloxmlJson));
      assert!(!available.contains(&OutputSelection::MatPb));
      assert!(!available.contains(&OutputSelection::MatJson));
    }
  }
}
