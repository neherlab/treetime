#[cfg(test)]
mod tests {
  use crate::{ExampleConfig, discover_datasets, parse_example_config};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::path::Path;
  use treetime_utils::o;

  const COMMANDS: &[&str] = &["timetree", "optimize", "prune", "ancestral", "clock", "mugration"];

  #[test]
  fn test_examples_discovery_reads_command_and_title_of_ebola_parsimony() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data");
    let catalog = discover_datasets(&data, COMMANDS).unwrap();
    let example = catalog
      .examples
      .iter()
      .find(|example| example.path == "ebola/20/ancestral-parsimony.yaml")
      .unwrap();
    assert_eq!(
      ("ancestral", "Ancestral reconstruction by Fitch parsimony (ebola/20)."),
      (example.command.as_str(), example.title.as_str())
    );
  }

  #[test]
  fn test_examples_discovery_skips_pipeline_configs_and_covers_every_app_command() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data");
    let catalog = discover_datasets(&data, COMMANDS).unwrap();
    assert!(
      catalog
        .examples
        .iter()
        .all(|example| !example.path.ends_with("pipeline.yaml")),
      "pipeline configs are not offered: {:?}",
      catalog.examples
    );
    let mut commands: Vec<&str> = catalog
      .examples
      .iter()
      .map(|example| example.command.as_str())
      .collect();
    commands.sort_unstable();
    commands.dedup();
    assert_eq!(
      vec!["ancestral", "clock", "mugration", "optimize", "prune", "timetree"],
      commands
    );
  }

  #[test]
  fn test_examples_discovery_lists_datasets_with_a_tree() {
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data");
    let catalog = discover_datasets(&data, COMMANDS).unwrap();
    let zika = catalog
      .datasets
      .iter()
      .find(|dataset| dataset.name == "zika/86")
      .unwrap();
    assert_eq!(
      vec!["aln.fasta.xz", "metadata.tsv", "tree.nwk", "zika.phylip.xz"],
      zika.files
    );
  }

  #[test]
  fn test_examples_discovery_names_the_data_dir_as_given() {
    let catalog = discover_datasets(Path::new("../../data"), COMMANDS).unwrap();
    assert_eq!("../../data", catalog.data_dir);
  }

  #[test]
  fn test_examples_parse_takes_the_first_comment_after_the_directive() {
    let content = indoc! {r#"

      # yaml-language-server: $schema=https://example.org/schemas/input-config-clock.schema.json

      #
      # Root-to-tip regression of dengue/100.
      # More text.
      tree: "data/dengue/100/tree.nwk"
    "#};
    assert_eq!(
      Some(ExampleConfig {
        path: o!("dengue/100/clock.yaml"),
        command: o!("clock"),
        title: o!("Root-to-tip regression of dengue/100."),
        content: content.to_owned(),
      }),
      parse_example_config("dengue/100/clock.yaml", content, COMMANDS)
    );
  }

  #[test]
  fn test_examples_parse_rejects_files_without_the_directive() {
    let content = indoc! {r#"
      # Root-to-tip regression.
      tree: "tree.nwk"
    "#};
    assert_eq!(None, parse_example_config("x.yaml", content, COMMANDS));
  }

  #[test]
  fn test_examples_parse_rejects_commands_the_app_does_not_run() {
    let content = indoc! {r#"
      # yaml-language-server: $schema=https://example.org/input-config-pipeline.schema.json
      # Two steps.
      steps: []
    "#};
    assert_eq!(None, parse_example_config("pipeline.yaml", content, COMMANDS));
  }
}
