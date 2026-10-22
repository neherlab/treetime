#[cfg(test)]
mod tests {
  use crate::{
    DatasetFiles, ExampleFile, discover_datasets, example_configs, group_datasets, parse_example_config,
    schema_directive, text_schema_command,
  };
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::path::Path;
  use treetime_utils::{o, vec_of_owned};

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
  fn test_examples_discovery_names_the_examples_dir_as_given() {
    let catalog = discover_datasets(Path::new("../../data"), COMMANDS).unwrap();
    assert_eq!("../../data", catalog.examples_dir);
  }

  #[test]
  fn test_examples_grouping_of_a_file_list_keeps_folders_with_a_tree() {
    let paths = [
      "README.md",
      "zika/country_weights.csv",
      "zika/20/tree.nwk",
      "zika/20/metadata.tsv",
      "zika/20/.hidden.nwk",
      "usher/.gitignore",
      "usher/global.pb",
      "tree.nwk",
      "sc2/4500/tree.nwk",
      "sc2/4500/metadata.tsv.xz",
      "sc2/4500/translations/F.fasta.xz",
    ];
    assert_eq!(
      vec![
        DatasetFiles {
          name: o!(""),
          files: vec_of_owned!["README.md", "tree.nwk"],
        },
        DatasetFiles {
          name: o!("sc2/4500"),
          files: vec_of_owned!["metadata.tsv.xz", "tree.nwk"],
        },
        DatasetFiles {
          name: o!("zika/20"),
          files: vec_of_owned!["metadata.tsv", "tree.nwk"],
        },
      ],
      group_datasets(paths)
    );
  }

  #[test]
  fn test_examples_configs_of_a_file_list_read_only_yaml_files() {
    let directive = schema_directive("clock");
    let paths = [
      "b/clock.yaml",
      "a/notes.txt",
      "a/.hidden.yaml",
      "a/clock.yml",
      "a/plain.yaml",
    ];
    let mut read = vec![];
    let examples = example_configs(
      paths,
      |path| {
        read.push(path.to_owned());
        Ok::<_, Report>(match path {
          "a/plain.yaml" => o!("tree: tree.nwk\n"),
          _ => format!("{directive}\n# Clock of {path}.\ntree: tree.nwk\n"),
        })
      },
      COMMANDS,
    )
    .unwrap();
    assert_eq!(
      (
        vec_of_owned!["a/clock.yml", "b/clock.yaml"],
        vec_of_owned!["a/clock.yml", "a/plain.yaml", "b/clock.yaml"]
      ),
      (
        examples.into_iter().map(|example| example.path).collect::<Vec<_>>(),
        read
      )
    );
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
      Some(ExampleFile {
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

  #[test]
  fn test_examples_schema_directive_names_the_published_schema_of_the_command() {
    assert_eq!(
      "# yaml-language-server: $schema=https://raw.githubusercontent.com/neherlab/treetime/rust/packages/schemas/input-config-clock.schema.json",
      schema_directive("clock")
    );
  }

  #[test]
  fn test_examples_schema_directive_is_read_back_as_its_command() {
    let content = format!("\n{}\ntree: tree.nwk\n", schema_directive("prune"));
    assert_eq!(Some("prune"), text_schema_command(&content));
  }

  #[test]
  fn test_examples_text_without_a_leading_directive_names_no_command() {
    let content = format!("tree: tree.nwk\n{}\n", schema_directive("prune"));
    assert_eq!(None, text_schema_command(&content));
  }
}
