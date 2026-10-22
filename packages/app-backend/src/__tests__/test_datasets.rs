#[cfg(test)]
mod tests {
  use crate::check_inputs::InputKind;
  use crate::datasets::{
    DataIndex, DatasetInput, IndexedConfig, IndexedDataset, IndexedFile, IndexedInput, IndexedPipeline, IndexedStep,
    dataset_catalog,
  };
  use app_commands::command::AppCommand;
  use helpers::{config_text, index_of};
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::path::Path;
  use treetime_utils::{assert_error, o, vec_of_owned};

  #[test]
  fn test_datasets_zika_86_inputs_by_kind() {
    let catalog = dataset_catalog(Path::new("../../data")).unwrap();
    let zika = catalog
      .datasets
      .iter()
      .find(|dataset| dataset.name == "zika/86")
      .unwrap();
    assert_eq!(
      vec![
        DatasetInput {
          kind: InputKind::Tree,
          file: o!("zika/86/tree.nwk"),
          path: o!("../../data/zika/86/tree.nwk"),
        },
        DatasetInput {
          kind: InputKind::Alignment,
          file: o!("zika/86/aln.fasta.xz"),
          path: o!("../../data/zika/86/aln.fasta.xz"),
        },
        DatasetInput {
          kind: InputKind::Metadata,
          file: o!("zika/86/metadata.tsv"),
          path: o!("../../data/zika/86/metadata.tsv"),
        },
      ],
      zika.inputs
    );
  }

  #[test]
  fn test_datasets_examples_name_their_command() {
    let catalog = dataset_catalog(Path::new("../../data")).unwrap();
    let example = catalog
      .examples
      .iter()
      .find(|example| example.path == "ebola/20/ancestral-parsimony.yaml")
      .unwrap();
    assert_eq!(AppCommand::Ancestral, example.command);
  }

  #[test]
  fn test_datasets_index_maps_compressed_inputs_and_lists_every_file_but_itself() {
    let index = index_of(
      &[
        "index.json",
        "sc2/4500/tree.nwk",
        "sc2/4500/metadata.tsv.xz",
        "sc2/4500/aln.fasta.xz",
        "sc2/4500/annotation.gff3",
        "zika/86/tree.nwk",
        "zika/86/zika.phylip.xz",
      ],
      &[],
    )
    .unwrap();
    assert_eq!(
      DataIndex {
        files: vec![
          IndexedFile {
            path: o!("sc2/4500/aln.fasta.xz"),
            size: 21,
          },
          IndexedFile {
            path: o!("sc2/4500/annotation.gff3"),
            size: 24,
          },
          IndexedFile {
            path: o!("sc2/4500/metadata.tsv.xz"),
            size: 24,
          },
          IndexedFile {
            path: o!("sc2/4500/tree.nwk"),
            size: 17,
          },
          IndexedFile {
            path: o!("zika/86/tree.nwk"),
            size: 16,
          },
          IndexedFile {
            path: o!("zika/86/zika.phylip.xz"),
            size: 22,
          },
        ],
        datasets: vec![
          IndexedDataset {
            name: o!("sc2/4500"),
            inputs: vec![
              IndexedInput {
                kind: InputKind::Tree,
                path: o!("sc2/4500/tree.nwk"),
              },
              IndexedInput {
                kind: InputKind::Alignment,
                path: o!("sc2/4500/aln.fasta.xz"),
              },
              IndexedInput {
                kind: InputKind::Metadata,
                path: o!("sc2/4500/metadata.tsv.xz"),
              },
            ],
          },
          IndexedDataset {
            name: o!("zika/86"),
            inputs: vec![IndexedInput {
              kind: InputKind::Tree,
              path: o!("zika/86/tree.nwk"),
            }],
          },
        ],
        configs: vec![],
        pipelines: vec![],
      },
      index
    );
  }

  #[test]
  fn test_datasets_index_leaves_out_hidden_files_and_hidden_folders() {
    let index = index_of(&["usher/.gitignore", "usher/.cache/tree.nwk", "usher/README.md"], &[]).unwrap();
    assert_eq!(
      DataIndex {
        files: vec![IndexedFile {
          path: o!("usher/README.md"),
          size: 15,
        }],
        datasets: vec![],
        configs: vec![],
        pipelines: vec![],
      },
      index
    );
  }

  #[test]
  fn test_datasets_index_resolves_config_inputs_from_the_config_folder() {
    let content = config_text(
      "mugration",
      indoc! {r#"
        tree: "tree.nwk"
        metadata: "./metadata.tsv"
        attribute: "country"
        weights: "../country_weights.csv"
      "#},
    );
    let index = index_of(
      &[
        "zika/country_weights.csv",
        "zika/20/tree.nwk",
        "zika/20/metadata.tsv",
        "zika/20/mugration.yaml",
      ],
      &[("zika/20/mugration.yaml", &content)],
    )
    .unwrap();
    assert_eq!(
      vec![IndexedConfig {
        path: o!("zika/20/mugration.yaml"),
        command: AppCommand::Mugration,
        title: o!("Example of mugration."),
        folder: o!("zika/20"),
        content,
        inputs: btreemap! {
          o!("metadata") => vec_of_owned!["zika/20/metadata.tsv"],
          o!("tree") => vec_of_owned!["zika/20/tree.nwk"],
          o!("weights") => vec_of_owned!["zika/country_weights.csv"],
        },
      }],
      index.configs
    );
  }

  #[test]
  fn test_datasets_index_lists_every_file_of_a_list_setting() {
    let content = config_text(
      "ancestral",
      indoc! {r#"
        tree: "tree.nwk"
        alignment:
          - "aln.fasta.xz"
          - "../shared/extra.fasta"
      "#},
    );
    let index = index_of(
      &["a/tree.nwk", "a/aln.fasta.xz", "a/ancestral.yaml", "shared/extra.fasta"],
      &[("a/ancestral.yaml", &content)],
    )
    .unwrap();
    assert_eq!(
      vec![btreemap! {
        o!("alignment") => vec_of_owned!["a/aln.fasta.xz", "shared/extra.fasta"],
        o!("tree") => vec_of_owned!["a/tree.nwk"],
      }],
      index
        .configs
        .into_iter()
        .map(|config| config.inputs)
        .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_datasets_index_rejects_a_config_input_without_a_tracked_file() {
    let content = config_text("clock", "tree: \"tree.nwk\"\nmetadata: \"metadata.tsv\"\n");
    assert_error!(
      index_of(&["a/tree.nwk", "a/clock.yaml"], &[("a/clock.yaml", &content)]),
      "When indexing the input files of example configuration 'a/clock.yaml': When reading the input setting `metadata`: the path 'a/metadata.tsv' names no tracked file of the data folder"
    );
  }

  #[test]
  fn test_datasets_index_rejects_a_config_input_outside_the_data_folder() {
    let content = config_text("clock", "tree: \"../../tree.nwk\"\n");
    assert_error!(
      index_of(&["a/tree.nwk", "a/clock.yaml"], &[("a/clock.yaml", &content)]),
      "When indexing the input files of example configuration 'a/clock.yaml': When reading the input setting `tree`: the path 'a/../../tree.nwk' leaves the data folder"
    );
  }

  #[test]
  fn test_datasets_index_rejects_an_absolute_config_input() {
    let content = config_text("clock", "tree: \"/data/a/tree.nwk\"\n");
    assert_error!(
      index_of(&["a/tree.nwk", "a/clock.yaml"], &[("a/clock.yaml", &content)]),
      "When indexing the input files of example configuration 'a/clock.yaml': When reading the input setting `tree`: the path '/data/a/tree.nwk' is outside the data folder"
    );
  }

  #[test]
  fn test_datasets_index_rejects_a_config_with_an_input_template() {
    let content = config_text(
      "ancestral",
      "tree: \"tree.nwk\"\ntranslations: \"translations/{cds}.fasta\"\n",
    );
    assert_error!(
      index_of(&["a/tree.nwk", "a/ancestral.yaml"], &[("a/ancestral.yaml", &content)]),
      "When indexing the input files of example configuration 'a/ancestral.yaml': the input template setting `translations` names files by a placeholder, which the index cannot list"
    );
  }

  #[test]
  fn test_datasets_index_reads_only_the_listed_files() {
    let config = config_text("clock", "tree: \"tree.nwk\"\n");
    let index = index_of(
      &["a/tree.nwk", "a/clock.yaml"],
      &[("a/clock.yaml", &config), ("usher/clock.yaml", &config)],
    )
    .unwrap();
    assert_eq!(
      (
        vec_of_owned!["a/clock.yaml", "a/tree.nwk"],
        vec_of_owned!["a"],
        vec_of_owned!["a/clock.yaml"]
      ),
      (
        index.files.into_iter().map(|file| file.path).collect::<Vec<_>>(),
        index
          .datasets
          .into_iter()
          .map(|dataset| dataset.name)
          .collect::<Vec<_>>(),
        index.configs.into_iter().map(|config| config.path).collect::<Vec<_>>()
      )
    );
  }

  #[rstest]
  #[case::parent("../secret.txt")]
  #[case::absolute("/etc/passwd")]
  #[case::not_normalized("a/./tree.nwk")]
  #[case::empty("")]
  #[trace]
  fn test_datasets_index_rejects_a_file_list_path_outside_the_data_folder(#[case] path: &str) {
    assert_error!(
      index_of(&["a/tree.nwk", path], &[]),
      format!("the file list names '{path}', which is not a normalized path relative to the data folder")
    );
  }

  #[test]
  fn test_datasets_index_lists_pipelines_with_their_steps() {
    let content = config_text("pipeline", "steps: []\n");
    let index = DataIndex::new(
      vec_of_owned!["a/pipeline.yaml"],
      |_| Ok(1),
      |_| Ok(content.clone()),
      |text| {
        assert_eq!(content, text);
        Ok(vec![IndexedStep {
          name: o!("first"),
          command: AppCommand::Timetree,
        }])
      },
    )
    .unwrap();
    assert_eq!(
      (
        vec![IndexedPipeline {
          path: o!("a/pipeline.yaml"),
          title: o!("Example of pipeline."),
          steps: vec![IndexedStep {
            name: o!("first"),
            command: AppCommand::Timetree,
          }],
        }],
        vec![]
      ),
      (index.pipelines, index.configs)
    );
  }

  mod helpers {
    use crate::datasets::DataIndex;
    use app_datasets::schema_directive;
    use eyre::Report;
    use treetime_utils::{make_error, make_internal_report};

    pub(super) fn index_of(paths: &[&str], contents: &[(&str, &str)]) -> Result<DataIndex, Report> {
      DataIndex::new(
        paths.iter().map(|path| (*path).to_owned()).collect(),
        |path| Ok(u64::try_from(path.len()).unwrap()),
        |path| {
          contents
            .iter()
            .find(|(name, _)| *name == path)
            .map(|(_, content)| (*content).to_owned())
            .ok_or_else(|| make_internal_report!("the test has no content for '{path}'"))
        },
        |_| make_error!("the test reads no pipeline"),
      )
    }

    pub(super) fn config_text(command: &str, settings: &str) -> String {
      let directive = schema_directive(command);
      format!("{directive}\n# Example of {command}.\n{settings}")
    }
  }
}
