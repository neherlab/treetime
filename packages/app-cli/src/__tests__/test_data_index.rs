#[cfg(test)]
mod tests {
  use crate::data_index::{generate_data_index, pipeline_steps};
  use app_backend::check_inputs::InputKind;
  use app_backend::datasets::{DataIndex, IndexedDataset, IndexedFile, IndexedInput, IndexedStep};
  use app_commands::command::AppCommand;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::fs;
  use tempfile::tempdir;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_data_index_regenerates_a_tracked_index_missing_on_disk() {
    let data_dir = tempdir().unwrap();
    fs::create_dir_all(data_dir.path().join("a")).unwrap();
    fs::write(data_dir.path().join("a/tree.nwk"), "(A,B);\n").unwrap();
    let index = generate_data_index(data_dir.path(), "a/tree.nwk\0index.json\0").unwrap();
    assert_eq!(
      DataIndex {
        files: vec![IndexedFile {
          path: o!("a/tree.nwk"),
          size: 7,
        }],
        datasets: vec![IndexedDataset {
          name: o!("a"),
          inputs: vec![IndexedInput {
            kind: InputKind::Tree,
            path: o!("a/tree.nwk"),
          }],
        }],
        configs: vec![],
        pipelines: vec![],
      },
      index
    );
  }

  #[test]
  fn test_data_index_pipeline_steps_name_the_command_of_each_step() {
    let steps = pipeline_steps(include_str!("../../../../data/zika/20/pipeline.yaml")).unwrap();
    assert_eq!(
      vec![
        IndexedStep {
          name: o!("timetree"),
          command: AppCommand::Timetree,
        },
        IndexedStep {
          name: o!("mugration"),
          command: AppCommand::Mugration,
        },
      ],
      steps
    );
  }

  #[test]
  fn test_data_index_pipeline_steps_keep_environment_variables_unresolved() {
    let text = indoc! {r#"
      vars:
        dataset: "{{ env.TREETIME_DATASET }}"
      steps:
        - name: "clock"
          clock:
            tree: "{{ vars.dataset }}/tree.nwk"
    "#};
    assert_eq!(
      vec![IndexedStep {
        name: o!("clock"),
        command: AppCommand::Clock,
      }],
      pipeline_steps(text).unwrap()
    );
  }

  #[test]
  fn test_data_index_pipeline_steps_reject_a_duplicate_step_name() {
    let text = indoc! {r#"
      steps:
        - name: "a"
          clock:
            tree: "tree.nwk"
        - name: "a"
          prune:
            tree: "tree.nwk"
    "#};
    assert_error!(
      pipeline_steps(text),
      "invalid pipeline configuration: duplicate step name `a`; step names must be unique"
    );
  }
}
