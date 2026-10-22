#[cfg(test)]
mod tests {
  use crate::data_index::pipeline_steps;
  use app_backend::datasets::IndexedStep;
  use app_commands::command::AppCommand;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use treetime_utils::{assert_error, o};

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
