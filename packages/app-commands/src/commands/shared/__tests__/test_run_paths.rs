#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use crate::commands::shared::run_paths::{RunPath, check_run_paths, run_inputs, run_outputs};
  use crate::config::properties::PathRole;
  use app_output::output_plan::OutputSelection;
  use helpers::{error_message, run_path, settings};
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Map, Value, json};
  use std::env;
  use std::fs;
  use std::os::unix::fs::symlink;
  use std::path::PathBuf;
  use tempfile::tempdir;

  #[rustfmt::skip]
  #[rstest]
  #[case::different_files(      (PathRole::Input,         "x/tree.nwk"),      "x/out.nwk",            None)]
  #[case::same_file(            (PathRole::Input,         "x/tree.nwk"),      "x/tree.nwk",           Some("names"))]
  #[case::same_file_spelled(    (PathRole::Input,         "./x//tree.nwk"),   "x/tree.nwk",           Some("names"))]
  #[case::relative_to_absolute( (PathRole::Input,         "x/tree.nwk"),      "CWD/x/tree.nwk"  ,     Some("names"))]
  #[case::input_template(       (PathRole::InputTemplate, "aa/{cds}.fasta"),  "aa/S.fasta",           Some("can name"))]
  #[case::augur_placeholder(    (PathRole::InputTemplate, "aa/%GENE.fasta"),  "aa/S.fasta",           Some("can name"))]
  #[case::output_template(      (PathRole::Input,         "aa/S.fasta"),      "aa/{cds}.fasta",       Some("can name"))]
  #[case::both_templates(       (PathRole::InputTemplate, "aa/{cds}.fasta"),  "aa/S.{cds}.fasta",     Some("can name"))]
  #[case::empty_placeholder(    (PathRole::InputTemplate, "aa/x{cds}.fasta"), "aa/x.fasta",           Some("can name"))]
  #[case::other_folder(         (PathRole::InputTemplate, "aa/{cds}.fasta"),  "aa/sub/S.fasta",       None)]
  #[case::parent_folder(        (PathRole::InputTemplate, "aa/{cds}.fasta"),  "aa/sub/../S.fasta",    Some("can name"))]
  #[case::parent_in_input(      (PathRole::Input,         "aa/x/../S.fasta"), "aa/S.fasta",           Some("names"))]
  #[case::other_extension(      (PathRole::InputTemplate, "aa/{cds}.fasta"),  "aa/S.fasta.xz",        None)]
  #[case::alignment_output(     (PathRole::Input,         "seq.2.fasta"),     "seq.{alignment}.fasta", Some("can name"))]
  #[case::literal_braces(       (PathRole::Input,         "aa/{cds}.fasta"),  "aa/S.fasta",           None)]
  #[case::standard_streams(     (PathRole::Input,         "-"),               "-",                    None)]
  #[case::standard_devices(     (PathRole::Input,         "/dev/stdin"),      "/dev/stdout",          None)]
  #[trace]
  fn test_run_paths_reject_an_output_that_can_overwrite_an_input(
    #[case] (role, input): (PathRole, &str),
    #[case] output: &str,
    #[case] expected_verb: Option<&str>,
  ) {
    let cwd = env::current_dir().unwrap();
    let output = output.replace("CWD", cwd.to_str().unwrap());
    let inputs = [run_path("translations", input, role)];
    let outputs = [run_path("output_sequences", &output, PathRole::Output)];
    let expected = expected_verb.map(|verb| {
      format!(
        "--translations '{input}' {verb} the same file as the output --output-sequences '{output}': the run would \
         overwrite its own input"
      )
    });
    assert_eq!(expected, error_message(check_run_paths(&inputs, &outputs)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::once(            &[("tree", "-")],                                         None)]
  #[case::once_device(     &[("tree", "/dev/stdin")],                                None)]
  #[case::same_setting(    &[("alignment", "-"), ("alignment", "-")],                Some("Standard input can be read only once, but --alignment lists it 2 times"))]
  #[case::same_device(     &[("alignment", "/dev/stdin"), ("alignment", "-")],       Some("Standard input can be read only once, but --alignment lists it 2 times"))]
  #[case::two_settings(    &[("tree", "-"), ("remove_file", "-")],                   Some("Standard input can be read only once, but --tree and --remove-file read it"))]
  #[case::mixed_spellings( &[("tree", "-"), ("alignment", "/dev/stdin")],            Some("Standard input can be read only once, but --tree and --alignment read it"))]
  #[case::three(           &[("tree", "-"), ("alignment", "-"), ("drms", "-")],      Some("Standard input can be read only once, but --tree, --alignment and --drms read it"))]
  #[trace]
  fn test_run_paths_read_standard_input_once(#[case] inputs: &[(&str, &str)], #[case] expected: Option<&str>) {
    let inputs = inputs
      .iter()
      .map(|(setting, path)| run_path(setting, path, PathRole::Input))
      .collect_vec();
    assert_eq!(
      expected.map(ToOwned::to_owned),
      error_message(check_run_paths(&inputs, &[]))
    );
  }

  #[test]
  fn test_run_paths_compare_the_files_behind_symbolic_links() {
    let dir = tempdir().unwrap();
    let real = dir.path().join("real");
    fs::create_dir_all(real.join("sub")).unwrap();
    symlink(real.join("sub"), dir.path().join("link")).unwrap();
    let output = dir.path().join("link/../S.fasta").to_str().unwrap().to_owned();
    let lexical_twin = dir.path().join("S.fasta").to_str().unwrap().to_owned();
    let opened_file = real.join("S.fasta").to_str().unwrap().to_owned();
    let outputs = [run_path("output_sequences", &output, PathRole::Output)];
    assert_eq!(
      (
        None,
        Some(format!(
          "--translations '{opened_file}' names the same file as the output --output-sequences '{output}': the run \
           would overwrite its own input"
        ))
      ),
      (
        error_message(check_run_paths(
          &[run_path("translations", &lexical_twin, PathRole::Input)],
          &outputs
        )),
        error_message(check_run_paths(
          &[run_path("translations", &opened_file, PathRole::Input)],
          &outputs
        ))
      )
    );
  }

  #[test]
  fn test_run_paths_inputs_follow_the_path_roles_of_the_schema() {
    let settings = settings(json!({
      "tree": "t.nwk",
      "alignment": ["a.fasta", "b.fasta"],
      "translations": "aa/{cds}.fasta",
      "annotation": "genes.gff",
      "aa_root_sequence": "root.fasta",
      "topology_order_list": "order.txt",
      "output_sequences": "out.fasta"
    }));
    assert_eq!(
      vec![
        run_path("aa_root_sequence", "root.fasta", PathRole::Input),
        run_path("alignment", "a.fasta", PathRole::Input),
        run_path("alignment", "b.fasta", PathRole::Input),
        run_path("annotation", "genes.gff", PathRole::Input),
        run_path("topology_order_list", "order.txt", PathRole::Input),
        run_path("translations", "aa/{cds}.fasta", PathRole::InputTemplate),
        run_path("tree", "t.nwk", PathRole::Input),
      ],
      run_inputs(AppCommand::Ancestral, &settings)
        .unwrap()
        .into_iter()
        .sorted_by(|a, b| (&a.setting, &a.path).cmp(&(&b.setting, &b.path)))
        .collect_vec()
    );
  }

  #[test]
  fn test_run_paths_outputs_name_the_setting_that_planned_them() {
    let settings = settings(json!({ "output_all": "out", "output_sequences": "seq.fasta" }));
    let by_selection = btreemap! {
      OutputSelection::Sequences => vec![PathBuf::from("seq.fasta")],
      OutputSelection::Model => vec![PathBuf::from("out/ancestral.model.json")],
    };
    assert_eq!(
      vec![
        run_path("output_sequences", "seq.fasta", PathRole::Output),
        run_path("output_all", "out/ancestral.model.json", PathRole::Output),
      ]
      .into_iter()
      .sorted_by(|a, b| a.path.cmp(&b.path))
      .collect_vec(),
      run_outputs(&settings, &by_selection)
        .into_iter()
        .sorted_by(|a, b| a.path.cmp(&b.path))
        .collect_vec()
    );
  }

  mod helpers {
    use super::*;
    use eyre::Report;
    use treetime_utils::error::report_to_string;

    pub(super) fn run_path(setting: &str, path: &str, role: PathRole) -> RunPath {
      RunPath {
        setting: setting.to_owned(),
        path: path.to_owned(),
        role,
      }
    }

    pub(super) fn settings(value: Value) -> Map<String, Value> {
      let Value::Object(settings) = value else {
        panic!("settings must be a JSON object")
      };
      settings
    }

    pub(super) fn error_message(result: Result<(), Report>) -> Option<String> {
      result.err().map(|report| report_to_string(&report))
    }
  }
}
