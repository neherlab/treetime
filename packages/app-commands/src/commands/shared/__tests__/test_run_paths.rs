#[cfg(test)]
mod tests {
  use crate::command::{AppCommand, CommandArgs};
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::shared::reconstructed_fasta::ReconstructedFastaOutputArgs;
  use crate::commands::shared::run_paths::{RunPath, RunPaths, run_inputs, run_outputs};
  use crate::config::properties::PathRole;
  use app_output::output_plan::OutputSelection;
  use helpers::{ancestral_raw, error_message, run_path, run_paths, settings};
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
  use treetime_utils::{o, vec_of_owned};

  #[rustfmt::skip]
  #[rstest]
  #[case::different_files(           (PathRole::Input,         "x/tree.nwk"),      ("output_sequences",    "x/out.nwk"),                         &[],                          None)]
  #[case::same_file(                 (PathRole::Input,         "x/tree.nwk"),      ("output_sequences",    "x/tree.nwk"),                        &[],                          Some(("x/tree.nwk",                   "x/tree.nwk")))]
  #[case::same_file_spelled(         (PathRole::Input,         "./x//tree.nwk"),   ("output_sequences",    "x/tree.nwk"),                        &[],                          Some(("./x//tree.nwk",                "x/tree.nwk")))]
  #[case::relative_to_absolute(      (PathRole::Input,         "x/tree.nwk"),      ("output_sequences",    "CWD/x/tree.nwk"),                    &[],                          Some(("x/tree.nwk",                   "CWD/x/tree.nwk")))]
  #[case::input_template(            (PathRole::InputTemplate, "aa/{cds}.fasta"),  ("output_sequences",    "aa/S.fasta"),                        &["S"],                       Some(("aa/S.fasta",                   "aa/S.fasta")))]
  #[case::input_template_other_cds(  (PathRole::InputTemplate, "aa/{cds}.fasta"),  ("output_sequences",    "aa/S.fasta"),                        &["E"],                       None)]
  #[case::augur_placeholder(         (PathRole::InputTemplate, "aa/%GENE.fasta"),  ("output_sequences",    "aa/S.fasta"),                        &["E", "S"],                  Some(("aa/S.fasta",                   "aa/S.fasta")))]
  #[case::output_template(           (PathRole::Input,         "aa/S.fasta"),      ("output_translations", "aa/{cds}.fasta"),                    &["S"],                       Some(("aa/S.fasta",                   "aa/S.fasta")))]
  #[case::output_template_other_cds( (PathRole::Input,         "aa/S.fasta"),      ("output_translations", "aa/{cds}.fasta"),                    &["E"],                       None)]
  #[case::same_templates(            (PathRole::InputTemplate, "aa/{cds}.fasta"),  ("output_translations", "aa/%GENE.fasta"),                    &["S"],                       Some(("aa/S.fasta",                   "aa/S.fasta")))]
  #[case::both_templates(            (PathRole::InputTemplate, "aa/{cds}.fasta"),  ("output_translations", "aa/S.{cds}.fasta"),                  &["E", "S"],                  None)]
  #[case::both_templates_meet(       (PathRole::InputTemplate, "aa/{cds}.fasta"),  ("output_translations", "aa/S.{cds}.fasta"),                  &["E", "S.E"],                Some(("aa/S.E.fasta",                 "aa/S.E.fasta")))]
  #[case::augur_templates(           (PathRole::InputTemplate, "aa/%GENE.fasta"),  ("output_translations", "aa/%GENE_withInternalNodes.fasta"),  &["G", "N"],                  None)]
  #[case::augur_templates_meet(      (PathRole::InputTemplate, "aa/%GENE.fasta"),  ("output_translations", "aa/%GENE_withInternalNodes.fasta"),  &["X", "X_withInternalNodes"], Some(("aa/X_withInternalNodes.fasta", "aa/X_withInternalNodes.fasta")))]
  #[case::other_folder(              (PathRole::InputTemplate, "aa/{cds}.fasta"),  ("output_sequences",    "aa/sub/S.fasta"),                    &["S"],                       None)]
  #[case::parent_folder(             (PathRole::InputTemplate, "aa/{cds}.fasta"),  ("output_sequences",    "aa/sub/../S.fasta"),                 &["S"],                       Some(("aa/S.fasta",                   "aa/sub/../S.fasta")))]
  #[case::parent_in_input(           (PathRole::Input,         "aa/x/../S.fasta"), ("output_sequences",    "aa/S.fasta"),                        &[],                          Some(("aa/x/../S.fasta",              "aa/S.fasta")))]
  #[case::other_extension(           (PathRole::InputTemplate, "aa/{cds}.fasta"),  ("output_sequences",    "aa/S.fasta.xz"),                     &["S"],                       None)]
  #[case::placeholder_in_fixed_path( (PathRole::Input,         "aa/S.fasta"),      ("output_sequences",    "aa/{cds}.fasta"),                    &["S"],                       None)]
  #[case::literal_braces(            (PathRole::Input,         "aa/{cds}.fasta"),  ("output_sequences",    "aa/S.fasta"),                        &["S"],                       None)]
  #[case::standard_streams(          (PathRole::Input,         "-"),               ("output_sequences",    "-"),                                 &[],                          None)]
  #[case::standard_devices(          (PathRole::Input,         "/dev/stdin"),      ("output_sequences",    "/dev/stdout"),                       &[],                          None)]
  #[trace]
  fn test_run_paths_reject_an_output_that_overwrites_an_input_of_the_cdses_of_the_run(
    #[case] input: (PathRole, &str),
    #[case] output: (&str, &str),
    #[case] cdses: &[&str],
    #[case] expected: Option<(&str, &str)>,
  ) {
    let cwd = env::current_dir().unwrap();
    let cwd = cwd.to_str().unwrap();
    let (output_setting, output) = output;
    let output = output.replace("CWD", cwd);
    let paths = run_paths(input, (output_setting, &output));
    let expected = expected.map(|(input, output)| {
      format!(
        "--translations '{input}' names the same file as the output --{} '{}': the run would overwrite its own input",
        output_setting.replace('_', "-"),
        output.replace("CWD", cwd)
      )
    });
    let cdses = cdses.iter().map(|cds| (*cds).to_owned()).collect_vec();
    assert_eq!(expected, error_message(paths.check_with_cds_names(&cdses)));
  }

  #[test]
  fn test_run_paths_compare_cds_templates_only_with_cds_names() {
    let paths = RunPaths {
      inputs: vec![run_path("translations", "aa/{cds}.fasta", PathRole::InputTemplate)],
      outputs: vec![run_path("output_sequences", "aa/S.fasta", PathRole::Output)],
      cds_outputs: vec![run_path("output_translations", "aa/{cds}.fasta", PathRole::Output)],
    };
    assert_eq!(
      (
        None,
        Some(o!(
          "--translations 'aa/S.fasta' names the same file as the output --output-sequences 'aa/S.fasta': the run \
           would overwrite its own input"
        ))
      ),
      (
        error_message(paths.check_without_cds_names()),
        error_message(paths.check_with_cds_names(&vec_of_owned!["S"]))
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::other_names(   &["G", "N"],                   None)]
  #[case::suffixed_name( &["X", "X_withInternalNodes"], Some("--translations 'aa/X_withInternalNodes.fasta' names the same file as the output --output-translations 'aa/X_withInternalNodes.fasta': the run would overwrite its own input"))]
  #[trace]
  fn test_run_paths_accept_augur_translation_templates_in_one_folder(
    #[case] cdses: &[&str],
    #[case] expected: Option<&str>,
  ) {
    let raw = TreetimeAncestralArgsRaw {
      translations: Some(o!("aa/%GENE.fasta")),
      cds_selection: cdses.iter().map(|cds| (*cds).to_owned()).collect(),
      output_translations: Some(PathBuf::from("aa/%GENE_withInternalNodes.fasta")),
      ..ancestral_raw()
    };
    assert_eq!(
      expected.map(ToOwned::to_owned),
      error_message(CommandArgs::try_from(raw).map(|_| ()))
    );
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
    let paths = RunPaths {
      inputs: inputs
        .iter()
        .map(|(setting, path)| run_path(setting, path, PathRole::Input))
        .collect_vec(),
      outputs: vec![],
      cds_outputs: vec![],
    };
    assert_eq!(
      expected.map(ToOwned::to_owned),
      error_message(paths.check_without_cds_names())
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
    let template = real.join("{cds}.fasta").to_str().unwrap().to_owned();
    let check = |input: (PathRole, &str)| {
      error_message(run_paths(input, ("output_sequences", &output)).check_with_cds_names(&vec_of_owned!["S"]))
    };
    let expected = Some(format!(
      "--translations '{opened_file}' names the same file as the output --output-sequences '{output}': the run would \
       overwrite its own input"
    ));
    assert_eq!(
      (None, expected.clone(), expected),
      (
        check((PathRole::Input, &lexical_twin)),
        check((PathRole::Input, &opened_file)),
        check((PathRole::InputTemplate, &template)),
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
    let raw = TreetimeAncestralArgsRaw {
      reconstructed_fasta_output: ReconstructedFastaOutputArgs {
        output_sequences: Some(PathBuf::from("seq.fasta")),
      },
      ..ancestral_raw()
    };
    let args = TreetimeAncestralArgs::try_from(raw).unwrap();
    let by_selection = btreemap! {
      OutputSelection::Sequences => vec![PathBuf::from("seq.fasta")],
      OutputSelection::Model => vec![PathBuf::from("out/ancestral.model.json")],
    };
    assert_eq!(
      vec![
        run_path("output_all", "out/ancestral.model.json", PathRole::Output),
        run_path("output_sequences", "seq.fasta", PathRole::Output),
      ],
      run_outputs(&args.outputs, &by_selection)
        .into_iter()
        .sorted_by(|a, b| a.path.cmp(&b.path))
        .collect_vec()
    );
  }

  mod helpers {
    use super::*;
    use crate::commands::shared::alignment::AlignmentArgs;
    use crate::commands::shared::output_args::OutputCoreArgs;
    use crate::commands::shared::tree_input::TreeArgs;
    use eyre::Report;
    use treetime_utils::error::report_to_string;

    pub(super) fn run_path(setting: &str, path: &str, role: PathRole) -> RunPath {
      RunPath {
        setting: setting.to_owned(),
        path: path.to_owned(),
        role,
      }
    }

    pub(super) fn run_paths((role, input): (PathRole, &str), (output_setting, output): (&str, &str)) -> RunPaths {
      let output = run_path(output_setting, output, PathRole::Output);
      let (outputs, cds_outputs) = if output_setting == "output_translations" {
        (vec![], vec![output])
      } else {
        (vec![output], vec![])
      };
      RunPaths {
        inputs: vec![run_path("translations", input, role)],
        outputs,
        cds_outputs,
      }
    }

    pub(super) fn ancestral_raw() -> TreetimeAncestralArgsRaw {
      TreetimeAncestralArgsRaw {
        tree: TreeArgs {
          tree: Some(PathBuf::from("tree.nwk")),
        },
        alignment: AlignmentArgs {
          alignment: vec![PathBuf::from("aln.fasta")],
        },
        output: OutputCoreArgs {
          output_all: Some(PathBuf::from("out")),
          ..OutputCoreArgs::default()
        },
        ..TreetimeAncestralArgsRaw::default()
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
