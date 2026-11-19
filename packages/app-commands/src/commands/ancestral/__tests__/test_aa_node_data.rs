#[cfg(test)]
mod tests {
  use crate::commands::ancestral::aa_node_data::{
    cds_names, cds_output_paths, translation_path, validate_aa_args, validate_aa_root_sequence_cdses,
  };
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use std::fs;
  use std::path::{Path, PathBuf};
  use tempfile::tempdir;
  use treetime_primitives::Seq;
  use treetime_utils::{assert_error, o, vec_of_owned};

  const CDS_TEMPLATE: &str = "out/{cds}.fasta";

  #[rustfmt::skip]
  #[rstest]
  #[case::no_aa_settings(        (None,                    vec![],             None))]
  #[case::gene_placeholder(      (Some("out/%GENE.fasta"), vec_of_owned!["S"], None))]
  #[case::cds_placeholder(       (Some(CDS_TEMPLATE),      vec_of_owned!["S"], None))]
  #[case::cdses_from_annotation( (Some(CDS_TEMPLATE),      vec![],             Some(Path::new("genes.gff"))))]
  #[trace]
  fn test_validate_aa_args_accepts(
    #[case] (translations, cdses, annotation): (Option<&str>, Vec<String>, Option<&Path>),
  ) -> Result<(), Report> {
    validate_aa_args(translations, &cdses, annotation, None, 1, None)
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::cdses_without_translations(       (None,                       vec_of_owned!["S"],      1), "--translations is required when using --cds-selection, --annotation, or --aa-root-sequence")]
  #[case::two_alignments(                   (Some(CDS_TEMPLATE),         vec_of_owned!["S"],      2), "--translations and --annotation need exactly one --alignment file, but 2 were given: CDS coordinates refer to the positions of one nucleotide alignment")]
  #[case::no_placeholder(                   (Some("translations.fasta"), vec_of_owned!["S"],      1), "--translations must contain a CDS placeholder ('{cds}' or '%GENE')")]
  #[case::empty_selection_without_annotation((Some(CDS_TEMPLATE),        vec![],                  1), "--cds-selection must list at least one CDS, or pass --annotation to derive the CDS set")]
  #[case::empty_cds_name(                   (Some(CDS_TEMPLATE),         vec_of_owned!["S", ""],  1), "--cds-selection must not contain an empty CDS name")]
  #[case::duplicate_cds(                    (Some(CDS_TEMPLATE),         vec_of_owned!["S", "S"], 1), "--cds-selection contains duplicate CDS 'S'")]
  #[trace]
  fn test_validate_aa_args_rejects(
    #[case] (translations, cdses, n_alignments): (Option<&str>, Vec<String>, usize),
    #[case] expected: &str,
  ) {
    assert_error!(
      validate_aa_args(translations, &cdses, None, None, n_alignments, None),
      expected
    );
  }

  #[test]
  fn test_validate_aa_args_rejects_output_translations_without_translations() {
    assert_error!(
      validate_aa_args(None, &[], None, None, 1, Some(Path::new("aa.fasta"))),
      "--output-translations requires --translations, which names the amino-acid alignments to reconstruct"
    );
  }

  #[expect(clippy::literal_string_with_formatting_args, reason = "the braces are a path template placeholder, not a format argument")]
  #[rustfmt::skip]
  #[rstest]
  #[case::cds_placeholder( "out/{cds}.fasta",  "S",  "out/S.fasta")]
  #[case::gene_placeholder("out/%GENE.fasta",   "S",  "out/S.fasta")]
  #[expect(clippy::literal_string_with_formatting_args, reason = "the braces are a path template placeholder, not a format argument")]
  #[case::both_placeholders("out/{cds}/%GENE.fasta", "ORF1a", "out/ORF1a/ORF1a.fasta")]
  fn test_translation_path_expands_placeholders(
    #[case] template: &str,
    #[case] cds: &str,
    #[case] expected: &str,
  ) {
    assert_eq!(PathBuf::from(expected), translation_path(template, cds));
  }

  #[test]
  fn test_validate_aa_root_sequence_cdses_requires_every_cds() {
    let by_cds = btreemap! {
      o!("S") => Seq::try_from_str("AC").unwrap(),
    };

    assert_error!(
      validate_aa_root_sequence_cdses(Path::new("roots.fasta"), &by_cds, &[o!("S"), o!("M")]),
      "--aa-root-sequence 'roots.fasta' does not contain a FASTA record for CDS 'M'"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::cds_placeholder(       "out/{cds}.fasta",         "out/S.fasta")]
  #[case::gene_placeholder(      "out/%GENE.fasta",         "out/S.fasta")]
  #[case::folder_keeps_template( "out/{cds}/%GENE.fasta",   "out/{cds}/S.fasta")]
  #[case::folder_only(           "out/{cds}/aa.fasta",      "out/{cds}/aa.fasta")]
  #[case::no_placeholder(        "out/aa.fasta",            "out/aa.fasta")]
  #[trace]
  fn test_cds_output_paths_expand_placeholders_in_file_name_only(
    #[case] template: &str,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let actual = cds_output_paths(Path::new(template), &[o!("S")])?;

    assert_eq!(btreemap! { o!("S") => PathBuf::from(expected) }, actual);
    Ok(())
  }

  #[test]
  fn test_cds_output_paths_one_path_per_cds() -> Result<(), Report> {
    let actual = cds_output_paths(Path::new("out/aa.{cds}.fasta"), &vec_of_owned!["S", "ORF1a"])?;

    let expected = btreemap! {
      o!("ORF1a") => PathBuf::from("out/aa.ORF1a.fasta"),
      o!("S") => PathBuf::from("out/aa.S.fasta"),
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_placeholder(    "out/aa.fasta")]
  #[case::folder_placeholder("out/{cds}/aa.fasta")]
  #[trace]
  fn test_cds_output_paths_two_cdses_to_one_path_fail(#[case] template: &str) {
    assert_error!(
      cds_output_paths(Path::new(template), &vec_of_owned!["S", "ORF1a"]),
      "--output-translations needs the placeholder '{cds}' or '%GENE' in its file name when \
       reconstructing several CDSes, otherwise the output of each CDS overwrites the same file"
    );
  }

  #[test]
  fn test_cds_output_paths_cds_with_path_separator_fails() {
    assert_error!(
      cds_output_paths(Path::new("out/aa.{cds}.fasta"), &vec_of_owned!["S", "../S"]),
      "CDS '../S' contains a path separator and cannot be part of the file name of --output-translations"
    );
  }

  #[test]
  fn test_cds_output_paths_cds_with_path_separator_kept_when_not_in_file_name() -> Result<(), Report> {
    let actual = cds_output_paths(Path::new("out/aa.fasta"), &[o!("a/b")])?;

    assert_eq!(btreemap! { o!("a/b") => PathBuf::from("out/aa.fasta") }, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::listed_names_first(  json!({ "cds_selection": ["G", "F"], "annotation": "genes.gff" }), vec!["G", "F"])]
  #[case::annotation_names(    json!({ "annotation": "genes.gff" }),                             vec!["F", "G"])]
  #[case::empty_selection(     json!({ "cds_selection": [], "annotation": "genes.gff" }),        vec!["F", "G"])]
  #[case::neither(             json!({}),                                                       vec![])]
  #[trace]
  fn test_cds_names_of_settings(#[case] settings: Value, #[case] expected: Vec<&str>) -> Result<(), Report> {
    let dir = tempdir()?;
    let gff = dir.path().join("genes.gff");
    fs::write(
      &gff,
      "##gff-version 3\nref\t.\tCDS\t1\t9\t.\t+\t0\tName=G\nref\t.\tCDS\t10\t18\t.\t+\t0\tName=F\n",
    )?;
    let mut settings = settings;
    if settings.get("annotation").is_some() {
      settings["annotation"] = json!([gff]);
    }

    let actual = cds_names(settings.as_object().unwrap())?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_cds_names_listed_cds_absent_from_annotation_fails() {
    let dir = tempdir().unwrap();
    let gff = dir.path().join("genes.gff");
    fs::write(&gff, "##gff-version 3\nref\t.\tCDS\t1\t9\t.\t+\t0\tName=F\n").unwrap();
    let settings = json!({ "cds_selection": ["S"], "annotation": [gff] });

    assert_error!(
      cds_names(settings.as_object().unwrap()),
      format!(
        "--annotation '{}' does not contain a CDS feature for CDS 'S'",
        gff.display()
      )
    );
  }
}
