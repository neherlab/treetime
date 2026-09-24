#[cfg(test)]
mod tests {
  use crate::commands::ancestral::aa_node_data::{translation_path, validate_aa_args, validate_aa_root_sequence_cdses};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::path::{Path, PathBuf};
  use treetime_primitives::Seq;
  use treetime_utils::o;

  const CDS_TEMPLATE: &str = "out/{cds}.fasta";

  #[test]
  fn test_validate_aa_args_requires_cds_placeholder() {
    let err = validate_aa_args(Some("translations.fasta"), &["S".to_owned()], None, None).unwrap_err();
    assert!(err.to_string().contains("CDS placeholder"));
  }

  #[test]
  fn test_validate_aa_args_accepts_percent_gene_placeholder() {
    let result = validate_aa_args(Some("out/%GENE.fasta"), &["S".to_owned()], None, None);
    result.unwrap();
  }

  #[test]
  fn test_validate_aa_args_accepts_cds_placeholder() {
    let result = validate_aa_args(
      Some(CDS_TEMPLATE),
      &["S".to_owned()],
      None,
      None,
    );
    result.unwrap();
  }

  #[test]
  fn test_validate_aa_args_empty_cdses_with_annotation_ok() {
    let result = validate_aa_args(
      Some(CDS_TEMPLATE),
      &[],
      Some(Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/Cargo.toml"))),
      None,
    );
    result.unwrap();
  }

  #[test]
  fn test_validate_aa_args_empty_cdses_no_annotation_errors() {
    let err = validate_aa_args(Some(CDS_TEMPLATE), &[], None, None).unwrap_err();
    assert!(err.to_string().contains("--cdses"));
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

    let err = validate_aa_root_sequence_cdses(Path::new("roots.fasta"), &by_cds, &[o!("S"), o!("M")]).unwrap_err();

    assert!(err.to_string().contains("CDS 'M'"));
  }
}
