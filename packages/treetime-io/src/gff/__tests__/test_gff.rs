#[cfg(test)]
mod tests {
  use crate::gff::*;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::o;

  #[test]
  fn test_parse_gff3_cds_features_basic() {
    let input = "\
##gff-version 3
MN908947.3\tNextclade\tCDS\t21563\t25384\t.\t+\t0\tID=cds-S;Name=S
MN908947.3\tNextclade\tCDS\t266\t13468\t.\t+\t0\tID=cds-ORF1a;Name=ORF%31a
MN908947.3\tNextclade\tCDS\t13468\t21555\t.\t+\t0\tID=cds-ORF1a-2;Name=ORF%31a
MN908947.3\tNextclade\tgene\t100\t200\t.\t+\t.\tName=ignored
";
    let features = parse_gff3_cds_features(input, Path::new("test.gff3")).unwrap();

    assert_eq!(2, features.len());
    assert_eq!("ORF1a", features[0].name);
    assert_eq!(2, features[0].segments.len());
    assert_eq!("S", features[1].name);
    assert_eq!(1, features[1].segments.len());
    assert_eq!(21563, features[1].segments[0].start);
  }

  #[test]
  fn test_parse_gff3_cds_features_minus_strand_5_to_3() {
    let input = "\
##gff-version 3
seq1\tNextclade\tCDS\t100\t150\t.\t-\t0\tName=ORF
seq1\tNextclade\tCDS\t200\t250\t.\t-\t0\tName=ORF
";
    let features = parse_gff3_cds_features(input, Path::new("test.gff3")).unwrap();

    assert_eq!(1, features.len());
    assert_eq!(200, features[0].segments[0].start);
    assert_eq!(100, features[0].segments[1].start);
  }

  #[test]
  fn test_parse_gff3_cds_features_rejects_non_multiple_of_three() {
    let input = "\
##gff-version 3
seq1\tNextclade\tCDS\t1\t5\t.\t+\t0\tName=BAD
";
    let err = parse_gff3_cds_features(input, Path::new("test.gff3")).unwrap_err();
    assert!(err.to_string().contains("not a multiple of 3"));
  }

  #[test]
  fn test_parse_gff3_cds_features_stops_at_fasta_directive() {
    let input = "\
##gff-version 3
seq1\tNextclade\tCDS\t1\t6\t.\t+\t0\tName=A
##FASTA
>seq1
ACGTAC
";
    let features = parse_gff3_cds_features(input, Path::new("test.gff3")).unwrap();
    assert_eq!(1, features.len());
    assert_eq!("A", features[0].name);
  }

  #[test]
  fn test_parse_gff3_cds_features_rejects_mixed_seqids() {
    let input = "\
##gff-version 3
chr1\tNextclade\tCDS\t1\t3\t.\t+\t0\tName=X
chr2\tNextclade\tCDS\t4\t6\t.\t+\t0\tName=X
";
    let err = parse_gff3_cds_features(input, Path::new("test.gff3")).unwrap_err();
    assert!(err.to_string().contains("multiple sequence IDs"));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::name(        btreemap! { o!("Name") => o!("S") },                          Some(o!("S")))]
  #[case::locus_tag(   btreemap! { o!("locus_tag") => o!("LT1"), o!("Parent") => o!("p") }, Some(o!("LT1")))]
  #[case::gene_over_id(btreemap! { o!("gene") => o!("G"), o!("ID") => o!("cds-G") },  Some(o!("G")))]
  #[case::parent_only( btreemap! { o!("Parent") => o!("gene-Y") },                    None)]
  #[case::id_fallback( btreemap! { o!("ID") => o!("cds-X") },                         Some(o!("cds-X")))]
  fn test_resolve_cds_name_priority(#[case] attrs: BTreeMap<String, String>, #[case] expected: Option<String>) {
    assert_eq!(expected, resolve_cds_name(&attrs));
  }
}
