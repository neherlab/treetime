#[cfg(test)]
mod tests {
  use crate::gff::*;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_utils::assert_error;
  use treetime_utils::o;

  #[test]
  fn test_gff3_cds_read_basic() {
    let input = indoc! {"
      ##gff-version 3
      ##sequence-region MN908947.3 1 29903
      MN908947.3\tNextclade\tCDS\t21563\t25384\t.\t+\t0\tID=cds-S;Name=S
      MN908947.3\tNextclade\tCDS\t266\t13468\t.\t+\t0\tID=cds-ORF1a;Name=ORF%31a
      MN908947.3\tNextclade\tCDS\t13468\t21555\t.\t+\t2\tID=cds-ORF1a;Name=ORF%31a
      MN908947.3\tNextclade\tgene\t100\t200\t.\t.\t.\tName=ignored
    "};
    let actual = gff3_cds_read(input.as_bytes()).unwrap();

    let expected = GffCdsAnnotation {
      features: vec![
        GffCdsFeature {
          label: o!("ORF1a"),
          id: Some(o!("cds-ORF1a")),
          name: o!("ORF1a"),
          seqid: o!("MN908947.3"),
          strand: GffStrand::Forward,
          segments: vec![
            GffCdsSegment {
              start: 266,
              end: 13468,
              phase: GffPhase::Zero,
            },
            GffCdsSegment {
              start: 13468,
              end: 21555,
              phase: GffPhase::Two,
            },
          ],
        },
        GffCdsFeature {
          label: o!("S"),
          id: Some(o!("cds-S")),
          name: o!("S"),
          seqid: o!("MN908947.3"),
          strand: GffStrand::Forward,
          segments: vec![GffCdsSegment {
            start: 21563,
            end: 25384,
            phase: GffPhase::Zero,
          }],
        },
      ],
      sequence_regions: vec![GffSequenceRegion {
        seqid: o!("MN908947.3"),
        start: 1,
        end: 29903,
      }],
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_gff3_cds_read_minus_strand_5_to_3() {
    let input = indoc! {"
      ##gff-version 3
      seq1\tNextclade\tCDS\t100\t150\t.\t-\t1\tID=orf;Name=ORF
      seq1\tNextclade\tCDS\t200\t250\t.\t-\t0\tID=orf;Name=ORF
    "};
    let actual = gff3_cds_read(input.as_bytes()).unwrap().features;

    let expected = vec![GffCdsFeature {
      label: o!("ORF"),
      id: Some(o!("orf")),
      name: o!("ORF"),
      seqid: o!("seq1"),
      strand: GffStrand::Reverse,
      segments: vec![
        GffCdsSegment {
          start: 200,
          end: 250,
          phase: GffPhase::Zero,
        },
        GffCdsSegment {
          start: 100,
          end: 150,
          phase: GffPhase::One,
        },
      ],
    }];
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_gff3_cds_read_rows_without_id_are_separate_features() {
    let input = indoc! {"
      ##gff-version 3
      seq1\tsrc\tCDS\t1\t6\t.\t+\t0\tName=A
      seq1\tsrc\tCDS\t10\t15\t.\t+\t0\tName=B
    "};
    let actual = gff3_cds_read(input.as_bytes())
      .unwrap()
      .features
      .into_iter()
      .map(|feature| (feature.label, feature.id, feature.segments.len()))
      .collect::<Vec<_>>();

    let expected = vec![(o!("A"), None, 1), (o!("B"), None, 1)];
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_gff3_cds_read_same_name_different_ids_stay_two_cdses() {
    let input = indoc! {"
      ##gff-version 3
      seq1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds-1;Name=ORF
      seq1\tsrc\tCDS\t10\t15\t.\t+\t0\tID=cds-2;Name=ORF
    "};
    let actual = gff3_cds_read(input.as_bytes())
      .unwrap()
      .features
      .into_iter()
      .map(|feature| (feature.label, feature.name))
      .collect::<Vec<_>>();

    let expected = vec![(o!("cds-1"), o!("ORF")), (o!("cds-2"), o!("ORF"))];
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_gff3_cds_read_rejects_repeated_name_without_id() {
    let input = indoc! {"
      ##gff-version 3
      seq1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds-1;Name=ORF
      seq1\tsrc\tCDS\t10\t15\t.\t+\t0\tName=ORF
    "};
    assert_error!(
      gff3_cds_read(input.as_bytes()),
      "When reading GFF3: CDS name 'ORF' is used by several CDS features (GFF3 lines 2, 3), and the CDS on GFF3 line 3 has no ID to tell them apart"
    );
  }

  #[test]
  fn test_gff3_cds_read_rejects_non_multiple_of_three() {
    let input = indoc! {"
      ##gff-version 3
      seq1\tNextclade\tCDS\t1\t5\t.\t+\t0\tName=BAD
    "};
    assert_error!(
      gff3_cds_read(input.as_bytes()),
      "When reading GFF3: CDS feature for CDS 'BAD' has nucleotide length 5, which is not a multiple of 3. Check that the annotation matches the translations."
    );
  }

  #[test]
  fn test_gff3_cds_read_stops_at_fasta_directive() {
    let input = indoc! {"
      ##gff-version 3
      seq1\tNextclade\tCDS\t1\t6\t.\t+\t0\tName=A
      ##FASTA
      >seq1
      ACGTAC
    "};
    let features = gff3_cds_read(input.as_bytes()).unwrap().features;
    assert_eq!(
      vec![o!("A")],
      features.into_iter().map(|feature| feature.label).collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_gff3_cds_read_rejects_mixed_seqids() {
    let input = indoc! {"
      ##gff-version 3
      chr1\tNextclade\tCDS\t1\t3\t.\t+\t0\tID=x;Name=X
      chr2\tNextclade\tCDS\t4\t6\t.\t+\t0\tID=x;Name=X
    "};
    assert_error!(
      gff3_cds_read(input.as_bytes()),
      "When reading GFF3: CDS feature rows for CDS 'X' span multiple sequence IDs"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::unstranded(".", "When reading GFF3: Invalid CDS strand '.' for CDS 'NS1' on GFF3 line 5: a CDS needs strand '+' or '-'")]
  #[case::unknown(   "?", "When reading GFF3: Invalid CDS strand '?' for CDS 'NS1' on GFF3 line 5: a CDS needs strand '+' or '-'")]
  #[trace]
  fn test_gff3_cds_read_rejects_strand(#[case] strand: &str, #[case] expected: &str) {
    let input = format!("##gff-version 3\n#\n#\n#\nseq1\tsrc\tCDS\t1\t6\t.\t{strand}\t0\tName=NS1\n");
    assert_error!(gff3_cds_read(input.as_bytes()), expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::missing(".", "When reading GFF3: Invalid CDS phase '.' for CDS 'A' on GFF3 line 2: a CDS needs phase 0, 1 or 2")]
  #[case::three(  "3", "When reading GFF3: Invalid CDS phase '3' for CDS 'A' on GFF3 line 2: a CDS needs phase 0, 1 or 2")]
  #[trace]
  fn test_gff3_cds_read_rejects_phase(#[case] phase: &str, #[case] expected: &str) {
    let input = format!("##gff-version 3\nseq1\tsrc\tCDS\t1\t6\t.\t+\t{phase}\tName=A\n");
    assert_error!(gff3_cds_read(input.as_bytes()), expected);
  }

  #[test]
  fn test_gff3_cds_read_rejects_file_without_cds_rows() {
    let input = indoc! {"
      ##gff-version 3
      seq1\tsrc\tgene\t1\t6\t.\t+\t.\tName=A
    "};
    assert_error!(
      gff3_cds_read(input.as_bytes()),
      "When reading GFF3: GFF3 file has no CDS rows: TreeTime reads the coding regions from CDS rows, and ignores gene and other feature rows"
    );
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
