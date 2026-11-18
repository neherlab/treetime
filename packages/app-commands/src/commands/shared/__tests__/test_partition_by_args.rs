#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::argument_error;
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::optimize::args::{TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use clap::Parser;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  const OPTIMIZE: [&str; 5] = ["optimize", "--tree", "t.nwk", "--alignment", "a.fasta"];
  const ANCESTRAL: [&str; 5] = ["ancestral", "--tree", "t.nwk", "--alignment", "a.fasta"];

  #[rustfmt::skip]
  #[rstest]
  #[case::scheme_without_annotation(&["--partition-by", "codon-positions"],                                "--partition-by codon-positions needs one --annotation file per --alignment file, matched in list order, but the run has 1 --alignment file and 0 --annotation files", &["annotation", "partition_by", "alignment"])]
  #[case::two_annotations(          &["--partition-by", "genes", "--annotation", "a.gff", "--annotation", "b.gff"], "--partition-by genes needs one --annotation file per --alignment file, matched in list order, but the run has 1 --alignment file and 2 --annotation files", &["annotation", "partition_by", "alignment"])]
  #[case::rate_categories_annotation(&["--partition-by", "rate-categories", "--annotation", "a.gff"],       "--annotation applies only with --partition-by codon-positions or genes, which read the CDSs of each alignment from it", &["annotation", "partition_by"])]
  #[case::annotation_alone(         &["--annotation", "a.gff"],                                             "--annotation applies only with --partition-by codon-positions or genes, which read the CDSs of each alignment from it", &["annotation", "partition_by"])]
  #[case::rule_without_scheme(      &["--rate-category-rule", "k-means"],                                   "--rate-category-rule applies only with --partition-by rate-categories", &["rate_category_rule", "partition_by"])]
  #[case::rule_with_genes(          &["--partition-by", "genes", "--annotation", "a.gff", "--rate-category-rule", "median-split"], "--rate-category-rule applies only with --partition-by rate-categories", &["rate_category_rule", "partition_by"])]
  #[case::amino_acids(              &["--partition-by", "rate-categories", "--alphabet", "aa"],             "--partition-by needs the alphabet nuc, because CDS positions are nucleotide positions and the rate categories count nucleotide changes, but the alphabet is aa", &["partition_by", "alphabet"])]
  #[trace]
  fn test_optimize_scheme_argument_errors(#[case] args: &[&str], #[case] message: &str, #[case] settings: &[&str]) {
    let raw = TreetimeOptimizeArgsRaw::try_parse_from([&OPTIMIZE[..], args].concat()).unwrap();
    let expected = (message.to_owned(), settings.iter().map(|setting| (*setting).to_owned()).collect());
    assert_eq!(Some(expected), argument_error(TreetimeOptimizeArgs::try_from(raw)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::annotation_alone(  &["--annotation", "a.gff"],                                   "--annotation applies only with --partition-by codon-positions or genes, which read the CDSs of each alignment from it, or with --translations", &["annotation", "partition_by"])]
  #[case::parsimony(         &["--method-anc", "parsimony", "--partition-by", "rate-categories"], "--partition-by applies only with --method-anc marginal: parsimony reconstruction uses no substitution model and treats every site alone, so partitions do not change its result", &["partition_by", "method_anc"])]
  #[case::translations_two(  &["--translations", "t/{cds}.fasta", "--annotation", "a.gff", "--annotation", "b.gff"], "--translations takes one --annotation file, the CDS coordinates of its one --alignment file, but 2 were given", &["annotation", "translations"])]
  #[trace]
  fn test_ancestral_scheme_argument_errors(#[case] args: &[&str], #[case] message: &str, #[case] settings: &[&str]) {
    let raw = TreetimeAncestralArgsRaw::try_parse_from([&ANCESTRAL[..], args].concat()).unwrap();
    let expected = (message.to_owned(), settings.iter().map(|setting| (*setting).to_owned()).collect());
    assert_eq!(Some(expected), argument_error(TreetimeAncestralArgs::try_from(raw)));
  }

  #[test]
  fn test_timetree_scheme_without_alignment_is_an_argument_error() {
    let raw = TreetimeTimetreeArgsRaw::try_parse_from([
      "timetree",
      "--tree",
      "t.nwk",
      "--metadata",
      "m.tsv",
      "--partition-by",
      "rate-categories",
    ])
    .unwrap();
    assert_eq!(
      Some((
        "--partition-by needs at least one --alignment file, whose sites it puts into partitions".to_owned(),
        vec!["partition_by".to_owned(), "alignment".to_owned()]
      )),
      argument_error(TreetimeTimetreeArgs::try_from(raw))
    );
  }

  #[test]
  fn test_scheme_counts_as_several_partitions_for_partition_rates() -> Result<(), Report> {
    let raw = TreetimeOptimizeArgsRaw::try_parse_from(
      [
        &OPTIMIZE[..],
        &["--partition-by", "rate-categories", "--partition-rates", "equal"],
      ]
      .concat(),
    )?;
    assert_eq!(None, argument_error(TreetimeOptimizeArgs::try_from(raw)));
    Ok(())
  }
}
