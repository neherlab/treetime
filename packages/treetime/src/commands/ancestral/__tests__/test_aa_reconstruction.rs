#[cfg(test)]
mod tests {
  use crate::ancestral::params::MethodAncestral;
  use crate::commands::ancestral::args::TreetimeAncestralArgs;
  use crate::commands::ancestral::run::run_ancestral_reconstruction;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::output::OutputArgs;
  use crate::progress::NoopProgress;
  use std::fs;
  use tempfile::tempdir;
  use treetime_utils::io::json::json_read_str;
  use util_augur_node_data_json::AugurNodeDataJsonAncestral;

  /// End-to-end amino-acid reconstruction with the default `infer` model over the stop-inclusive
  /// alphabet. The peptides contain a stop codon `*`, which must be carried as a real state (not
  /// rejected as out-of-alphabet, the bug when reconstruction used the 20-state no-stop alphabet).
  #[test]
  fn test_aa_reconstruction_infer_handles_stop_codon_end_to_end() {
    let dir = tempdir().unwrap();
    let tree_path = dir.path().join("tree.nwk");
    let nuc_path = dir.path().join("aln.fasta");
    let aa_path = dir.path().join("aa.S.fasta");
    let gff_path = dir.path().join("annotation.gff3");

    fs::write(&tree_path, "(A:0.1,B:0.1)root;").unwrap();
    fs::write(&nuc_path, ">A\nACGT\n>B\nACGT\n").unwrap();
    fs::write(&aa_path, ">A\nMC*\n>B\nMA*\n").unwrap();
    fs::write(&gff_path, "##gff-version 3\nref\ttest\tCDS\t1\t9\t.\t+\t0\tName=S\n").unwrap();

    // The `{cds}` placeholder is a literal template token (substituted per CDS at runtime), built via
    // a named argument so it is not mistaken for a format-string argument.
    let translations = format!("{base}/aa.{token}.fasta", base = dir.path().display(), token = "{cds}");

    let args = TreetimeAncestralArgs {
      alignment: AlignmentArgs {
        alignment: vec![nuc_path],
      },
      tree: tree_path,
      method_anc: MethodAncestral::Marginal,
      dense: Some(false),
      translations: Some(translations),
      genes: vec!["S".to_owned()],
      annotation_gff: Some(gff_path),
      output: OutputArgs {
        outdir: dir.path().to_path_buf(),
        ..Default::default()
      },
      ..TreetimeAncestralArgs::default()
    };

    run_ancestral_reconstruction(&args, &NoopProgress).unwrap();

    let json = fs::read_to_string(dir.path().join("ancestral.augur-node-data.json")).unwrap();
    let data: AugurNodeDataJsonAncestral = json_read_str(&json).unwrap();

    // Per-CDS reference and annotation are present and keyed by the CDS name.
    let reference = data.metadata.reference.expect("reference present");
    assert!(reference.contains_key("S"));
    assert!(reference.contains_key("nuc"));
    let annotations = data.metadata.annotations.expect("annotations present");
    let cds = annotations.other.get("S").expect("CDS S annotation present");
    assert_eq!(Some(1), cds.start);
    assert_eq!(Some(9), cds.end);

    // The root carries the reconstructed AA sequence (root-only), of the peptide length, with the
    // stop codon preserved at the final position.
    let root = &data.nodes["root"];
    let root_aa = root.aa_sequences.as_ref().expect("root aa_sequences present");
    let root_s = &root_aa["S"];
    assert_eq!(3, root_s.chars().count());
    assert_eq!(Some('*'), root_s.chars().last());
    assert_eq!(Some('M'), root_s.chars().next());

    // Every node carries aa_muts for the CDS; non-root nodes never carry aa_sequences.
    for (name, node) in &data.nodes {
      assert!(node.aa_muts.as_ref().expect("aa_muts present").contains_key("S"));
      if name != "root" {
        assert!(node.aa_sequences.is_none());
      }
    }
  }
}
