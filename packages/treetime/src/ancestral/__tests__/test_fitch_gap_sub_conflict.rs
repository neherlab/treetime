#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::fitch::compress_sequences;
  use crate::partition::fitch::partition::PartitionFitch;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use indoc::indoc;
  use itertools::Itertools;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  type EdgeReport = (String, Vec<(String, usize)>, Vec<(usize, usize)>);

  /// Run Fitch compression and return, per edge, the substitutions and deletion ranges.
  fn compress(nwk: &str, fasta: &str) -> Result<Vec<EdgeReport>, Report> {
    let alphabet = Alphabet::default();
    let aln = read_many_fasta_str(fasta, &alphabet)?;
    let graph: GraphAncestral = nwk_read_str(nwk)?;
    let partitions = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];
    compress_sequences(&graph, &partitions, &aln)?;

    let partition = partitions[0].read_arc();
    let name = |key| -> String {
      graph
        .get_node(key)
        .unwrap()
        .read_arc()
        .payload()
        .read_arc()
        .name
        .clone()
        .unwrap_or_default()
    };

    Ok(
      graph
        .get_edges()
        .iter()
        .map(|edge| {
          let edge = edge.read_arc();
          let data = &partition.edges[&edge.key()];
          (
            format!("{}->{}", name(edge.source()), name(edge.target())),
            data
              .fitch_subs()
              .iter()
              .map(|sub| (sub.to_string(), sub.pos()))
              .collect_vec(),
            data
              .indels
              .iter()
              .filter(|i| i.is_deletion())
              .map(|i| i.range)
              .collect_vec(),
          )
        })
        .collect_vec(),
    )
  }

  fn conflicts(edges: &[EdgeReport]) -> Vec<String> {
    edges
      .iter()
      .flat_map(|(name, subs, dels)| {
        subs
          .iter()
          .cartesian_product(dels)
          .filter(|((_, pos), (lo, hi))| pos >= lo && pos < hi)
          .map(move |((sub, pos), (lo, hi))| format!("{name}: {sub} (pos {pos}) inside deletion {lo}..{hi}"))
      })
      .collect_vec()
  }

  /// Field-data shape: a single determined residue stranded inside a run of unknowns, with sibling
  /// clades gapped across the whole region. This is the `PP_0012SSJ` / `PP_000ZH4A` pattern.
  ///
  /// `stray` has `G` at column 6 surrounded by `N`, so its `non_char` (union of unknown and gap
  /// ranges) has a one-column hole at 6. That hole survives every ancestral intersection. At `X`
  /// the column is only `variable_indel`; at `Y` a gapped sibling turns the whole block into a
  /// resolved gap -- `variable_indel` counts as gap-compatible -- while `non_char` keeps the hole.
  /// So column 6 is inside `Y`'s gap yet still carries `G`, and `G` is what gets substituted.
  #[test]
  fn stray_residue_in_unknown_run_does_not_create_sub_in_deletion() -> Result<(), Report> {
    let nwk = "(((stray:0.01,gapped:0.01)X:0.01,gapped2:0.01)Y:0.01,out:0.01)root:0.01;";
    let fasta = indoc! {r#"
      >stray
      AAAANNGNNAAA
      >gapped
      AAAA-----AAA
      >gapped2
      AAAA-----AAA
      >out
      AAAAAAAAAAAA
    "#};

    let edges = compress(nwk, fasta)?;
    for (name, subs, dels) in &edges {
      println!("{name}: subs={subs:?} deletions={dels:?}");
    }
    let conflicts = conflicts(&edges);
    assert!(conflicts.is_empty(), "conflicts:\n{}", conflicts.join("\n"));
    Ok(())
  }

  /// Gap-boundary variant: the descendants of `c2` disagree only about the *first* column of the
  /// deleted block (4), and agree that 5..8 is gapped. `X` resolves the whole block 4..8 to a
  /// gap, but `non_char` only covers 5..8, so column 4 keeps a canonical state and attracts a
  /// substitution -- landing exactly on the site where the deletion begins.
  #[test]
  fn sub_never_lands_at_deletion_start() -> Result<(), Report> {
    let nwk = "(((g1:0.01,g2:0.01,g3:0.01)c2:0.01,c1:0.01)X:0.01,out:0.01)root:0.01;";
    let fasta = indoc! {r#"
      >g1
      AAAA----AAAA
      >g2
      AAAAC---AAAA
      >g3
      AAAAT---AAAA
      >c1
      AAAA----AAAA
      >out
      AAAAAAAAAAAA
    "#};

    let edges = compress(nwk, fasta)?;
    for (name, subs, dels) in &edges {
      println!("{name}: subs={subs:?} deletions={dels:?}");
    }
    let conflicts = conflicts(&edges);
    assert!(conflicts.is_empty(), "conflicts:\n{}", conflicts.join("\n"));
    Ok(())
  }

  /// A Fitch substitution and a deletion on the same edge must not target the same site: the
  /// child either has a character there (substitution) or a gap (deletion), never both.
  ///
  /// Minimal construction. Node `X` has one child that is gapped over 4..8 (`c1`) and one child
  /// whose own descendants disagree about the gap (`c2`, giving it a `variable_indel` over
  /// 4..8). `resolve_indels_backward` treats "gapped + variable" as consensus gap, so `X` gets
  /// `gaps = [(4, 8)]` -- but `X.non_char` is the *intersection* of its children's `non_char`,
  /// which is empty over 4..8 because `c2` is not `non_char` there. So 4..8 is gapped but not
  /// masked, and the substitution machinery still sees canonical states there.
  #[test]
  fn sub_never_lands_inside_deletion_on_same_edge() -> Result<(), Report> {
    let nwk = "(((g1:0.01,g2:0.01,g3:0.01)c2:0.01,c1:0.01)X:0.01,out:0.01)root:0.01;";
    let fasta = indoc! {r#"
      >g1
      AAAA----AAAA
      >g2
      AAAACCCCAAAA
      >g3
      AAAATTTTAAAA
      >c1
      AAAA----AAAA
      >out
      AAAAAAAAAAAA
    "#};

    let edges = compress(nwk, fasta)?;
    for (name, subs, dels) in &edges {
      println!("{name}: subs={subs:?} deletions={dels:?}");
    }

    let conflicts = conflicts(&edges);
    assert!(conflicts.is_empty(), "conflicts:\n{}", conflicts.join("\n"));
    Ok(())
  }
}
