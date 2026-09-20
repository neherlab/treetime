#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::fitch::compress_sequences;
  use crate::partition::fitch::partition::PartitionFitch;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use eyre::Report;
  use indoc::indoc;
  use itertools::Itertools;
  use maplit::btreemap;
  use treetime_graph::graph::Graph;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  type EdgeReport = (String, Vec<(String, usize)>, Vec<(usize, usize)>);

  fn compress(nwk: &str, fasta: &str) -> Result<Vec<EdgeReport>, Report> {
    let alphabet = Alphabet::default();
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(fasta, &alphabet)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read_str(nwk)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let mut partition = PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    compress_sequences(&graph, &mut partition, &node_seq_inputs(&graph, &names, aln))?;

    let name = |key| -> String { names.get(&key).cloned().flatten().unwrap_or_default() };

    Ok(
      graph
        .get_edges()
        .map(|edge| {
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
