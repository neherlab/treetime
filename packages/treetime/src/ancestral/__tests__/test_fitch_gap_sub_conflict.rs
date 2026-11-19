#[cfg(test)]
mod tests {
  use eyre::Report;
  use helpers::{compress, conflicts};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::stray_residue_in_unknown_run(
    "(((stray:0.01,gapped:0.01)X:0.01,gapped2:0.01)Y:0.01,out:0.01)root:0.01;",
    indoc! {r#"
      >stray
      AAAANNGNNAAA
      >gapped
      AAAA-----AAA
      >gapped2
      AAAA-----AAA
      >out
      AAAAAAAAAAAA
    "#},
  )]
  #[case::residue_at_deletion_start(
    "(((g1:0.01,g2:0.01,g3:0.01)c2:0.01,c1:0.01)X:0.01,out:0.01)root:0.01;",
    indoc! {r#"
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
    "#},
  )]
  #[case::residues_across_deletion_on_same_edge(
    "(((g1:0.01,g2:0.01,g3:0.01)c2:0.01,c1:0.01)X:0.01,out:0.01)root:0.01;",
    indoc! {r#"
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
    "#},
  )]
  #[trace]
  fn test_fitch_places_no_substitution_inside_a_deletion_of_the_same_edge(
    #[case] nwk: &str,
    #[case] fasta: &str,
  ) -> Result<(), Report> {
    let edges = compress(nwk, fasta)?;
    assert_eq!(Vec::<String>::new(), conflicts(&edges));
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::partition::fitch::partition::PartitionFitch;
    use crate::partition::fitch::passes::compress_sequences;
    use crate::test_utils::{get_common_length, leaf_seq_inputs};
    use eyre::Report;
    use itertools::Itertools;
    use treetime_graph::keymap;
    use treetime_graph::scope::Inheritance;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::AlignmentRecord;

    pub(super) type EdgeReport = (String, Vec<(String, usize)>, Vec<(usize, usize)>);

    pub(super) fn compress(nwk: &str, fasta: &str) -> Result<Vec<EdgeReport>, Report> {
      let alphabet = Alphabet::default();
      let aln: Vec<AlignmentRecord> = fasta_read(fasta.as_bytes(), &alphabet)?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let nwk_parsed = nwk_read(nwk.as_bytes())?;
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let partition = PartitionFitch {
        alphabet,
        length: get_common_length(&aln)?,
        nodes: keymap! {},
        edges: keymap! {},
        inheritance: Inheritance::ALL,
      };
      let partition = compress_sequences(&graph, partition, leaf_seq_inputs(&graph, &names, aln))?;

      let name = |key| -> String { names[key].to_owned() };

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

    pub(super) fn conflicts(edges: &[EdgeReport]) -> Vec<String> {
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
  }
}
