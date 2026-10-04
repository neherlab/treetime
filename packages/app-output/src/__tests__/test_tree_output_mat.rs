#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, and scratch collections"
)]

#[cfg(test)]
mod tests {

  use crate::__tests__::test_tree_output::tests::helpers::{Mutations, all_mat_documents, branch_length, c};
  use crate::tree_output::{MatGapCounts, mat_mutation};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
  use treetime::seq::mutation::{Mutation, MutationTrack, Sub};
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read;
  use treetime_io::usher_mat::UsherMutation;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_tree_output_mat_leaves_out_amino_acid_mutations() -> Result<(), Report> {
    let expected = helpers::written_ancestral_mat(Mutations::NucleotideSubstitution)?;
    let actual = helpers::written_ancestral_mat(Mutations::NucleotideSubstitutionAndAminoAcid)?;
    assert_eq!(expected, actual);
    let positions: Vec<i32> = actual
      .node_mutations
      .iter()
      .flat_map(|mutations| &mutations.mutation)
      .map(|mutation| mutation.position)
      .collect();
    assert_eq!(vec![1], positions);
    Ok(())
  }

  #[test]
  fn test_tree_output_mat_rejects_amino_acid_mutation_as_internal_error() -> Result<(), Report> {
    let mutation = Mutation::substitution(MutationTrack::AminoAcid(o!("S")), Sub::new(c(b'A'), 0_usize, c(b'T'))?);
    assert_error!(
      mat_mutation(&mutation, Some("A"), &Alphabet::new(AlphabetName::Nuc)?, "A"),
      "Node 'A' has an amino-acid mutation, but UShER MAT stores nucleotide mutations only. This is an internal error. Please report it to developers."
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::deletion_is_missing_data(
    "ACGT", btreemap! { "B" => vec![helpers::del((1, 3), "CG")] },
    (btreemap! {}, MatGapCounts { deletions: 1, insertions: 0, substitutions: 0 }),
  )]
  #[case::reinsertion_restores_state_before_deletion(
    "ACGT", btreemap! { "B" => vec![helpers::del((1, 2), "C")], "C" => vec![helpers::ins((1, 2), "T")], "D" => vec![helpers::ins((1, 2), "C")] },
    (btreemap! { "C" => vec![helpers::usher(2, b'C', b'C', b'T')] }, MatGapCounts { deletions: 1, insertions: 0, substitutions: 0 }),
  )]
  #[case::reinsertion_sorted_with_substitutions(
    "ACGT", btreemap! { "B" => vec![helpers::del((1, 2), "C")], "C" => vec![helpers::sub(b'G', 2, b'A'), helpers::ins((1, 2), "T")] },
    (btreemap! { "C" => vec![helpers::usher(2, b'C', b'C', b'T'), helpers::usher(3, b'G', b'G', b'A')] }, MatGapCounts { deletions: 1, insertions: 0, substitutions: 0 }),
  )]
  #[case::insertion_in_root_gap_is_left_out(
    "A-GT", btreemap! { "B" => vec![helpers::ins((1, 2), "C")], "C" => vec![helpers::sub(b'C', 1, b'T')] },
    (btreemap! {}, MatGapCounts { deletions: 0, insertions: 1, substitutions: 1 }),
  )]
  #[case::reinsertion_in_root_gap_is_left_out(
    "A-GT", btreemap! { "B" => vec![helpers::ins((1, 2), "C")], "C" => vec![helpers::del((1, 2), "C")] },
    (btreemap! {}, MatGapCounts { deletions: 1, insertions: 1, substitutions: 0 }),
  )]
  #[case::deleted_unknown_base_hides_no_state(
    "ANGT", btreemap! { "B" => vec![helpers::del((1, 2), "N")], "C" => vec![helpers::ins((1, 2), "T")] },
    (btreemap! {}, MatGapCounts { deletions: 1, insertions: 0, substitutions: 0 }),
  )]
  #[trace]
  fn test_tree_output_mat_writes_gaps_as_missing_data(
    #[case] reference: &str,
    #[case] edge_mutations: BTreeMap<&str, Vec<Mutation>>,
    #[case] expected: (BTreeMap<&str, Vec<UsherMutation>>, MatGapCounts),
  ) -> Result<(), Report> {
    let actual = helpers::mat_mutations("((C:1,D:1)B:1,E:1)root;", reference, &edge_mutations)?;
    let expected = (
      expected.0.into_iter().map(|(name, mutations)| (name.to_owned(), mutations)).collect(),
      expected.1,
    );
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nothing_dropped(  MatGapCounts { deletions: 0, insertions: 0, substitutions: 0 }, None)]
  #[case::deletions(        MatGapCounts { deletions: 2, insertions: 0, substitutions: 0 }, Some("UShER MAT has no gap state: wrote 2 deletion(s) as missing data (N)"))]
  #[case::root_gap_columns( MatGapCounts { deletions: 0, insertions: 0, substitutions: 3 }, Some("UShER MAT has no gap state: left out 0 insertion(s) and 3 substitution(s) in alignment columns where the root sequence, which is the MAT reference, has a gap"))]
  #[case::all(              MatGapCounts { deletions: 1, insertions: 2, substitutions: 3 }, Some("UShER MAT has no gap state: wrote 1 deletion(s) as missing data (N); left out 2 insertion(s) and 3 substitution(s) in alignment columns where the root sequence, which is the MAT reference, has a gap"))]
  #[trace]
  fn test_tree_output_mat_gap_warning(#[case] gaps: MatGapCounts, #[case] expected: Option<&str>) {
    assert_eq!(expected.map(str::to_owned), gaps.warning());
  }

  #[test]
  fn test_tree_output_mat_uses_one_global_reference_for_recurrent_mutations() -> Result<(), Report> {
    let first = Mutation::substitution(MutationTrack::Nucleotide, Sub::new(c(b'A'), 0_usize, c(b'T'))?);
    let recurrent = Mutation::substitution(MutationTrack::Nucleotide, Sub::new(c(b'T'), 0_usize, c(b'C'))?);

    let first = mat_mutation(&first, Some("A"), &Alphabet::new(AlphabetName::Nuc)?, "inner")?;
    let recurrent = mat_mutation(&recurrent, Some("A"), &Alphabet::new(AlphabetName::Nuc)?, "leaf")?;
    assert_eq!((0, 0, vec![3]), (first.ref_nuc, first.par_nuc, first.mut_nuc));
    assert_eq!(
      (0, 3, vec![1]),
      (recurrent.ref_nuc, recurrent.par_nuc, recurrent.mut_nuc)
    );
    Ok(())
  }

  #[test]
  fn test_tree_output_mat_rejects_missing_reference() -> Result<(), Report> {
    let mutation = Mutation::substitution(MutationTrack::Nucleotide, Sub::new(c(b'A'), 0_usize, c(b'T'))?);
    let error = mat_mutation(&mutation, None, &Alphabet::new(AlphabetName::Nuc)?, "A")
      .expect_err("MAT must require a global reference");
    assert!(error.to_string().contains("requires a root nucleotide reference"));
    Ok(())
  }

  #[test]
  fn test_tree_output_mat_rejects_reference_lookup_out_of_range() -> Result<(), Report> {
    let mutation = Mutation::substitution(MutationTrack::Nucleotide, Sub::new(c(b'A'), 1_usize, c(b'T'))?);
    let error = mat_mutation(&mutation, Some("A"), &Alphabet::new(AlphabetName::Nuc)?, "A")
      .expect_err("MAT must check the reference length");
    assert!(error.to_string().contains("outside the root nucleotide reference"));
    Ok(())
  }

  #[test]
  fn test_tree_output_mat_rejects_coordinate_above_i32() -> Result<(), Report> {
    let position = usize::try_from(i32::MAX)?;
    let mutation = Mutation::substitution(MutationTrack::Nucleotide, Sub::new(c(b'A'), position, c(b'T'))?);
    let error = mat_mutation(&mutation, Some("A"), &Alphabet::new(AlphabetName::Nuc)?, "A")
      .expect_err("MAT must check its coordinate range");
    assert!(error.to_string().contains("exceeds the UShER MAT i32 coordinate range"));
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::root_reference(("N", b'A', b'T'), "root reference nucleotide 'N'")]
  #[case::parent(        ("A", b'N', b'T'), "parent nucleotide 'N'")]
  #[case::child_gap(     ("A", b'A', b'.'), "child nucleotide '.'")]
  #[trace]
  fn test_tree_output_mat_rejects_noncanonical_nucleotide(
    #[case] (reference, parent, child): (&str, u8, u8),
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let mutation = Mutation::substitution(
      MutationTrack::Nucleotide,
      Sub::new(c(parent), 0_usize, c(child))?,
    );
    let error = mat_mutation(&mutation, Some(reference), &Alphabet::new(AlphabetName::Nuc)?, "A").expect_err("MAT must reject a reference, parent, or child nucleotide it cannot encode");
    assert!(error.to_string().contains(expected));
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::canonical(b'T', vec![3])]
  #[case::ambiguous(b'K', vec![2, 3])]
  #[case::unknown(  b'N', vec![0, 1, 2, 3])]
  #[trace]
  fn test_tree_output_mat_encodes_child_state_as_its_canonical_nucleotides(
    #[case] child: u8,
    #[case] expected: Vec<i32>,
  ) -> Result<(), Report> {
    let mutation = Mutation::substitution(
      MutationTrack::Nucleotide,
      Sub::new(c(b'C'), 0_usize, c(child))?,
    );
    let actual = mat_mutation(&mutation, Some("C"), &Alphabet::new(AlphabetName::Nuc)?, "A")?;
    assert_eq!((1, 1, expected), (actual.ref_nuc, actual.par_nuc, actual.mut_nuc));
    Ok(())
  }

  #[test]
  fn test_tree_output_all_mat_models_preserve_embedded_newick_lengths() -> Result<(), Report> {
    let documents = all_mat_documents()?;
    assert_eq!(6, documents.len());
    assert!(documents.iter().all(|document| {
      document.node_mutations.len() == 4
        && document
          .node_mutations
          .iter()
          .all(|mutations| mutations.mutation.is_empty())
    }));

    for (command, document) in ["ancestral", "optimize", "prune", "clock", "mugration", "timetree"]
      .into_iter()
      .zip(documents)
    {
      let nwk_parsed = nwk_read(document.newick.as_bytes())?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;
      assert_eq!(
        None,
        branch_length(&graph, &names, &branch_lengths, "A")?,
        "{command}: {}",
        document.newick
      );
      assert_eq!(
        Some(0.0),
        branch_length(&graph, &names, &branch_lengths, "B")?,
        "{command}: {}",
        document.newick
      );
      assert_eq!(
        Some(0.5),
        branch_length(&graph, &names, &branch_lengths, "C")?,
        "{command}: {}",
        document.newick
      );
    }

    Ok(())
  }

  mod helpers {
    use crate::__tests__::test_tree_output::tests::helpers::{Mutations, ancestral_graph, ancestral_nodes};
    use crate::ancestral_tree_output::write_ancestral_tree_outputs;
    use crate::tree_output::{MatGapCounts, MatOutput, mat_from_graph};
    use eyre::{Report, WrapErr};
    use maplit::btreemap;
    use std::collections::BTreeMap;
    use tempfile::TempDir;
    use treetime::progress::NoopProgress;
    use treetime::seq::mutation::{AlignedMutation, Mutation, MutationEvent, MutationTrack, Sub};
    use treetime_io::graph::TreeWriteKind;
    use treetime_io::nwk::{CommentProviders, nwk_read};
    use treetime_io::usher_mat::{UsherMutation, UsherTree};
    use treetime_primitives::{AsciiChar, Seq};
    use treetime_utils::io::json::json_read_file;

    pub(super) fn written_ancestral_mat(mutations: Mutations) -> Result<UsherTree, Report> {
      let (graph, names, branch_lengths, maps, aa_node_data, aa_annotations) = ancestral_graph(mutations)?;
      let dir = TempDir::new().wrap_err("When creating a temporary directory")?;
      let path = dir.path().join("tree.mat.json");
      write_ancestral_tree_outputs(
        &graph,
        &ancestral_nodes(&names, &graph, &btreemap! {}),
        &branch_lengths,
        &maps,
        aa_node_data.as_ref(),
        &aa_annotations,
        &btreemap! { TreeWriteKind::MatJson => path.clone() },
        &CommentProviders::new(),
        &NoopProgress,
      )?;
      json_read_file(&path)
    }

    pub(super) fn mat_mutations(
      nwk: &str,
      reference: &str,
      edge_mutations: &BTreeMap<&str, Vec<Mutation>>,
    ) -> Result<(BTreeMap<String, Vec<UsherMutation>>, MatGapCounts), Report> {
      let parsed = nwk_read(nwk.as_bytes())?;
      let names = parsed.names();
      let MatOutput { tree, gaps } = mat_from_graph(
        &parsed.graph,
        &names,
        &parsed.branch_lengths,
        Some(reference),
        |node_key, _edge_key| {
          Ok(
            names[&node_key]
              .as_deref()
              .and_then(|name| edge_mutations.get(name))
              .cloned()
              .unwrap_or_default(),
          )
        },
      )?;
      let mutations = tree
        .condensed_nodes
        .into_iter()
        .zip(tree.node_mutations)
        .filter(|(_, mutations)| !mutations.mutation.is_empty())
        .map(|(node, mutations)| (node.node_name, mutations.mutation))
        .collect();
      Ok((mutations, gaps))
    }

    pub(super) fn sub(reff: u8, pos: usize, qry: u8) -> Mutation {
      Mutation::substitution(
        MutationTrack::Nucleotide,
        Sub::new(
          AsciiChar::from_byte_unchecked(reff),
          pos,
          AsciiChar::from_byte_unchecked(qry),
        )
        .unwrap(),
      )
    }

    pub(super) fn del(range: (usize, usize), sequence: &str) -> Mutation {
      aligned(MutationEvent::Deletion, range, sequence)
    }

    pub(super) fn ins(range: (usize, usize), sequence: &str) -> Mutation {
      aligned(MutationEvent::Insertion, range, sequence)
    }

    fn aligned(event: fn(AlignedMutation) -> MutationEvent, range: (usize, usize), sequence: &str) -> Mutation {
      Mutation {
        track: MutationTrack::Nucleotide,
        event: event(AlignedMutation {
          range,
          sequence: Seq::try_from_str(sequence).unwrap(),
        }),
      }
    }

    pub(super) fn usher(position: i32, reff: u8, par: u8, qry: u8) -> UsherMutation {
      UsherMutation {
        position,
        ref_nuc: nuc(reff),
        par_nuc: nuc(par),
        mut_nuc: vec![nuc(qry)],
        chromosome: String::new(),
      }
    }

    fn nuc(state: u8) -> i32 {
      b"ACGT".iter().position(|&nuc| nuc == state).unwrap() as i32
    }
  }
}
