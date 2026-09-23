#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {

  use crate::ancestral_tree_output::ancestral_to_mat;

  use crate::tree_output::mat_mutation;

  use eyre::Report;

  use pretty_assertions::assert_eq;
  use rstest::rstest;

  use treetime::seq::mutation::{Mutation, MutationTrack, Sub};
  use treetime_graph::graph::Graph;

  use treetime_io::nwk::nwk_read_str;

  use crate::__tests__::test_tree_output::tests::helpers;

  #[test]
  fn test_tree_output_mat_rejects_unsupported_events() -> Result<(), Report> {
    let (graph, names, branch_lengths, partition, aa_node_data, _aa_annotations) =
      helpers::ancestral_graph(helpers::Mutations::Indel)?;
    let error = ancestral_to_mat(
      &graph,
      &names,
      &branch_lengths,
      &helpers::ancestral_maps(&graph, partition.as_ref()),
      aa_node_data.as_ref(),
    )
    .expect_err("MAT must reject indels");
    assert!(error.to_string().contains("insertion or deletion"));

    let (graph, names, branch_lengths, partition, aa_node_data, _aa_annotations) =
      helpers::ancestral_graph(helpers::Mutations::AminoAcid)?;
    let error = ancestral_to_mat(
      &graph,
      &names,
      &branch_lengths,
      &helpers::ancestral_maps(&graph, partition.as_ref()),
      aa_node_data.as_ref(),
    )
    .expect_err("MAT must reject amino-acid mutations");
    assert!(error.to_string().contains("amino-acid mutation"));

    Ok(())
  }

  #[test]
  fn test_tree_output_mat_uses_one_global_reference_for_recurrent_mutations() -> Result<(), Report> {
    let first = Mutation::substitution(
      MutationTrack::Nucleotide,
      Sub::new(helpers::c(b'A'), 0_usize, helpers::c(b'T'))?,
    );
    let recurrent = Mutation::substitution(
      MutationTrack::Nucleotide,
      Sub::new(helpers::c(b'T'), 0_usize, helpers::c(b'C'))?,
    );

    let first = mat_mutation(&first, Some("A"), "inner")?;
    let recurrent = mat_mutation(&recurrent, Some("A"), "leaf")?;
    assert_eq!((0, 0, vec![3]), (first.ref_nuc, first.par_nuc, first.mut_nuc));
    assert_eq!(
      (0, 3, vec![1]),
      (recurrent.ref_nuc, recurrent.par_nuc, recurrent.mut_nuc)
    );
    Ok(())
  }

  #[test]
  fn test_tree_output_mat_rejects_missing_reference() -> Result<(), Report> {
    let mutation = Mutation::substitution(
      MutationTrack::Nucleotide,
      Sub::new(helpers::c(b'A'), 0_usize, helpers::c(b'T'))?,
    );
    let error = mat_mutation(&mutation, None, "A").expect_err("MAT must require a global reference");
    assert!(error.to_string().contains("requires a root nucleotide reference"));
    Ok(())
  }

  #[test]
  fn test_tree_output_mat_rejects_reference_lookup_out_of_range() -> Result<(), Report> {
    let mutation = Mutation::substitution(
      MutationTrack::Nucleotide,
      Sub::new(helpers::c(b'A'), 1_usize, helpers::c(b'T'))?,
    );
    let error = mat_mutation(&mutation, Some("A"), "A").expect_err("MAT must check the reference length");
    assert!(error.to_string().contains("outside the root nucleotide reference"));
    Ok(())
  }

  #[test]
  fn test_tree_output_mat_rejects_coordinate_above_i32() -> Result<(), Report> {
    let position = usize::try_from(i32::MAX)?;
    let mutation = Mutation::substitution(
      MutationTrack::Nucleotide,
      Sub::new(helpers::c(b'A'), position, helpers::c(b'T'))?,
    );
    let error = mat_mutation(&mutation, Some("A"), "A").expect_err("MAT must check its coordinate range");
    assert!(error.to_string().contains("exceeds the UShER MAT i32 coordinate range"));
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::root_reference(("N", b'A', b'T'), "root reference nucleotide 'N'")]
  #[case::parent(        ("A", b'N', b'T'), "parent nucleotide 'N'")]
  #[case::child(         ("A", b'A', b'N'), "child nucleotide 'N'")]
  #[trace]
  fn test_tree_output_mat_rejects_noncanonical_nucleotide(
    #[case] (reference, parent, child): (&str, u8, u8),
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let mutation = Mutation::substitution(
      MutationTrack::Nucleotide,
      Sub::new(helpers::c(parent), 0_usize, helpers::c(child))?,
    );
    let error = mat_mutation(&mutation, Some(reference), "A").expect_err("MAT must accept only A, C, G, or T");
    assert!(error.to_string().contains(expected));
    Ok(())
  }

  #[test]
  fn test_tree_output_mutation_free_mat_needs_no_reference() -> Result<(), Report> {
    let (graph, names, branch_lengths) = helpers::ancestral_graph_without_partition()?;
    let mat = ancestral_to_mat(
      &graph,
      &names,
      &branch_lengths,
      &helpers::ancestral_maps(&graph, None),
      None,
    )?;
    assert!(mat.node_mutations.iter().all(|mutations| mutations.mutation.is_empty()));
    Ok(())
  }

  #[test]
  fn test_tree_output_all_mat_models_preserve_embedded_newick_lengths() -> Result<(), Report> {
    let documents = helpers::all_mat_documents()?;
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
      let nwk_parsed = nwk_read_str(&document.newick)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;
      assert_eq!(
        None,
        helpers::branch_length(&graph, &names, &branch_lengths, "A")?,
        "{command}: {}",
        document.newick
      );
      assert_eq!(
        Some(0.0),
        helpers::branch_length(&graph, &names, &branch_lengths, "B")?,
        "{command}: {}",
        document.newick
      );
      assert_eq!(
        Some(0.5),
        helpers::branch_length(&graph, &names, &branch_lengths, "C")?,
        "{command}: {}",
        document.newick
      );
    }

    Ok(())
  }
}
