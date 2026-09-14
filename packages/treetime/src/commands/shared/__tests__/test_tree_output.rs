#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::pipeline::AncestralPartition;
  use crate::commands::ancestral::aa_node_data::AaNodeData;
  use crate::commands::ancestral::result::AncestralNodeOut;
  use crate::commands::ancestral::tree_output::{
    ancestral_to_auspice, ancestral_to_mat, ancestral_to_phyloxml, write_ancestral_tree_outputs,
  };
  use crate::commands::clock::tree_output::{clock_to_auspice, clock_to_mat, clock_to_phyloxml};
  use crate::commands::mugration::tree_output::{mugration_to_auspice, mugration_to_mat, mugration_to_phyloxml};
  use crate::commands::optimize::tree_output::{optimize_to_auspice, optimize_to_mat, optimize_to_phyloxml};
  use crate::commands::prune::tree_output::{prune_to_auspice, prune_to_mat, prune_to_phyloxml};
  use crate::commands::shared::tree_output::{format_number, group_mutations, mat_mutation};
  use crate::commands::timetree::tree_output::{timetree_to_auspice, timetree_to_mat, timetree_to_phyloxml};
  use crate::partition::fitch::partition::PartitionFitch;
  use crate::partition::storage::sparse::{FitchNodeData, SparseEdgeObs};
  use crate::seq::indel::InDel;
  use crate::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub};
  use approx::assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::Value;
  use tempfile::TempDir;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::graph::TreeWriteKind;
  use treetime_io::nwk::{CommentProviders, NwkStyle, nwk_read_str};
  use treetime_primitives::{AsciiChar, LogLh, Seq};
  use treetime_utils::io::json::{JsonPretty, json_read_file, json_read_str, json_write_str};

  #[test]
  fn test_tree_output_ancestral_models_preserve_semantics() -> Result<(), Report> {
    let (graph, names, branch_lengths, partition, aa_node_data) =
      helpers::ancestral_graph(helpers::Mutations::NucleotideSubstitution)?;

    let nodes = helpers::ancestral_nodes(&names, &graph, &helpers::ancestral_confidences(&names, &graph));
    let auspice = ancestral_to_auspice(
      &graph,
      &nodes,
      &branch_lengths,
      &helpers::ancestral_maps(&graph, partition.as_ref()),
      aa_node_data.as_ref(),
      "2026-07-19",
    )?;
    let child = helpers::auspice_child(&auspice, "A");
    assert_eq!(Some("2026-07-19"), auspice.data.meta.updated.as_deref());
    assert_eq!(vec!["tree".to_owned()], auspice.data.meta.panels);
    assert_eq!(Some(0.5), child.node_attrs.div);
    assert_eq!(vec!["A1T".to_owned()], child.branch_attrs.mutations["nuc"]);

    let phyloxml = ancestral_to_phyloxml(
      &graph,
      &nodes,
      &branch_lengths,
      &helpers::ancestral_maps(&graph, partition.as_ref()),
      aa_node_data.as_ref(),
    )?;
    let child = helpers::phyloxml_child(&phyloxml, "A");
    assert_eq!(Some(0.5), child.branch_length_elem);
    assert_eq!(Some(0.9), child.confidence.first().map(|confidence| confidence.value));
    assert!(
      child
        .property
        .iter()
        .any(|property| { property.ref_ == "treetime:mutation" && property.value == "nuc:sub:A1T" })
    );
    assert_eq!(
      Some("TCG"),
      child
        .sequence
        .iter()
        .find(|sequence| sequence.name.as_deref() == Some("nuc"))
        .and_then(|sequence| sequence.mol_seq.as_ref())
        .map(|sequence| sequence.sequence.as_str())
    );

    let mat = ancestral_to_mat(
      &graph,
      &names,
      &branch_lengths,
      &helpers::ancestral_maps(&graph, partition.as_ref()),
      aa_node_data.as_ref(),
    )?;
    let mutation = mat
      .node_mutations
      .iter()
      .flat_map(|mutations| &mutations.mutation)
      .next()
      .expect("fixture must contain one MAT mutation");
    assert_eq!(1, mutation.position);
    assert_eq!(0, mutation.ref_nuc);
    assert_eq!(0, mutation.par_nuc);
    assert_eq!(vec![3], mutation.mut_nuc);

    Ok(())
  }

  #[test]
  fn test_tree_output_phyloxml_encodes_aa_track_and_grouped_indel() -> Result<(), Report> {
    let (graph, names, branch_lengths, partition, aa_node_data) =
      helpers::ancestral_graph(helpers::Mutations::IndelAndAminoAcid)?;

    let phyloxml = ancestral_to_phyloxml(
      &graph,
      &helpers::ancestral_nodes(&names, &graph, &btreemap! {}),
      &branch_lengths,
      &helpers::ancestral_maps(&graph, partition.as_ref()),
      aa_node_data.as_ref(),
    )?;
    let child = helpers::phyloxml_child(&phyloxml, "A");
    let properties = child
      .property
      .iter()
      .map(|property| property.value.as_str())
      .collect::<Vec<_>>();
    assert!(properties.contains(&"nuc:del:2-3:CG"));
    assert!(properties.contains(&"aa:S%2F1%3Aweird:sub:A2T"));

    // Nucleotide indels are dropped from the Auspice nuc mutation list, which mirrors the
    // substitution-only augur node-data muts. A branch whose only nucleotide change is a
    // deletion therefore has no `nuc` entry (phyloxml above still encodes it as `nuc:del:2-3:CG`).
    let (graph, names, branch_lengths, partition, aa_node_data) = helpers::ancestral_graph(helpers::Mutations::Indel)?;
    let auspice = ancestral_to_auspice(
      &graph,
      &helpers::ancestral_nodes(&names, &graph, &btreemap! {}),
      &branch_lengths,
      &helpers::ancestral_maps(&graph, partition.as_ref()),
      aa_node_data.as_ref(),
      "2026-07-19",
    )?;
    let child = helpers::auspice_child(&auspice, "A");
    assert!(!child.branch_attrs.mutations.contains_key("nuc"));

    let (graph, names, branch_lengths, partition, aa_node_data) =
      helpers::ancestral_graph(helpers::Mutations::AminoAcid)?;
    let auspice = ancestral_to_auspice(
      &graph,
      &helpers::ancestral_nodes(&names, &graph, &btreemap! {}),
      &branch_lengths,
      &helpers::ancestral_maps(&graph, partition.as_ref()),
      aa_node_data.as_ref(),
      "2026-07-19",
    )?;
    let child = helpers::auspice_child(&auspice, "A");
    assert_eq!(vec!["A2T".to_owned()], child.branch_attrs.mutations["S"]);
    let annotations = auspice
      .data
      .meta
      .genome_annotations
      .as_ref()
      .expect("AA fixture must have genome annotations");
    assert!(annotations.nuc.is_some());
    assert!(annotations.cdses.contains_key("S"));

    Ok(())
  }

  #[test]
  fn test_tree_output_mat_rejects_unsupported_events() -> Result<(), Report> {
    let (graph, names, branch_lengths, partition, aa_node_data) = helpers::ancestral_graph(helpers::Mutations::Indel)?;
    let error = ancestral_to_mat(
      &graph,
      &names,
      &branch_lengths,
      &helpers::ancestral_maps(&graph, partition.as_ref()),
      aa_node_data.as_ref(),
    )
    .expect_err("MAT must reject indels");
    assert!(error.to_string().contains("insertion or deletion"));

    let (graph, names, branch_lengths, partition, aa_node_data) =
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
  fn test_tree_output_conversion_failure_does_not_create_target_and_keeps_prior_file() -> Result<(), Report> {
    let (graph, names, branch_lengths, partition, aa_node_data) = helpers::ancestral_graph(helpers::Mutations::Indel)?;
    let dir = TempDir::new()?;
    let nwk_path = dir.path().join("tree.nwk");
    let mat_path = dir.path().join("tree.mat.json");
    let outputs = btreemap! {
      TreeWriteKind::nwk(NwkStyle::Plain) => nwk_path.clone(),
      TreeWriteKind::MatJson => mat_path.clone(),
    };

    let error = write_ancestral_tree_outputs(
      &graph,
      &helpers::ancestral_nodes(&names, &graph, &btreemap! {}),
      &branch_lengths,
      &helpers::ancestral_maps(&graph, partition.as_ref()),
      aa_node_data.as_ref(),
      &outputs,
      &CommentProviders::new(),
    )
    .expect_err("MAT conversion must fail");
    assert!(error.to_string().contains("insertion or deletion"));
    assert!(nwk_path.is_file());
    assert!(!mat_path.exists());

    Ok(())
  }

  #[test]
  fn test_tree_output_graph_json_dumps_topology() -> Result<(), Report> {
    let (graph, names, branch_lengths, partition, aa_node_data) = helpers::ancestral_graph(helpers::Mutations::None)?;
    let dir = TempDir::new()?;
    let path = dir.path().join("graph.json");
    let outputs = btreemap! { TreeWriteKind::GraphJson => path.clone() };

    write_ancestral_tree_outputs(
      &graph,
      &helpers::ancestral_nodes(&names, &graph, &btreemap! {}),
      &branch_lengths,
      &helpers::ancestral_maps(&graph, partition.as_ref()),
      aa_node_data.as_ref(),
      &outputs,
      &CommentProviders::new(),
    )?;
    let actual: Value = json_read_file(&path)?;
    // GraphJson serializes pure topology: the node and edge sets, with no command data slot.
    let nodes = actual["nodes"]
      .as_array()
      .expect("graph.json must carry the node topology");
    assert!(!nodes.is_empty());
    let edges = actual["edges"]
      .as_array()
      .expect("graph.json must carry the edge topology");
    assert!(!edges.is_empty());

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

  // Oracle: Nextstrain Augur's dataset v2 schema at
  // d8e38736037ba9474a809f9a5a63bc2b279d2407.
  #[test]
  fn test_tree_output_all_auspice_models_match_augur_v2_schema() -> Result<(), Report> {
    let documents = helpers::all_auspice_documents()?;
    let validator = helpers::auspice_validator()?;

    for (command, document) in ["ancestral", "optimize", "prune", "clock", "mugration", "timetree"]
      .into_iter()
      .zip(&documents)
    {
      let errors = validator
        .iter_errors(document)
        .map(|error| error.to_string())
        .collect::<Vec<_>>()
        .join("\n");
      assert!(errors.is_empty(), "{command} Auspice schema errors:\n{errors}");
    }

    let mut malformed = documents[0].clone();
    malformed["meta"]
      .as_object_mut()
      .expect("Auspice meta must be an object")
      .remove("updated");
    assert!(!validator.is_valid(&malformed));

    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_rejects_node_without_divergence_or_date() -> Result<(), Report> {
    let error = helpers::optimize_auspice_without_required_node_data()
      .expect_err("Auspice must reject a node with neither divergence nor numerical date");
    assert!(error.to_string().contains("requires divergence or numerical date"));
    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_rejects_invalid_amino_acid_track_name() -> Result<(), Report> {
    let (graph, names, branch_lengths, partition, aa_node_data) =
      helpers::ancestral_graph(helpers::Mutations::IndelAndAminoAcid)?;
    let error = ancestral_to_auspice(
      &graph,
      &helpers::ancestral_nodes(&names, &graph, &btreemap! {}),
      &branch_lengths,
      &helpers::ancestral_maps(&graph, partition.as_ref()),
      aa_node_data.as_ref(),
      "2026-07-19",
    )
    .expect_err("Auspice must reject an amino-acid track outside its schema grammar");
    assert!(error.to_string().contains("cannot represent amino-acid mutation track"));
    Ok(())
  }

  #[test]
  fn test_tree_output_all_phyloxml_models_have_one_rooted_phylogeny() -> Result<(), Report> {
    let documents = helpers::all_phyloxml_documents()?;
    assert_eq!(6, documents.len());
    assert!(documents.iter().all(|document| {
      document.phylogeny.len() == 1 && document.phylogeny[0].rooted && document.phylogeny[0].clade.is_some()
    }));
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

  // Oracle: augur export_v2.format_number keeps `precision` significant figures in the
  // fractional part while preserving integer digits.
  #[test]
  fn test_tree_output_format_number_fractional_precision() {
    assert_ulps_eq!(0.123457, format_number(0.12345678, 6), max_ulps = 0);
    assert_ulps_eq!(123.456789, format_number(123.456789, 6), max_ulps = 0);
    assert_ulps_eq!(0.0, format_number(0.0, 6), max_ulps = 0);
    assert_ulps_eq!(2020.123, format_number(2020.1234567, 3), max_ulps = 0);
  }

  #[test]
  fn test_tree_output_group_mutations_drops_nucleotide_indels_keeps_amino_acid_indels() -> Result<(), Report> {
    // Auspice v2 mutation lists mirror the augur node-data `muts`: substitution-only for
    // the nucleotide track (augur export copies node-data nuc muts verbatim), indels retained
    // for amino-acid tracks (aa node-data emits them). Deletion of range (1, 3) over "CG"
    // expands to per-position tokens "C2-", "G3-".
    let mutations = vec![
      Mutation::substitution(
        MutationTrack::Nucleotide,
        Sub::new(helpers::c(b'A'), 0_usize, helpers::c(b'T'))?,
      ),
      Mutation::indel(
        MutationTrack::Nucleotide,
        &InDel::del((1, 3), Seq::try_from_str("CG")?)?,
      )?,
      Mutation::substitution(
        MutationTrack::AminoAcid("GENE".to_owned()),
        Sub::new(helpers::c(b'K'), 4_usize, helpers::c(b'R'))?,
      ),
      Mutation::indel(
        MutationTrack::AminoAcid("GENE".to_owned()),
        &InDel::del((1, 3), Seq::try_from_str("CG")?)?,
      )?,
    ];

    let grouped = group_mutations(mutations)?;

    let expected = btreemap! {
      "GENE".to_owned() => vec!["K5R".to_owned(), "C2-".to_owned(), "G3-".to_owned()],
      "nuc".to_owned() => vec!["A1T".to_owned()],
    };
    assert_eq!(expected, grouped);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::commands::ancestral::result::AncestralOutputMaps;
    use crate::commands::ancestral::run::gather_ancestral_output_maps;
    use crate::commands::clock::run::ClockNodeOut;
    use crate::commands::optimize::result::{OptimizeNodeOut, OptimizeOutputMaps};
    use crate::commands::optimize::run::gather_optimize_output_maps;
    use crate::commands::prune::result::{PruneNodeOut, PruneOutputMaps};
    use crate::commands::prune::run::gather_prune_output_maps;
    use crate::commands::timetree::result::{TimetreeEdgeOut, TimetreeNodeOut, TimetreeOutputMaps};
    use crate::commands::timetree::run::gather_timetree_output_maps;
    use crate::gtr::gtr::{GTR, GTRParams};
    use crate::mugration::result::gather_mugration_output_maps;
    use crate::mugration::result::{MugrationNodeOut, MugrationOutputMaps, MugrationResult};
    use crate::partition::marginal::discrete::partition::PartitionMarginalDiscrete;
    use crate::partition::storage::dense::{DenseNodeState, DenseSeqDistribution, DenseSeqInfo};
    use crate::partition::storage::discrete::DiscreteStates;
    use jsonschema::{Retrieve, Uri, Validator};
    use ndarray::array;
    use serde::Serialize;
    use std::collections::BTreeMap;
    use std::error::Error as StdError;
    use std::io;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_io::auspice_types::{AuspiceTree, AuspiceTreeNode};
    use treetime_io::phyloxml::{Phyloxml, PhyloxmlClade};
    use treetime_io::usher_mat::UsherTree;
    use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

    const AUSPICE_SCHEMA: &str = include_str!("schemas/auspice/schema-export-v2.json");
    const AUSPICE_CONFIG_SCHEMA: &str = include_str!("schemas/auspice/schema-auspice-config-v2.json");
    const ANNOTATIONS_SCHEMA: &str = include_str!("schemas/auspice/schema-annotations.json");
    const ROOT_SEQUENCE_SCHEMA: &str = include_str!("schemas/auspice/schema-export-root-sequence.json");
    const MODEL_TREE: &str = "(A:0.1,B:0,C:0.5)root;";

    #[derive(Clone, Copy)]
    pub enum Mutations {
      None,
      NucleotideSubstitution,
      Indel,
      AminoAcid,
      IndelAndAminoAcid,
    }

    pub fn ancestral_maps(graph: &Graph, partition: Option<&AncestralPartition>) -> AncestralOutputMaps {
      gather_ancestral_output_maps(graph, partition).unwrap()
    }

    pub fn optimize_maps(graph: &Graph) -> OptimizeOutputMaps {
      gather_optimize_output_maps(graph, &[], &[]).unwrap()
    }

    pub fn prune_maps(graph: &Graph) -> PruneOutputMaps {
      gather_prune_output_maps(graph, &[]).unwrap()
    }

    pub fn timetree_maps(graph: &Graph) -> TimetreeOutputMaps {
      gather_timetree_output_maps(graph, &[]).unwrap()
    }

    pub fn mugration_maps(
      graph: &Graph,
      partition: &PartitionMarginalDiscrete,
      node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    ) -> MugrationOutputMaps {
      gather_mugration_output_maps(graph, partition, node_states)
    }

    type AncestralGraphSetup = (
      Graph,
      BTreeMap<GraphNodeKey, Option<String>>,
      BTreeMap<GraphEdgeKey, Option<f64>>,
      Option<AncestralPartition>,
      Option<AaNodeData>,
    );

    pub fn ancestral_graph(mutations: Mutations) -> Result<AncestralGraphSetup, Report> {
      let nwk_parsed = nwk_read_str("(A:0.5,B:0)root;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;
      let root_key = node_key(&graph, &names, "root");
      let a_key = node_key(&graph, &names, "A");
      let b_key = node_key(&graph, &names, "B");
      let a_edge = graph.node_parent(a_key)?.unwrap().1;
      let b_edge = graph.node_parent(b_key)?.unwrap().1;
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let root_sequence = Seq::try_from_str("ACG")?;
      let a_sequence = Seq::try_from_str("TCG")?;
      let b_sequence = root_sequence.clone();

      let include_substitution = matches!(mutations, Mutations::NucleotideSubstitution);
      let include_indel = matches!(mutations, Mutations::Indel | Mutations::IndelAndAminoAcid);
      let include_aa = matches!(mutations, Mutations::AminoAcid | Mutations::IndelAndAminoAcid);
      let mut a_edge_data = if include_substitution {
        SparseEdgeObs::with_fitch_subs(vec![Sub::new(c(b'A'), 0_usize, c(b'T'))?])
      } else {
        SparseEdgeObs::default()
      };
      if include_indel {
        a_edge_data.indels = vec![InDel::del((1, 3), Seq::try_from_str("CG")?)?];
      }
      let partition = PartitionFitch {
        index: 0,
        alphabet: alphabet.clone(),
        length: 3,
        nodes: btreemap! {
          root_key => FitchNodeData::new(&root_sequence, &alphabet)?,
          a_key => FitchNodeData::new(&a_sequence, &alphabet)?,
          b_key => FitchNodeData::new(&b_sequence, &alphabet)?,
        },
        edges: btreemap! {
          a_edge => a_edge_data,
          b_edge => SparseEdgeObs::default(),
        },
      };

      let aa_node_data = include_aa.then(|| {
        let mut aa = AaNodeData::default();
        let track = if matches!(mutations, Mutations::IndelAndAminoAcid) {
          "S/1:weird"
        } else {
          "S"
        };
        aa.annotations.insert(
          "S".to_owned(),
          AugurNodeDataJsonAnnotationEntry {
            start: Some(1),
            end: Some(3),
            strand: Some("+".to_owned()),
            entry_type: Some("CDS".to_owned()),
            ..AugurNodeDataJsonAnnotationEntry::default()
          },
        );
        aa.root_aa_sequences.insert(track.to_owned(), "AA".to_owned());
        aa.node_aa_mutations.insert(
          a_key,
          btreemap! {
            track.to_owned() => vec![MutationEvent::Substitution(
              Sub::new(c(b'A'), 1_usize, c(b'T')).unwrap(),
            )],
          },
        );
        aa
      });
      let ancestral_partition = AncestralPartition::Fitch(partition);
      Ok((graph, names, branch_lengths, Some(ancestral_partition), aa_node_data))
    }

    pub fn ancestral_nodes(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    ) -> BTreeMap<GraphNodeKey, AncestralNodeOut> {
      graph
        .get_nodes()
        .iter()
        .map(|node| {
          let node = node.read_arc();
          let key = node.key();
          (
            key,
            AncestralNodeOut {
              name: names.get(&node.key()).cloned().flatten(),
              confidence: confidences.get(&key).copied().flatten(),
            },
          )
        })
        .collect()
    }

    /// Input-tree branch support for the ancestral fixture: node `A` carries 0.9, every other node
    /// none. Mirrors a Newick parse that annotated only `A`, so the output writers surface 0.9 on
    /// `A` and nothing elsewhere.
    pub fn ancestral_confidences(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
    ) -> BTreeMap<GraphNodeKey, Option<f64>> {
      graph
        .get_nodes()
        .iter()
        .filter_map(|node| {
          let node = node.read_arc();
          (names.get(&node.key()).and_then(|x| x.as_deref()) == Some("A")).then(|| (node.key(), Some(0.9)))
        })
        .collect()
    }

    pub fn ancestral_graph_without_partition() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str(MODEL_TREE)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      Ok((graph, names, branch_lengths))
    }

    pub fn all_auspice_documents() -> Result<Vec<Value>, Report> {
      let (ancestral_graph, ancestral_names, ancestral_bl, ancestral_partition, ancestral_aa) =
        ancestral_graph(Mutations::NucleotideSubstitution)?;
      let ancestral = ancestral_to_auspice(
        &ancestral_graph,
        &ancestral_nodes(&ancestral_names, &ancestral_graph, &btreemap! {}),
        &ancestral_bl,
        &ancestral_maps(&ancestral_graph, ancestral_partition.as_ref()),
        ancestral_aa.as_ref(),
        "2026-07-19",
      )?;
      let (optimize_graph, optimize_names, optimize_bl) = optimize_graph()?;
      let optimize = optimize_to_auspice(
        &optimize_graph,
        &optimize_nodes(&optimize_names, &optimize_graph, &btreemap! {}),
        &optimize_bl,
        &optimize_maps(&optimize_graph),
        "2026-07-19",
      )?;
      let (prune_graph, prune_names, prune_bl) = prune_graph()?;
      let prune = prune_to_auspice(
        &prune_graph,
        &prune_nodes(&prune_names, &prune_graph, &btreemap! {}),
        &prune_bl,
        &prune_maps(&prune_graph),
        "2026-07-19",
      )?;
      let (clock_graph, clock_names, _clock_bl) = clock_graph()?;
      let clock = clock_to_auspice(&clock_graph, &clock_nodes(&clock_names, &clock_graph), "2026-07-19")?;
      let (mugration_graph, mugration_names, mugration_bl, mugration_partition, mugration_node_states) =
        mugration_graph()?;
      let mugration = mugration_to_auspice(
        &mugration_graph,
        &mugration_nodes(&mugration_names, &mugration_graph, &btreemap! {}),
        &mugration_bl,
        &mugration_maps(&mugration_graph, &mugration_partition, &mugration_node_states),
        "country",
        "2026-07-19",
      )?;
      let (timetree_graph, timetree_names, _timetree_bl) = timetree_graph()?;
      let timetree = timetree_to_auspice(
        &timetree_graph,
        &timetree_nodes(&timetree_names, &timetree_graph, &btreemap! {}),
        &timetree_maps(&timetree_graph),
        None,
        None,
        "2026-07-19",
      )?;

      [ancestral, optimize, prune, clock, mugration, timetree]
        .iter()
        .map(json_value)
        .collect()
    }

    pub fn all_phyloxml_documents() -> Result<Vec<Phyloxml>, Report> {
      let (ancestral_graph, ancestral_names, ancestral_bl) = ancestral_graph_without_partition()?;
      Ok(vec![
        ancestral_to_phyloxml(
          &ancestral_graph,
          &ancestral_nodes(&ancestral_names, &ancestral_graph, &btreemap! {}),
          &ancestral_bl,
          &ancestral_maps(&ancestral_graph, None),
          None,
        )?,
        {
          let (optimize_graph, optimize_names, optimize_bl) = optimize_graph()?;
          optimize_to_phyloxml(
            &optimize_graph,
            &optimize_nodes(&optimize_names, &optimize_graph, &btreemap! {}),
            &optimize_bl,
            &optimize_maps(&optimize_graph),
          )?
        },
        {
          let (prune_graph, prune_names, prune_bl) = prune_graph()?;
          prune_to_phyloxml(
            &prune_graph,
            &prune_nodes(&prune_names, &prune_graph, &btreemap! {}),
            &prune_bl,
            &prune_maps(&prune_graph),
          )?
        },
        {
          let (clock_graph, clock_names, clock_bl) = clock_graph()?;
          clock_to_phyloxml(&clock_graph, &clock_nodes(&clock_names, &clock_graph), &clock_bl)?
        },
        {
          let (mugration_graph, mugration_names, mugration_bl, mugration_partition, mugration_node_states) =
            mugration_graph()?;
          mugration_to_phyloxml(
            &mugration_graph,
            &mugration_nodes(&mugration_names, &mugration_graph, &btreemap! {}),
            &mugration_bl,
            &mugration_maps(&mugration_graph, &mugration_partition, &mugration_node_states),
            "country",
          )?
        },
        {
          let (timetree_graph, timetree_names, timetree_bl) = timetree_graph()?;
          timetree_to_phyloxml(
            &timetree_graph,
            &timetree_nodes(&timetree_names, &timetree_graph, &btreemap! {}),
            &timetree_edges(&timetree_graph, &timetree_bl),
            &timetree_maps(&timetree_graph),
            None,
            None,
            None,
          )?
        },
      ])
    }

    pub fn all_mat_documents() -> Result<Vec<UsherTree>, Report> {
      let (ancestral, ancestral_names, mut ancestral_bl) = ancestral_graph_without_partition()?;
      set_mat_branch_lengths(&ancestral, &ancestral_names, &mut ancestral_bl)?;
      let (optimize, optimize_names, mut optimize_bl) = optimize_graph()?;
      set_mat_branch_lengths(&optimize, &optimize_names, &mut optimize_bl)?;
      let (prune, prune_names, mut prune_bl) = prune_graph()?;
      set_mat_branch_lengths(&prune, &prune_names, &mut prune_bl)?;
      let (clock, clock_names, mut clock_bl) = clock_graph()?;
      set_mat_branch_lengths(&clock, &clock_names, &mut clock_bl)?;
      let (mugration, mugration_names, mut mugration_bl, _mugration_partition, _mugration_node_states) =
        mugration_graph()?;
      set_mat_branch_lengths(&mugration, &mugration_names, &mut mugration_bl)?;
      let (timetree, timetree_names, _timetree_bl) = timetree_graph()?;
      let timetree_weights = timetree_mat_nwk_weights(&timetree, &timetree_names)?;

      Ok(vec![
        ancestral_to_mat(
          &ancestral,
          &ancestral_names,
          &ancestral_bl,
          &ancestral_maps(&ancestral, None),
          None,
        )?,
        optimize_to_mat(&optimize, &optimize_names, &optimize_bl, &optimize_maps(&optimize))?,
        prune_to_mat(&prune, &prune_names, &prune_bl, &prune_maps(&prune))?,
        clock_to_mat(&clock, &clock_names, &clock_bl)?,
        mugration_to_mat(&mugration, &mugration_names, &mugration_bl)?,
        timetree_to_mat(&timetree, &timetree_names, &timetree_weights, &timetree_maps(&timetree))?,
      ])
    }

    pub fn auspice_validator() -> Result<Validator, Report> {
      let schema = json_read_str(AUSPICE_SCHEMA)?;
      Ok(
        jsonschema::draft6::options()
          .with_retriever(AuspiceSchemaRetriever::new()?)
          .build(&schema)?,
      )
    }

    pub fn optimize_auspice_without_required_node_data() -> Result<AuspiceTree, Report> {
      let (graph, names, mut branch_lengths) = optimize_graph()?;
      set_branch_length(&graph, &names, &mut branch_lengths, "A", None)?;
      optimize_to_auspice(
        &graph,
        &optimize_nodes(&names, &graph, &btreemap! {}),
        &branch_lengths,
        &optimize_maps(&graph),
        "2026-07-19",
      )
    }

    pub fn auspice_child<'a>(tree: &'a AuspiceTree, name: &str) -> &'a AuspiceTreeNode {
      tree
        .tree
        .children
        .iter()
        .find(|child| child.name == name)
        .expect("fixture child must exist")
    }

    pub fn phyloxml_child<'a>(tree: &'a Phyloxml, name: &str) -> &'a PhyloxmlClade {
      tree.phylogeny[0]
        .clade
        .as_ref()
        .expect("PhyloXML fixture must have a root")
        .clade
        .iter()
        .find(|child| child.name.as_deref() == Some(name))
        .expect("fixture child must exist")
    }

    pub fn branch_length(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
      target_name: &str,
    ) -> Result<Option<f64>, Report> {
      let node_key = node_key(graph, names, target_name);
      let edge_key = graph
        .node_parent(node_key)?
        .expect("fixture target must not be the root")
        .1;
      Ok(branch_lengths.get(&edge_key).copied().flatten())
    }

    fn node_key(graph: &Graph, names: &BTreeMap<GraphNodeKey, Option<String>>, name: &str) -> GraphNodeKey {
      graph
        .get_nodes()
        .into_iter()
        .find_map(|node| {
          let node = node.read_arc();
          (names.get(&node.key()).and_then(|x| x.as_deref()) == Some(name)).then(|| node.key())
        })
        .expect("fixture node must exist")
    }

    pub fn c(value: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(value)
    }

    pub fn optimize_nodes(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    ) -> BTreeMap<GraphNodeKey, OptimizeNodeOut> {
      graph
        .get_nodes()
        .iter()
        .map(|node| {
          let node = node.read_arc();
          let key = node.key();
          (
            key,
            OptimizeNodeOut {
              name: names.get(&node.key()).cloned().flatten(),
              confidence: confidences.get(&key).copied().flatten(),
            },
          )
        })
        .collect()
    }

    fn optimize_graph() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str(MODEL_TREE)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      Ok((graph, names, branch_lengths))
    }

    pub fn prune_nodes(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    ) -> BTreeMap<GraphNodeKey, PruneNodeOut> {
      graph
        .get_nodes()
        .iter()
        .map(|node| {
          let node = node.read_arc();
          let key = node.key();
          (
            key,
            PruneNodeOut {
              name: names.get(&node.key()).cloned().flatten(),
              confidence: confidences.get(&key).copied().flatten(),
            },
          )
        })
        .collect()
    }

    fn prune_graph() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str(MODEL_TREE)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      Ok((graph, names, branch_lengths))
    }

    fn clock_graph() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str(MODEL_TREE)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      Ok((graph, names, branch_lengths))
    }

    pub fn clock_nodes(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
    ) -> BTreeMap<GraphNodeKey, ClockNodeOut> {
      graph
        .get_nodes()
        .iter()
        .enumerate()
        .map(|(index, node)| {
          let node = node.read_arc();
          let key = node.key();
          (
            key,
            ClockNodeOut {
              name: names.get(&node.key()).cloned().flatten(),
              div: index as f64 / 2.0,
              time: Some(2020.0 + index as f64),
              is_outlier: false,
              bad_branch: false,
            },
          )
        })
        .collect()
    }

    pub fn mugration_nodes(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    ) -> BTreeMap<GraphNodeKey, MugrationNodeOut> {
      graph
        .get_nodes()
        .iter()
        .map(|node| {
          let node = node.read_arc();
          let key = node.key();
          (
            key,
            MugrationNodeOut {
              name: names.get(&node.key()).cloned().flatten(),
              confidence: confidences.get(&key).copied().flatten(),
            },
          )
        })
        .collect()
    }

    fn mugration_graph() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
        PartitionMarginalDiscrete,
        BTreeMap<GraphNodeKey, DenseNodeState>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str(MODEL_TREE)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;
      let states = DiscreteStates::from_values(["CH", "US"].into_iter(), "?");
      let gtr = GTR::new(GTRParams {
        n_states: 2,
        mu: 1.0,
        W: None,
        pi: array![0.5, 0.5],
      })?;
      let partition = PartitionMarginalDiscrete::new(gtr, states, 1e-8, false);
      let node_states: BTreeMap<GraphNodeKey, DenseNodeState> = graph
        .get_nodes()
        .into_iter()
        .enumerate()
        .map(|(index, node)| {
          let key = node.read_arc().key();
          let profile = if index % 2 == 0 {
            array![[1.0, 0.0]]
          } else {
            array![[0.0, 1.0]]
          };
          (
            key,
            DenseNodeState {
              seq: DenseSeqInfo::default(),
              profile: DenseSeqDistribution::new(profile, LogLh::ZERO),
            },
          )
        })
        .collect();
      let result = MugrationResult::new(
        graph,
        &btreemap! {},
        &names,
        &branch_lengths,
        &partition,
        &node_states,
        "country",
      );
      Ok((result.graph, names, branch_lengths, partition, node_states))
    }

    fn timetree_graph() -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str(MODEL_TREE)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      Ok((graph, names, branch_lengths))
    }

    pub fn timetree_nodes(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    ) -> BTreeMap<GraphNodeKey, TimetreeNodeOut> {
      graph
        .get_nodes()
        .iter()
        .enumerate()
        .map(|(index, node)| {
          let node = node.read_arc();
          let key = node.key();
          (
            key,
            TimetreeNodeOut {
              name: names.get(&node.key()).cloned().flatten(),
              desc: None,
              confidence: confidences.get(&key).copied().flatten(),
              time: Some(2020.0 + index as f64),
              div: index as f64 / 2.0,
              is_outlier: false,
              bad_branch: false,
              // Rate-susceptibility dates are produced only by the confidence pass and threaded as a
              // value map; this fixture graph runs no such pass, so production surfaces None here too.
              rate_susceptibility_dates: None,
            },
          )
        })
        .collect()
    }

    pub fn timetree_edges(
      graph: &Graph,
      branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    ) -> BTreeMap<GraphEdgeKey, TimetreeEdgeOut> {
      graph
        .get_edges()
        .iter()
        .map(|edge| {
          let key = edge.read_arc().key();
          (
            key,
            TimetreeEdgeOut {
              branch_length: branch_lengths.get(&key).copied().flatten(),
              time_length: None,
              clock_branch_length: None,
              // Strict-clock test graph: the relaxed-clock multiplier is its default 1.0, matching
              // what production reads from the threaded edge state.
              gamma: 1.0,
            },
          )
        })
        .collect()
    }

    fn set_mat_branch_lengths(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
    ) -> Result<(), Report> {
      set_branch_length(graph, names, branch_lengths, "A", None)?;
      set_branch_length(graph, names, branch_lengths, "B", Some(0.0))?;
      set_branch_length(graph, names, branch_lengths, "C", Some(0.5))
    }

    fn set_branch_length(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
      name: &str,
      length: Option<f64>,
    ) -> Result<(), Report> {
      let key = graph
        .get_nodes()
        .into_iter()
        .find_map(|node| {
          let node = node.read_arc();
          (names.get(&node.key()).and_then(|x| x.as_deref()) == Some(name)).then(|| node.key())
        })
        .expect("fixture node must exist");
      let edge_key = graph.node_parent(key)?.expect("fixture node must have a parent").1;
      branch_lengths.insert(edge_key, length);
      Ok(())
    }

    fn timetree_mat_nwk_weights(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
    ) -> Result<BTreeMap<GraphEdgeKey, Option<f64>>, Report> {
      let mut weights: BTreeMap<GraphEdgeKey, Option<f64>> = graph
        .get_edges()
        .iter()
        .map(|edge| (edge.read_arc().key(), None))
        .collect();
      for (name, length) in [("A", None), ("B", Some(0.0)), ("C", Some(0.5))] {
        let key = graph
          .get_nodes()
          .into_iter()
          .find_map(|node| {
            let node = node.read_arc();
            (names.get(&node.key()).and_then(|x| x.as_deref()) == Some(name)).then(|| node.key())
          })
          .expect("fixture node must exist");
        let edge_key = graph.node_parent(key)?.expect("fixture node must have a parent").1;
        weights.insert(edge_key, length);
      }
      Ok(weights)
    }

    fn json_value(value: &impl Serialize) -> Result<Value, Report> {
      json_write_str(value, JsonPretty(false)).and_then(|json| json_read_str(&json))
    }

    struct AuspiceSchemaRetriever {
      schemas: BTreeMap<String, Value>,
    }

    impl AuspiceSchemaRetriever {
      fn new() -> Result<Self, Report> {
        Ok(Self {
          schemas: [AUSPICE_CONFIG_SCHEMA, ANNOTATIONS_SCHEMA, ROOT_SEQUENCE_SCHEMA]
            .into_iter()
            .map(|schema| {
              let schema: Value = json_read_str(schema)?;
              let id = schema["$id"]
                .as_str()
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Vendored schema has no $id"))?
                .to_owned();
              Ok((id, schema))
            })
            .collect::<Result<_, Report>>()?,
        })
      }
    }

    impl Retrieve for AuspiceSchemaRetriever {
      fn retrieve(&self, uri: &Uri<String>) -> Result<Value, Box<dyn StdError + Send + Sync>> {
        self
          .schemas
          .get(uri.as_str())
          .cloned()
          .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, format!("Schema not found: {uri}")).into())
      }
    }
  }
}
