#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::pipeline::AncestralPartition;
  use crate::commands::ancestral::aa_node_data::AaNodeData;
  use crate::commands::ancestral::result::AncestralGraphData;
  use crate::commands::shared::tree_output::{
    ancestral_to_auspice, ancestral_to_mat, ancestral_to_phyloxml, clock_to_auspice, clock_to_mat, clock_to_phyloxml,
    format_number, group_mutations, mat_mutation, mugration_to_auspice, mugration_to_mat, mugration_to_phyloxml,
    optimize_to_auspice,
    optimize_to_mat, optimize_to_phyloxml, prune_to_auspice, prune_to_mat, prune_to_phyloxml, timetree_to_auspice,
    timetree_to_mat, timetree_to_phyloxml, write_ancestral_tree_outputs,
  };
  use crate::gtr::get_gtr::GtrModelName;
  use crate::partition::fitch::PartitionFitch;
  use crate::partition::sparse::{SparseEdgePartition, SparseNodePartition};
  use crate::partition::traits::BranchTopology;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::indel::InDel;
  use crate::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub};
  use approx::assert_ulps_eq;
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::Value;
  use std::sync::Arc;
  use tempfile::TempDir;
  use treetime_graph::node::{GraphNodeKey, Named};
  use treetime_io::graph::TreeWriteKind;
  use treetime_io::nwk::{CommentProviders, NwkStyle, nwk_read_str};
  use treetime_primitives::{AsciiChar, Seq};
  use treetime_utils::io::json::{JsonPretty, json_read_file, json_read_str, json_write_str};

  #[test]
  fn test_tree_output_ancestral_models_preserve_semantics() -> Result<(), Report> {
    let graph = helpers::ancestral_graph(helpers::Mutations::NucleotideSubstitution)?;

    let auspice = ancestral_to_auspice(&graph, "2026-07-19")?;
    let child = helpers::auspice_child(&auspice, "A");
    assert_eq!(Some("2026-07-19"), auspice.data.meta.updated.as_deref());
    assert_eq!(vec!["tree".to_owned()], auspice.data.meta.panels);
    assert_eq!(Some(0.5), child.node_attrs.div);
    assert_eq!(vec!["A1T".to_owned()], child.branch_attrs.mutations["nuc"]);

    let phyloxml = ancestral_to_phyloxml(&graph)?;
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

    let mat = ancestral_to_mat(&graph)?;
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
    let graph = helpers::ancestral_graph(helpers::Mutations::IndelAndAminoAcid)?;

    let phyloxml = ancestral_to_phyloxml(&graph)?;
    let child = helpers::phyloxml_child(&phyloxml, "A");
    let properties = child
      .property
      .iter()
      .map(|property| property.value.as_str())
      .collect::<Vec<_>>();
    assert!(properties.contains(&"nuc:del:2-3:CG"));
    assert!(properties.contains(&"aa:S%2F1%3Aweird:sub:A2T"));

    let graph = helpers::ancestral_graph(helpers::Mutations::Indel)?;
    let auspice = ancestral_to_auspice(&graph, "2026-07-19")?;
    let child = helpers::auspice_child(&auspice, "A");
    assert_eq!(
      vec!["C2-".to_owned(), "G3-".to_owned()],
      child.branch_attrs.mutations["nuc"]
    );

    let graph = helpers::ancestral_graph(helpers::Mutations::AminoAcid)?;
    let auspice = ancestral_to_auspice(&graph, "2026-07-19")?;
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
    let graph = helpers::ancestral_graph(helpers::Mutations::Indel)?;
    let error = ancestral_to_mat(&graph).expect_err("MAT must reject indels");
    assert!(error.to_string().contains("insertion or deletion"));

    let graph = helpers::ancestral_graph(helpers::Mutations::AminoAcid)?;
    let error = ancestral_to_mat(&graph).expect_err("MAT must reject amino-acid mutations");
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
    let graph = helpers::ancestral_graph(helpers::Mutations::Indel)?;
    let dir = TempDir::new()?;
    let nwk_path = dir.path().join("tree.nwk");
    let mat_path = dir.path().join("tree.mat.json");
    let outputs = btreemap! {
      TreeWriteKind::nwk(NwkStyle::Plain) => nwk_path.clone(),
      TreeWriteKind::MatJson => mat_path.clone(),
    };

    let error =
      write_ancestral_tree_outputs(&graph, &outputs, &CommentProviders::new()).expect_err("MAT conversion must fail");
    assert!(error.to_string().contains("insertion or deletion"));
    assert!(nwk_path.is_file());
    assert!(!mat_path.exists());

    Ok(())
  }

  #[test]
  fn test_tree_output_graph_json_dumps_concrete_graph_data() -> Result<(), Report> {
    let graph = helpers::ancestral_graph(helpers::Mutations::None)?;
    let dir = TempDir::new()?;
    let path = dir.path().join("graph.json");
    let outputs = btreemap! { TreeWriteKind::GraphJson => path.clone() };

    write_ancestral_tree_outputs(&graph, &outputs, &CommentProviders::new())?;
    let actual: Value = json_read_file(&path)?;
    assert_eq!(Value::String("JC69".to_owned()), actual["data"]["model_name"]);
    assert_eq!(Value::Array(vec![Value::Bool(false); 3]), actual["data"]["mask"]);
    assert!(actual["data"]["partition"]["Fitch"].is_object());
    assert_eq!(Value::from(0.0), helpers::edge_branch_length(&actual, "B"));

    Ok(())
  }

  #[test]
  fn test_tree_output_mutation_free_mat_needs_no_reference() -> Result<(), Report> {
    let graph = helpers::ancestral_graph_without_partition()?;
    let mat = ancestral_to_mat(&graph)?;
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
    let graph = helpers::ancestral_graph(helpers::Mutations::IndelAndAminoAcid)?;
    let error = ancestral_to_auspice(&graph, "2026-07-19")
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
      let graph: GraphAncestral = nwk_read_str(&document.newick)?;
      assert_eq!(
        None,
        helpers::branch_length(&graph, "A")?,
        "{command}: {}",
        document.newick
      );
      assert_eq!(
        Some(0.0),
        helpers::branch_length(&graph, "B")?,
        "{command}: {}",
        document.newick
      );
      assert_eq!(
        Some(0.5),
        helpers::branch_length(&graph, "C")?,
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
      Mutation::substitution(MutationTrack::Nucleotide, Sub::new(helpers::c(b'A'), 0_usize, helpers::c(b'T'))?),
      Mutation::indel(MutationTrack::Nucleotide, &InDel::del((1, 3), Seq::try_from_str("CG")?)?)?,
      Mutation::substitution(
        MutationTrack::AminoAcid("GENE".to_owned()),
        Sub::new(helpers::c(b'K'), 4_usize, helpers::c(b'R'))?,
      ),
      Mutation::indel(MutationTrack::AminoAcid("GENE".to_owned()), &InDel::del((1, 3), Seq::try_from_str("CG")?)?)?,
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
    use crate::clock::clock_graph::GraphClock;
    use crate::clock::clock_model::ClockModel;
    use crate::commands::clock::run::ClockGraphData;
    use crate::commands::optimize::result::OptimizeGraphData;
    use crate::commands::prune::result::PruneGraphData;
    use crate::commands::timetree::result::TimetreeGraphData;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::gtr::gtr::{GTR, GTRParams};
    use crate::mugration::result::{MugrationGraphData, MugrationResult};
    use crate::partition::dense::{DenseNodePartition, DenseSeqDistribution, DenseSeqInfo};
    use crate::partition::discrete_states::DiscreteStates;
    use crate::partition::marginal_discrete::PartitionMarginalDiscrete;
    use crate::partition::timetree::GraphTimetree;
    use crate::payload::clock_set::ClockSet;
    use jsonschema::{Retrieve, Uri, Validator};
    use ndarray::array;
    use serde::Serialize;
    use std::collections::BTreeMap;
    use std::error::Error as StdError;
    use std::io;
    use treetime_graph::edge::{GraphEdge, HasBranchLength};
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNode;
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

    pub fn ancestral_graph(mutations: Mutations) -> Result<GraphAncestral<AncestralGraphData>, Report> {
      let graph: GraphAncestral = nwk_read_str("(A:0.5,B:0)root;")?;
      let root_key = node_key(&graph, "root");
      let a_key = node_key(&graph, "A");
      let b_key = node_key(&graph, "B");
      graph
        .get_node(a_key)
        .unwrap()
        .write_arc()
        .payload()
        .write_arc()
        .confidence = Some(0.9);
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
        SparseEdgePartition::with_fitch_subs(vec![Sub::new(c(b'A'), 0_usize, c(b'T'))?])
      } else {
        SparseEdgePartition::default()
      };
      if include_indel {
        a_edge_data.indels = vec![InDel::del((1, 3), Seq::try_from_str("CG")?)?];
      }
      let partition = PartitionFitch {
        index: 0,
        alphabet: alphabet.clone(),
        length: 3,
        nodes: btreemap! {
          root_key => SparseNodePartition::new(&root_sequence, &alphabet)?,
          a_key => SparseNodePartition::new(&a_sequence, &alphabet)?,
          b_key => SparseNodePartition::new(&b_sequence, &alphabet)?,
        },
        edges: btreemap! {
          a_edge => a_edge_data,
          b_edge => SparseEdgePartition::default(),
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
      let data = AncestralGraphData::new(
        Some(AncestralPartition::Fitch(Arc::new(RwLock::new(partition)))),
        None,
        GtrModelName::JC69,
        vec![false; 3],
        aa_node_data,
      );
      Ok(graph.map_data(data))
    }

    pub fn ancestral_graph_without_partition() -> Result<GraphAncestral<AncestralGraphData>, Report> {
      let graph: GraphAncestral = nwk_read_str(MODEL_TREE)?;
      Ok(graph.map_data(AncestralGraphData::new(None, None, GtrModelName::JC69, vec![], None)))
    }

    pub fn all_auspice_documents() -> Result<Vec<Value>, Report> {
      let ancestral = ancestral_to_auspice(&ancestral_graph(Mutations::NucleotideSubstitution)?, "2026-07-19")?;
      let optimize = optimize_to_auspice(&optimize_graph()?, "2026-07-19")?;
      let prune = prune_to_auspice(&prune_graph()?, "2026-07-19")?;
      let clock = clock_to_auspice(&clock_graph()?, "2026-07-19")?;
      let mugration = mugration_to_auspice(&mugration_graph()?, "2026-07-19")?;
      let timetree = timetree_to_auspice(&timetree_graph()?, "2026-07-19")?;

      [ancestral, optimize, prune, clock, mugration, timetree]
        .iter()
        .map(json_value)
        .collect()
    }

    pub fn all_phyloxml_documents() -> Result<Vec<Phyloxml>, Report> {
      Ok(vec![
        ancestral_to_phyloxml(&ancestral_graph_without_partition()?)?,
        optimize_to_phyloxml(&optimize_graph()?)?,
        prune_to_phyloxml(&prune_graph()?)?,
        clock_to_phyloxml(&clock_graph()?)?,
        mugration_to_phyloxml(&mugration_graph()?)?,
        timetree_to_phyloxml(&timetree_graph()?)?,
      ])
    }

    pub fn all_mat_documents() -> Result<Vec<UsherTree>, Report> {
      let ancestral = ancestral_graph_without_partition()?;
      set_mat_branch_lengths(&ancestral)?;
      let optimize = optimize_graph()?;
      set_mat_branch_lengths(&optimize)?;
      let prune = prune_graph()?;
      set_mat_branch_lengths(&prune)?;
      let clock = clock_graph()?;
      set_mat_branch_lengths(&clock)?;
      let mugration = mugration_graph()?;
      set_mat_branch_lengths(&mugration)?;
      let timetree = timetree_graph()?;
      set_timetree_mat_branch_lengths(&timetree)?;

      Ok(vec![
        ancestral_to_mat(&ancestral)?,
        optimize_to_mat(&optimize)?,
        prune_to_mat(&prune)?,
        clock_to_mat(&clock)?,
        mugration_to_mat(&mugration)?,
        timetree_to_mat(&timetree)?,
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
      let graph = optimize_graph()?;
      set_branch_length(&graph, "A", None)?;
      optimize_to_auspice(&graph, "2026-07-19")
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

    pub fn edge_branch_length(graph: &Value, target_name: &str) -> Value {
      let target_key = graph["nodes"]
        .as_array()
        .unwrap()
        .iter()
        .position(|node| node["data"]["name"] == target_name)
        .unwrap();
      graph["edges"]
        .as_array()
        .unwrap()
        .iter()
        .find(|edge| edge["target"] == target_key)
        .unwrap()["data"]["branch_length"]
        .clone()
    }

    pub fn branch_length(graph: &GraphAncestral, target_name: &str) -> Result<Option<f64>, Report> {
      let node_key = node_key(graph, target_name);
      let edge_key = graph
        .node_parent(node_key)?
        .expect("fixture target must not be the root")
        .1;
      Ok(
        graph
          .get_edge(edge_key)
          .expect("fixture edge must exist")
          .read_arc()
          .payload()
          .read_arc()
          .branch_length(),
      )
    }

    fn node_key(graph: &GraphAncestral, name: &str) -> GraphNodeKey {
      graph
        .get_nodes()
        .into_iter()
        .find_map(|node| {
          let node = node.read_arc();
          let payload = node.payload().read_arc();
          (payload.name().as_ref().map(AsRef::as_ref) == Some(name)).then(|| node.key())
        })
        .expect("fixture node must exist")
    }

    pub fn c(value: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(value)
    }

    fn fixed_clock_model() -> Result<ClockModel, Report> {
      ClockModel::with_fixed_rate(&ClockSet::leaf_contribution(Some(2020.0)), 1.0)
    }

    fn optimize_graph() -> Result<GraphAncestral<OptimizeGraphData>, Report> {
      let graph: GraphAncestral = nwk_read_str(MODEL_TREE)?;
      Ok(graph.map_data(OptimizeGraphData::new(
        jc69(JC69Params::default())?,
        GtrModelName::JC69,
        vec![],
        vec![],
      )))
    }

    fn prune_graph() -> Result<GraphAncestral<PruneGraphData>, Report> {
      let graph: GraphAncestral = nwk_read_str(MODEL_TREE)?;
      Ok(graph.map_data(PruneGraphData::new(Some(jc69(JC69Params::default())?), vec![])))
    }

    fn clock_graph() -> Result<GraphClock<ClockGraphData>, Report> {
      let graph: GraphClock = nwk_read_str(MODEL_TREE)?;
      for (index, node) in graph.get_nodes().into_iter().enumerate() {
        let node = node.write_arc();
        let mut payload = node.payload().write_arc();
        payload.div = index as f64 / 2.0;
        payload.time = Some(2020.0 + index as f64);
      }
      Ok(graph.map_data(ClockGraphData::new(fixed_clock_model()?, vec![])))
    }

    fn mugration_graph() -> Result<GraphAncestral<MugrationGraphData>, Report> {
      let graph: GraphAncestral = nwk_read_str(MODEL_TREE)?;
      let states = DiscreteStates::from_values(["CH", "US"].into_iter(), "?");
      let gtr = GTR::new(GTRParams {
        n_states: 2,
        mu: 1.0,
        W: None,
        pi: array![0.5, 0.5],
      })?;
      let mut partition = PartitionMarginalDiscrete::new(gtr, states, 1e-8, false);
      partition.data.nodes = graph
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
            DenseNodePartition {
              seq: DenseSeqInfo::default(),
              profile: DenseSeqDistribution::new(profile, 0.0),
            },
          )
        })
        .collect();
      Ok(MugrationResult::new(graph, partition, "country", 0.0).graph)
    }

    fn timetree_graph() -> Result<GraphTimetree<TimetreeGraphData>, Report> {
      let graph: GraphTimetree = nwk_read_str(MODEL_TREE)?;
      for (index, node) in graph.get_nodes().into_iter().enumerate() {
        let node = node.write_arc();
        let mut payload = node.payload().write_arc();
        payload.div = index as f64 / 2.0;
        payload.time = Some(2020.0 + index as f64);
      }
      Ok(graph.map_data(TimetreeGraphData::new(
        fixed_clock_model()?,
        None,
        vec![],
        None,
        None,
        None,
        None,
      )))
    }

    fn set_mat_branch_lengths<N, E, D>(graph: &Graph<N, E, D>) -> Result<(), Report>
    where
      N: GraphNode + Named,
      E: GraphEdge + HasBranchLength,
      D: Send + Sync,
    {
      set_branch_length(graph, "A", None)?;
      set_branch_length(graph, "B", Some(0.0))?;
      set_branch_length(graph, "C", Some(0.5))
    }

    fn set_branch_length<N, E, D>(graph: &Graph<N, E, D>, name: &str, length: Option<f64>) -> Result<(), Report>
    where
      N: GraphNode + Named,
      E: GraphEdge + HasBranchLength,
      D: Send + Sync,
    {
      let key = graph
        .get_nodes()
        .into_iter()
        .find_map(|node| {
          let node = node.read_arc();
          let payload = node.payload().read_arc();
          (payload.name().as_ref().map(AsRef::as_ref) == Some(name)).then(|| node.key())
        })
        .expect("fixture node must exist");
      let edge_key = graph.node_parent(key)?.expect("fixture node must have a parent").1;
      graph
        .get_edge(edge_key)
        .expect("fixture edge must exist")
        .write_arc()
        .payload()
        .write_arc()
        .set_branch_length(length);
      Ok(())
    }

    fn set_timetree_mat_branch_lengths(graph: &GraphTimetree<TimetreeGraphData>) -> Result<(), Report> {
      for (name, length) in [("A", None), ("B", Some(0.0)), ("C", Some(0.5))] {
        let key = graph
          .get_nodes()
          .into_iter()
          .find_map(|node| {
            let node = node.read_arc();
            let payload = node.payload().read_arc();
            (payload.name().as_ref().map(AsRef::as_ref) == Some(name)).then(|| node.key())
          })
          .expect("fixture node must exist");
        let edge_key = graph.node_parent(key)?.expect("fixture node must have a parent").1;
        graph
          .get_edge(edge_key)
          .expect("fixture edge must exist")
          .write_arc()
          .payload()
          .write_arc()
          .time_length = length;
      }
      Ok(())
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
