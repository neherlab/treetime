#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, and scratch collections"
)]

#[cfg(test)]
pub(super) mod tests {
  use crate::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, TreeSequences, TreeTraits};
  use crate::auspice::{format_number, group_mutations};
  use crate::output_plan::{CommandKind, ResolvedOutputs, TreeWriteKind};
  use crate::tree_output::{tree_view_for_outputs, write_graph_outputs, write_tree_outputs};
  use approx::assert_ulps_eq;
  use eyre::{Report, WrapErr};
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::Value;
  use serde_json::json;
  use std::collections::BTreeMap;
  use std::iter::once;
  use std::path::PathBuf;
  use tempfile::TempDir;
  use treetime::progress::NoopProgress;
  use treetime::seq::mutation::{Mutation, MutationTrack, Sub};
  use treetime_graph::graph::Graph;
  use treetime_io::auspice_types::AuspiceGenomeAnnotationNuc;
  use treetime_io::nwk::{NwkStyle, nwk_read};
  use treetime_primitives::Seq;
  use treetime_utils::io::fs::read_file_to_string;
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_tree_output_ancestral_models_preserve_semantics() -> Result<(), Report> {
    let setup = helpers::ancestral_setup(helpers::Mutations::NucleotideSubstitution)?;
    let auspice = helpers::auspice(&helpers::ancestral_graph(&setup), CommandKind::Ancestral)?;
    let child = helpers::auspice_child(&auspice, "A");
    assert_eq!(Some(helpers::UPDATED), auspice.data.meta.updated.as_deref());
    assert_eq!(vec![o!("tree"), o!("entropy")], auspice.data.meta.panels);
    assert_eq!(
      Some(3),
      auspice
        .data
        .meta
        .genome_annotations
        .as_ref()
        .and_then(|annotations| annotations.nuc.as_ref())
        .map(|nuc| nuc.end)
    );
    assert_eq!(Some(0.5), child.node_attrs.div);
    assert_eq!(vec![o!("A1T")], child.branch_attrs.mutations["nuc"]);

    let mat = helpers::mat(&helpers::ancestral_graph(&setup))?;
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
  fn test_tree_output_auspice_drops_nucleotide_indel_and_encodes_amino_acid() -> Result<(), Report> {
    let setup = helpers::ancestral_setup(helpers::Mutations::Indel)?;
    let auspice = helpers::auspice(&helpers::ancestral_graph(&setup), CommandKind::Ancestral)?;
    let child = helpers::auspice_child(&auspice, "A");
    assert!(!child.branch_attrs.mutations.contains_key("nuc"));

    let setup = helpers::ancestral_setup(helpers::Mutations::AminoAcid)?;
    let auspice = helpers::auspice(&helpers::ancestral_graph(&setup), CommandKind::Ancestral)?;
    let child = helpers::auspice_child(&auspice, "A");
    assert_eq!(vec![o!("A2T")], child.branch_attrs.mutations["S"]);
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
  fn test_tree_output_auspice_omits_amino_acid_tracks_without_mutations() -> Result<(), Report> {
    let mut setup = helpers::ancestral_setup(helpers::Mutations::AminoAcid)?;
    let a_key = helpers::node_key(&setup.topology, "A");
    setup
      .aa_node_data
      .as_mut()
      .expect("AA fixture must have amino-acid data")
      .node_aa_mutations
      .get_mut(&a_key)
      .expect("AA fixture must have mutations on A")
      .insert(o!("E"), vec![]);

    let auspice = helpers::auspice(&helpers::ancestral_graph(&setup), CommandKind::Ancestral)?;

    let expected = btreemap! { o!("S") => vec![o!("A2T")] };
    assert_eq!(expected, helpers::auspice_child(&auspice, "A").branch_attrs.mutations);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nucleotide_substitution(helpers::Mutations::NucleotideSubstitution, true)]
  #[case::amino_acid_only(        helpers::Mutations::AminoAcid,              true)]
  #[case::nucleotide_indel_only(  helpers::Mutations::Indel,                  false)]
  #[case::no_mutation(            helpers::Mutations::None,                   false)]
  #[trace]
  fn test_tree_output_auspice_lists_genotype_coloring_only_when_a_branch_shows_a_mutation(
    #[case] mutations: helpers::Mutations,
    #[case] expected: bool,
  ) -> Result<(), Report> {
    let setup = helpers::ancestral_setup(mutations)?;

    let auspice = helpers::auspice(&helpers::ancestral_graph(&setup), CommandKind::Ancestral)?;

    let actual = auspice.data.meta.colorings.iter().any(|coloring| coloring.key == "gt");
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_orders_colorings_date_excluded_trait_genotype() -> Result<(), Report> {
    let dated = helpers::dated_setup()?;
    let mugration = helpers::mugration_setup()?;
    let root_sequence = Seq::try_from_str("ACGT")?;
    let edge_mutations = helpers::mutation_on(&dated.topology, "A")?;
    let graph = AnnotatedGraph {
      sequences: Some(TreeSequences {
        root_sequence: &root_sequence,
        edge_mutations: &edge_mutations,
        mutation_counts: None,
        amino_acids: None,
      }),
      traits: Some(TreeTraits {
        attribute: "country",
        states: &mugration.states,
        values: &mugration.values,
        profiles: &mugration.profiles,
      }),
      ..helpers::dated_graph(&dated, None)
    };

    let auspice = helpers::auspice(&graph, CommandKind::Timetree)?;

    let colorings: Vec<(&str, &str)> = auspice
      .data
      .meta
      .colorings
      .iter()
      .map(|coloring| (coloring.key.as_str(), coloring.title.as_str()))
      .collect();
    let expected = vec![
      ("num_date", "Date"),
      ("bad_branch", "Excluded"),
      ("country", "country"),
      ("gt", "Genotype"),
    ];
    assert_eq!(expected, colorings);
    assert_eq!(vec![o!("bad_branch"), o!("country")], auspice.data.meta.filters);
    assert_eq!(Some(o!("country")), auspice.data.meta.display_defaults.color_by);
    Ok(())
  }

  #[rstest]
  #[trace]
  fn test_tree_output_every_fact_combination_gives_valid_auspice_and_newick(
    #[values(false, true)] with_sequences: bool,
    #[values(false, true)] with_dates: bool,
    #[values(false, true)] with_traits: bool,
  ) -> Result<(), Report> {
    let dated = helpers::dated_setup()?;
    let mugration = helpers::mugration_setup()?;
    let root_sequence = Seq::try_from_str("ACGT")?;
    let edge_mutations = helpers::mutation_on(&dated.topology, "A")?;
    let dates = helpers::dated_graph(&dated, None).dates;
    let graph = AnnotatedGraph {
      sequences: with_sequences.then_some(TreeSequences {
        root_sequence: &root_sequence,
        edge_mutations: &edge_mutations,
        mutation_counts: None,
        amino_acids: None,
      }),
      dates: if with_dates { dates } else { None },
      traits: with_traits.then_some(TreeTraits {
        attribute: "country",
        states: &mugration.states,
        values: &mugration.values,
        profiles: &mugration.profiles,
      }),
      ..helpers::annotated(&dated.topology)
    };
    let dir = TempDir::new().wrap_err("When creating a temporary directory")?;
    let auspice_path = dir.path().join("tree.auspice.json");
    let nwk_path = dir.path().join("tree.nwk");
    let outputs = btreemap! {
      TreeWriteKind::Auspice => auspice_path.clone(),
      TreeWriteKind::Nwk(NwkStyle::Beast) => nwk_path.clone(),
    };

    write_tree_outputs(
      &AnnotatedTreeView::new(&graph)?,
      &outputs,
      CommandKind::Timetree,
      &NoopProgress,
    )?;

    let document: Value = json_read_file(&auspice_path)?;
    let errors = helpers::auspice_validator()?
      .iter_errors(&document)
      .map(|error| error.to_string())
      .collect::<Vec<_>>()
      .join("\n");
    assert!(errors.is_empty(), "Auspice schema errors:\n{errors}");
    let reparsed = nwk_read(read_file_to_string(&nwk_path)?.as_bytes())?;
    assert_eq!(dated.topology.names, reparsed.names());
    Ok(())
  }

  #[test]
  fn test_tree_output_conversion_failure_does_not_create_target_and_keeps_prior_file() -> Result<(), Report> {
    let mut setup = helpers::ancestral_setup(helpers::Mutations::NucleotideSubstitution)?;
    setup.root_sequence = Seq::try_from_str("NCG")?;
    let dir = TempDir::new().wrap_err("When creating a temporary directory")?;
    let nwk_path = dir.path().join("tree.nwk");
    let mat_path = dir.path().join("tree.mat.json");
    let outputs = btreemap! {
      TreeWriteKind::Nwk(NwkStyle::Plain) => nwk_path.clone(),
      TreeWriteKind::MatJson => mat_path.clone(),
    };
    let graph = helpers::ancestral_graph(&setup);

    assert_error!(
      write_tree_outputs(
        &AnnotatedTreeView::new(&graph)?,
        &outputs,
        CommandKind::Ancestral,
        &NoopProgress
      ),
      "When writing the tree outputs of ancestral: Node 'A' has root reference nucleotide 'N', but UShER MAT accepts only A, C, G, or T"
    );
    assert!(nwk_path.is_file());
    assert!(!mat_path.exists());

    Ok(())
  }

  #[test]
  fn test_tree_output_graph_json_dumps_topology() -> Result<(), Report> {
    let setup = helpers::ancestral_setup(helpers::Mutations::None)?;
    let dir = TempDir::new().wrap_err("When creating a temporary directory")?;
    let path = dir.path().join("graph.json");
    let outputs = btreemap! { TreeWriteKind::GraphJson => path.clone() };

    write_graph_outputs(&helpers::ancestral_graph(&setup), &outputs)?;

    let actual: Value = json_read_file(&path)?;
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
  fn test_tree_output_graph_outputs_are_written_for_a_graph_that_is_not_a_tree() -> Result<(), Report> {
    let mut graph = Graph::new();
    let [root, left, right, shared] = [graph.add_node(), graph.add_node(), graph.add_node(), graph.add_node()];
    let edges = [
      graph.add_edge(root, left)?,
      graph.add_edge(root, right)?,
      graph.add_edge(left, shared)?,
      graph.add_edge(right, shared)?,
    ];
    graph.build()?;
    let topology = helpers::Topology {
      graph,
      names: btreemap! {
        root => Some(o!("root")),
        left => Some(o!("L")),
        right => Some(o!("R")),
        shared => Some(o!("S")),
      },
      branch_lengths: edges.iter().map(|&edge| (edge, Some(0.1))).collect(),
    };
    let dir = TempDir::new().wrap_err("When creating a temporary directory")?;
    let dot_path = dir.path().join("graph.dot");
    let json_path = dir.path().join("graph.json");
    let nwk_path = dir.path().join("tree.nwk");
    let outputs = ResolvedOutputs {
      tree_outputs: btreemap! {
        TreeWriteKind::Dot => dot_path.clone(),
        TreeWriteKind::GraphJson => json_path.clone(),
        TreeWriteKind::Nwk(NwkStyle::Plain) => nwk_path.clone(),
      },
      non_tree_outputs: btreemap! {},
    };
    let annotated = helpers::annotated(&topology);

    write_graph_outputs(&annotated, &outputs.tree_outputs)?;
    let result = tree_view_for_outputs(&annotated, &outputs);

    assert_error!(
      result,
      format!(
        "These outputs need a tree and were not written: '{}': The graph is not a tree: node {shared} has 2 parents, but a node of a tree has at most one",
        nwk_path.display()
      )
    );
    assert!(dot_path.is_file());
    assert!(json_path.is_file());
    assert!(!nwk_path.exists());
    Ok(())
  }

  #[test]
  fn test_tree_output_tree_view_is_not_built_without_tree_based_outputs() -> Result<(), Report> {
    let topology = helpers::topology()?;
    let outputs = ResolvedOutputs {
      tree_outputs: btreemap! { TreeWriteKind::Dot => PathBuf::from("tree.dot") },
      non_tree_outputs: btreemap! {},
    };
    let annotated = helpers::annotated(&topology);

    let view = tree_view_for_outputs(&annotated, &outputs)?;

    assert!(view.is_none());
    Ok(())
  }

  #[test]
  fn test_tree_output_all_auspice_models_match_augur_v2_schema() -> Result<(), Report> {
    let documents = helpers::all_auspice_documents()?;
    let validator = helpers::auspice_validator()?;

    for (command, document) in &documents {
      let errors = validator
        .iter_errors(document)
        .map(|error| error.to_string())
        .collect::<Vec<_>>()
        .join("\n");
      assert!(errors.is_empty(), "{command} Auspice schema errors:\n{errors}");
    }

    let mut malformed = documents["ancestral"].clone();
    malformed["meta"]
      .as_object_mut()
      .expect("Auspice meta must be an object")
      .remove("updated");
    assert!(!validator.is_valid(&malformed));

    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_writes_a_trait_named_confidence_like_any_trait() -> Result<(), Report> {
    let setup = helpers::mugration_setup()?;
    let country = helpers::auspice(&helpers::mugration_graph(&setup, "country"), CommandKind::Mugration)?;
    let confidence = helpers::auspice(&helpers::mugration_graph(&setup, "confidence"), CommandKind::Mugration)?;

    let expected = json!({ "confidence": helpers::auspice_child(&country, "A").node_attrs.other["country"] });
    let actual = &helpers::auspice_child(&confidence, "A").node_attrs.other;
    assert_eq!(&expected, actual);
    Ok(())
  }

  #[test]
  fn test_tree_output_errors_name_the_command() -> Result<(), Report> {
    let mut topology = helpers::topology()?;
    helpers::set_branch_length(&mut topology, "A", Some(f64::NAN))?;
    let root_sequence = Seq::try_from_str("ACGT")?;
    let edge_mutations = helpers::no_mutations(&topology);
    let graph = helpers::optimize_graph(&topology, &root_sequence, &edge_mutations);
    let outputs = btreemap! { TreeWriteKind::Auspice => PathBuf::from("tree.auspice.json") };

    assert_error!(
      write_tree_outputs(
        &AnnotatedTreeView::new(&graph)?,
        &outputs,
        CommandKind::Optimize,
        &NoopProgress
      ),
      "When writing the tree outputs of optimize: Node 'A' has non-finite div=NaN"
    );
    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_root_sequence_adds_nuc_annotation_and_entropy_panel() -> Result<(), Report> {
    let setup = helpers::dated_setup()?;
    let root_sequence = Seq::try_from_str("ACGT")?;
    let auspice = helpers::auspice(
      &helpers::dated_graph(&setup, Some(&root_sequence)),
      CommandKind::Timetree,
    )?;

    assert_eq!(vec![o!("tree"), o!("entropy")], auspice.data.meta.panels);
    let annotations = auspice
      .data
      .meta
      .genome_annotations
      .as_ref()
      .expect("a root sequence must produce genome annotations");
    assert_eq!(
      Some(AuspiceGenomeAnnotationNuc {
        start: 1,
        end: 4,
        strand: Some(o!("+")),
        r#type: Some(o!("source")),
        other: Value::default(),
      }),
      annotations.nuc
    );
    assert!(annotations.cdses.is_empty());
    let document = helpers::json_value(&auspice)?;
    let errors = helpers::auspice_validator()?
      .iter_errors(&document)
      .map(|error| error.to_string())
      .collect::<Vec<_>>()
      .join("\n");
    assert!(errors.is_empty(), "timetree Auspice schema errors:\n{errors}");

    Ok(())
  }

  #[rstest]
  #[case::prune("prune")]
  #[case::clock("clock")]
  #[case::mugration("mugration")]
  #[trace]
  fn test_tree_output_auspice_without_sequences_has_tree_panel_only(#[case] command: &str) -> Result<(), Report> {
    let document = &helpers::all_auspice_documents()?[command];

    assert_eq!(json!(["tree"]), document["meta"]["panels"]);
    assert_eq!(None, document["meta"].get("genome_annotations"));

    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_rejects_node_without_divergence_or_date() -> Result<(), Report> {
    let mut topology = helpers::topology()?;
    helpers::set_branch_length(&mut topology, "A", None)?;
    let root_sequence = Seq::try_from_str("ACGT")?;
    let edge_mutations = helpers::no_mutations(&topology);
    assert_error!(
      helpers::auspice(
        &helpers::optimize_graph(&topology, &root_sequence, &edge_mutations),
        CommandKind::Optimize
      ),
      "Auspice v2 node 'A' requires divergence or numerical date data"
    );
    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_divergence_sums_from_the_root() -> Result<(), Report> {
    let topology = helpers::topology_from("((A:0.125,B:0.25)AB:0.5,C:0.25)root;")?;

    let auspice = helpers::auspice(&helpers::annotated(&topology), CommandKind::Ancestral)?;

    let ab = helpers::auspice_child(&auspice, "AB");
    let actual: BTreeMap<String, Option<f64>> = once(ab)
      .chain(&ab.children)
      .map(|node| (node.name.clone(), node.node_attrs.div))
      .collect();
    let expected = btreemap! { o!("AB") => Some(0.5), o!("A") => Some(0.625), o!("B") => Some(0.75) };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_rejects_invalid_amino_acid_track_name() -> Result<(), Report> {
    let setup = helpers::ancestral_setup(helpers::Mutations::IndelAndAminoAcid)?;
    assert_error!(
      helpers::auspice(&helpers::ancestral_graph(&setup), CommandKind::Ancestral),
      "Auspice v2 cannot represent amino-acid mutation track 'S/1:weird'"
    );
    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_marks_excluded_nodes() -> Result<(), Report> {
    let mut setup = helpers::dated_setup()?;
    setup.excluded = btreeset! { helpers::node_key(&setup.topology, "B") };

    let auspice = helpers::auspice(&helpers::dated_graph(&setup, None), CommandKind::Clock)?;

    let actual: BTreeMap<&str, &str> = auspice
      .tree
      .children
      .iter()
      .map(|child| {
        let bad_branch = child
          .node_attrs
          .bad_branch
          .as_ref()
          .expect("dated nodes carry bad_branch");
        (child.name.as_str(), bad_branch.value())
      })
      .collect();
    assert_eq!(btreemap! { "A" => "No", "B" => "Yes", "C" => "No" }, actual);
    Ok(())
  }

  #[test]
  fn test_tree_output_format_number_fractional_precision() {
    assert_ulps_eq!(0.123457, format_number(0.12345678, 6), max_ulps = 0);
    assert_ulps_eq!(123.456789, format_number(123.456789, 6), max_ulps = 0);
    assert_ulps_eq!(0.0, format_number(0.0, 6), max_ulps = 0);
    assert_ulps_eq!(2020.123, format_number(2020.1234567, 3), max_ulps = 0);
  }

  #[test]
  fn test_tree_output_group_mutations_drops_nucleotide_indels_keeps_amino_acid_indels() -> Result<(), Report> {
    let mutations = vec![
      Mutation::substitution(
        MutationTrack::Nucleotide,
        Sub::new(helpers::c(b'A'), 0_usize, helpers::c(b'T'))?,
      ),
      helpers::deletion(MutationTrack::Nucleotide, (1, 3), "CG")?,
      Mutation::substitution(
        MutationTrack::AminoAcid(o!("GENE")),
        Sub::new(helpers::c(b'K'), 4_usize, helpers::c(b'R'))?,
      ),
      helpers::deletion(MutationTrack::AminoAcid(o!("GENE")), (1, 3), "CG")?,
    ];

    let grouped = group_mutations(&mutations)?;

    let expected = btreemap! {
      o!("GENE") => vec![o!("K5R"), o!("C2-"), o!("G3-")],
      o!("nuc") => vec![o!("A1T")],
    };
    assert_eq!(expected, grouped);
    Ok(())
  }

  pub(crate) mod helpers {
    use crate::annotated_graph::{
      AnnotatedGraph, AnnotatedTreeView, Divergence, TreeAminoAcids, TreeDates, TreeSequences, TreeTraits,
    };
    use crate::auspice::auspice_tree;
    use crate::output_plan::CommandKind;
    use crate::usher_mat::mat_tree;
    use eyre::Report;
    use jsonschema::{Retrieve, Uri, Validator};
    use maplit::btreemap;
    use ndarray::{Array1, array};
    use serde::Serialize;
    use serde_json::Value;
    use std::collections::{BTreeMap, BTreeSet};
    use std::error::Error;
    use std::io;
    use treetime::ancestral::aa::AaNodeData;
    use treetime::partition::storage::discrete::DiscreteStates;
    use treetime::seq::mutation::{AlignedMutation, Mutation, MutationEvent, MutationTrack, Sub};
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::auspice_types::{AuspiceTree, AuspiceTreeNode};
    use treetime_io::nwk::nwk_read;
    use treetime_io::usher_mat::UsherTree;
    use treetime_primitives::{AsciiChar, Seq};
    use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
    use treetime_utils::{make_report, o};
    use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

    pub(crate) const UPDATED: &str = "2026-07-19";
    const AUSPICE_SCHEMA: &str = include_str!("schemas/auspice/schema-export-v2.json");
    const AUSPICE_CONFIG_SCHEMA: &str = include_str!("schemas/auspice/schema-auspice-config-v2.json");
    const ANNOTATIONS_SCHEMA: &str = include_str!("schemas/auspice/schema-annotations.json");
    const ROOT_SEQUENCE_SCHEMA: &str = include_str!("schemas/auspice/schema-export-root-sequence.json");
    const MODEL_TREE: &str = "(A:0.1,B:0,C:0.5)root;";

    #[derive(Clone, Copy, Debug)]
    pub(crate) enum Mutations {
      None,
      NucleotideSubstitution,
      NucleotideSubstitutionAndAminoAcid,
      Indel,
      AminoAcid,
      IndelAndAminoAcid,
    }

    pub(crate) struct Topology {
      pub(crate) graph: Graph,
      pub(crate) names: BTreeMap<GraphNodeKey, Option<String>>,
      pub(crate) branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
    }

    pub(crate) struct AncestralSetup {
      pub(crate) topology: Topology,
      pub(crate) root_sequence: Seq,
      pub(crate) edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
      pub(crate) aa_node_data: Option<AaNodeData>,
      pub(crate) aa_annotations: BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
    }

    pub(crate) struct MugrationSetup {
      pub(crate) topology: Topology,
      pub(crate) states: DiscreteStates,
      pub(crate) values: BTreeMap<GraphNodeKey, Option<String>>,
      pub(crate) profiles: BTreeMap<GraphNodeKey, Option<Array1<f64>>>,
    }

    pub(crate) struct DatedSetup {
      pub(crate) topology: Topology,
      pub(crate) div: BTreeMap<GraphNodeKey, f64>,
      pub(crate) num_date: BTreeMap<GraphNodeKey, Option<f64>>,
      pub(crate) excluded: BTreeSet<GraphNodeKey>,
      pub(crate) edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
    }

    pub(crate) fn topology() -> Result<Topology, Report> {
      topology_from(MODEL_TREE)
    }

    pub(crate) fn topology_from(nwk: &str) -> Result<Topology, Report> {
      let parsed = nwk_read(nwk.as_bytes())?;
      Ok(Topology {
        names: parsed.names(),
        graph: parsed.graph,
        branch_lengths: parsed.branch_lengths,
      })
    }

    pub(crate) fn annotated(topology: &Topology) -> AnnotatedGraph<'_> {
      AnnotatedGraph {
        graph: &topology.graph,
        names: &topology.names,
        divergence_branch_lengths: &topology.branch_lengths,
        time_branch_lengths: None,
        divergence: Divergence::CumulativeBranchLength,
        sequences: None,
        dates: None,
        traits: None,
      }
    }

    pub(crate) fn auspice(graph: &AnnotatedGraph<'_>, command: CommandKind) -> Result<AuspiceTree, Report> {
      auspice_tree(&AnnotatedTreeView::new(graph)?, command, UPDATED)
    }

    pub(crate) fn mat(graph: &AnnotatedGraph<'_>) -> Result<UsherTree, Report> {
      Ok(mat_tree(&AnnotatedTreeView::new(graph)?)?.tree)
    }

    pub(crate) fn no_mutations(topology: &Topology) -> BTreeMap<GraphEdgeKey, Vec<Mutation>> {
      topology.graph.get_edges().map(|edge| (edge.key(), vec![])).collect()
    }

    pub(crate) fn mutation_on(
      topology: &Topology,
      name: &str,
    ) -> Result<BTreeMap<GraphEdgeKey, Vec<Mutation>>, Report> {
      let mut edge_mutations = no_mutations(topology);
      edge_mutations.insert(parent_edge(topology, name)?, vec![substitution(b'A', 0, b'T')?]);
      Ok(edge_mutations)
    }

    pub(crate) fn ancestral_setup(mutations: Mutations) -> Result<AncestralSetup, Report> {
      let topology = topology_from("(A:0.5,B:0)root;")?;
      let a_edge = parent_edge(&topology, "A")?;
      let b_edge = parent_edge(&topology, "B")?;
      let a_key = node_key(&topology, "A");

      let mut a_mutations = Vec::new();
      if matches!(
        mutations,
        Mutations::NucleotideSubstitution | Mutations::NucleotideSubstitutionAndAminoAcid
      ) {
        a_mutations.push(substitution(b'A', 0, b'T')?);
      }
      if matches!(mutations, Mutations::Indel | Mutations::IndelAndAminoAcid) {
        a_mutations.push(deletion(MutationTrack::Nucleotide, (1, 3), "CG")?);
      }

      let include_aa = matches!(
        mutations,
        Mutations::AminoAcid | Mutations::NucleotideSubstitutionAndAminoAcid | Mutations::IndelAndAminoAcid
      );
      let aa_node_data = include_aa.then(|| {
        let track = if matches!(mutations, Mutations::IndelAndAminoAcid) {
          "S/1:weird"
        } else {
          "S"
        };
        let mut aa = AaNodeData::default();
        aa.root_aa_sequences.insert(o!(track), o!("AA"));
        aa.node_aa_mutations.insert(
          a_key,
          btreemap! {
            o!(track) => vec![MutationEvent::Substitution(Sub::new(c(b'A'), 1_usize, c(b'T')).unwrap())],
          },
        );
        aa
      });
      let aa_annotations = if include_aa {
        btreemap! {
          o!("S") => AugurNodeDataJsonAnnotationEntry {
            start: Some(1),
            end: Some(3),
            strand: Some(o!("+")),
            entry_type: Some(o!("CDS")),
            ..AugurNodeDataJsonAnnotationEntry::default()
          },
        }
      } else {
        btreemap! {}
      };

      Ok(AncestralSetup {
        topology,
        root_sequence: Seq::try_from_str("ACG")?,
        edge_mutations: btreemap! { a_edge => a_mutations, b_edge => vec![] },
        aa_node_data,
        aa_annotations,
      })
    }

    pub(crate) fn ancestral_graph(setup: &AncestralSetup) -> AnnotatedGraph<'_> {
      AnnotatedGraph {
        sequences: Some(TreeSequences {
          root_sequence: &setup.root_sequence,
          edge_mutations: &setup.edge_mutations,
          mutation_counts: None,
          amino_acids: setup.aa_node_data.as_ref().map(|node_data| TreeAminoAcids {
            node_data,
            cdses: &setup.aa_annotations,
          }),
        }),
        ..annotated(&setup.topology)
      }
    }

    pub(crate) fn optimize_graph<'a>(
      topology: &'a Topology,
      root_sequence: &'a Seq,
      edge_mutations: &'a BTreeMap<GraphEdgeKey, Vec<Mutation>>,
    ) -> AnnotatedGraph<'a> {
      AnnotatedGraph {
        sequences: Some(TreeSequences {
          root_sequence,
          edge_mutations,
          mutation_counts: None,
          amino_acids: None,
        }),
        ..annotated(topology)
      }
    }

    pub(crate) fn mugration_setup() -> Result<MugrationSetup, Report> {
      let topology = topology()?;
      let (values, profiles) = topology
        .graph
        .get_nodes()
        .enumerate()
        .map(|(index, node)| {
          let key = node.key();
          let (state, profile) = if index % 2 == 0 {
            ("CH", array![1.0, 0.0])
          } else {
            ("US", array![0.0, 1.0])
          };
          ((key, Some(o!(state))), (key, Some(profile)))
        })
        .unzip();
      Ok(MugrationSetup {
        topology,
        states: DiscreteStates::from_values(["CH", "US"].into_iter(), "?"),
        values,
        profiles,
      })
    }

    pub(crate) fn mugration_graph<'a>(setup: &'a MugrationSetup, attribute: &'a str) -> AnnotatedGraph<'a> {
      AnnotatedGraph {
        traits: Some(TreeTraits {
          attribute,
          states: &setup.states,
          values: &setup.values,
          profiles: &setup.profiles,
        }),
        ..annotated(&setup.topology)
      }
    }

    pub(crate) fn dated_setup() -> Result<DatedSetup, Report> {
      let topology = topology()?;
      let div = topology
        .graph
        .get_nodes()
        .enumerate()
        .map(|(index, node)| (node.key(), index as f64 / 2.0))
        .collect();
      let num_date = topology
        .graph
        .get_nodes()
        .enumerate()
        .map(|(index, node)| (node.key(), Some(2020.0 + index as f64)))
        .collect();
      let edge_mutations = no_mutations(&topology);
      Ok(DatedSetup {
        topology,
        div,
        num_date,
        excluded: BTreeSet::new(),
        edge_mutations,
      })
    }

    pub(crate) fn dated_graph<'a>(setup: &'a DatedSetup, root_sequence: Option<&'a Seq>) -> AnnotatedGraph<'a> {
      AnnotatedGraph {
        divergence: Divergence::Values(&setup.div),
        sequences: root_sequence.map(|root_sequence| TreeSequences {
          root_sequence,
          edge_mutations: &setup.edge_mutations,
          mutation_counts: None,
          amino_acids: None,
        }),
        dates: Some(TreeDates {
          num_date: &setup.num_date,
          confidence: None,
          excluded: &setup.excluded,
          input_dates: None,
        }),
        ..annotated(&setup.topology)
      }
    }

    pub(crate) fn all_auspice_documents() -> Result<BTreeMap<&'static str, Value>, Report> {
      let ancestral = ancestral_setup(Mutations::NucleotideSubstitution)?;
      let ancestral = auspice(&ancestral_graph(&ancestral), CommandKind::Ancestral)?;

      let optimize = topology()?;
      let optimize_reference = Seq::try_from_str("ACGT")?;
      let optimize_mutations = no_mutations(&optimize);
      let optimize = auspice(
        &optimize_graph(&optimize, &optimize_reference, &optimize_mutations),
        CommandKind::Optimize,
      )?;

      let prune = topology()?;
      let prune = auspice(&annotated(&prune), CommandKind::Prune)?;

      let clock = dated_setup()?;
      let clock = auspice(&dated_graph(&clock, None), CommandKind::Clock)?;

      let mugration = mugration_setup()?;
      let mugration = auspice(&mugration_graph(&mugration, "country"), CommandKind::Mugration)?;

      let timetree = dated_setup()?;
      let timetree = auspice(&dated_graph(&timetree, None), CommandKind::Timetree)?;

      Ok(btreemap! {
        "ancestral" => json_value(&ancestral)?,
        "optimize" => json_value(&optimize)?,
        "prune" => json_value(&prune)?,
        "clock" => json_value(&clock)?,
        "mugration" => json_value(&mugration)?,
        "timetree" => json_value(&timetree)?,
      })
    }

    pub(crate) fn all_mat_documents() -> Result<Vec<UsherTree>, Report> {
      let mut topology = topology()?;
      set_branch_length(&mut topology, "A", None)?;
      set_branch_length(&mut topology, "B", Some(0.0))?;
      set_branch_length(&mut topology, "C", Some(0.5))?;
      let no_mutations = no_mutations(&topology);
      let ancestral_reference = Seq::try_from_str("ACG")?;
      let optimize_reference = Seq::try_from_str("ACGT")?;
      [
        Some(&ancestral_reference),
        Some(&optimize_reference),
        None,
        None,
        None,
        None,
      ]
      .into_iter()
      .map(|root_sequence| {
        mat(&AnnotatedGraph {
          sequences: root_sequence.map(|root_sequence| TreeSequences {
            root_sequence,
            edge_mutations: &no_mutations,
            mutation_counts: None,
            amino_acids: None,
          }),
          ..annotated(&topology)
        })
      })
      .collect()
    }

    pub(crate) fn auspice_validator() -> Result<Validator, Report> {
      let schema = json_read_str(AUSPICE_SCHEMA)?;
      Ok(
        jsonschema::draft6::options()
          .with_retriever(AuspiceSchemaRetriever::new()?)
          .build(&schema)?,
      )
    }

    pub(crate) fn auspice_child<'a>(tree: &'a AuspiceTree, name: &str) -> &'a AuspiceTreeNode {
      tree
        .tree
        .children
        .iter()
        .find(|child| child.name == name)
        .expect("fixture child must exist")
    }

    pub(crate) fn branch_length(topology: &Topology, target_name: &str) -> Result<Option<f64>, Report> {
      Ok(topology.branch_lengths[&parent_edge(topology, target_name)?])
    }

    pub(crate) fn set_branch_length(topology: &mut Topology, name: &str, length: Option<f64>) -> Result<(), Report> {
      let edge_key = parent_edge(topology, name)?;
      topology.branch_lengths.insert(edge_key, length);
      Ok(())
    }

    pub(crate) fn node_key(topology: &Topology, name: &str) -> GraphNodeKey {
      topology
        .names
        .iter()
        .find_map(|(key, node_name)| (node_name.as_deref() == Some(name)).then_some(*key))
        .expect("fixture node must exist")
    }

    pub(crate) fn parent_edge(topology: &Topology, name: &str) -> Result<GraphEdgeKey, Report> {
      Ok(
        topology
          .graph
          .node_parent(node_key(topology, name))?
          .expect("fixture node must have a parent")
          .1,
      )
    }

    pub(crate) fn c(value: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(value)
    }

    pub(crate) fn substitution(reff: u8, pos: usize, qry: u8) -> Result<Mutation, Report> {
      Ok(Mutation::substitution(
        MutationTrack::Nucleotide,
        Sub::new(c(reff), pos, c(qry))?,
      ))
    }

    pub(crate) fn deletion(track: MutationTrack, range: (usize, usize), sequence: &str) -> Result<Mutation, Report> {
      Ok(Mutation {
        track,
        event: MutationEvent::Deletion(AlignedMutation {
          range,
          sequence: Seq::try_from_str(sequence)?,
        }),
      })
    }

    pub(crate) fn json_value(value: &impl Serialize) -> Result<Value, Report> {
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
                .ok_or_else(|| make_report!("Vendored schema has no $id"))?
                .to_owned();
              Ok((id, schema))
            })
            .collect::<Result<_, Report>>()?,
        })
      }
    }

    impl Retrieve for AuspiceSchemaRetriever {
      fn retrieve(&self, uri: &Uri<String>) -> Result<Value, Box<dyn Error + Send + Sync>> {
        self
          .schemas
          .get(uri.as_str())
          .cloned()
          .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, format!("Schema not found: {uri}")).into())
      }
    }
  }
}
