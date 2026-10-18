#[cfg(test)]
pub(super) mod tests {
  use crate::__tests__::test_support::tests::{
    Mutations, NUC_ALPHABET, Topology, UPDATED, all_auspice_documents, ancestral_graph, ancestral_setup, annotated,
    auspice_child, auspice_tree_of, auspice_validator, c, dated_graph, dated_setup, deletion, mat_tree_of,
    mugration_graph, mugration_setup, mutation_on, no_mutations, node_key, optimize_graph, set_branch_length, topology,
    topology_from,
  };
  use crate::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, TreeSequences, TreeTraits};
  use crate::auspice::{format_number, group_mutations};
  use crate::output_plan::{CommandKind, ResolvedOutputs, TreeWriteKind};
  use crate::tree_output::{tree_view_for_outputs, write_graph_outputs, write_tree_outputs};
  use approx::assert_ulps_eq;
  use deser_value::Map;
  use eyre::{Report, WrapErr};
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
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
  use treetime_utils::io::json::{json_value_read_file, to_json_value};
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_tree_output_ancestral_models_preserve_semantics() -> Result<(), Report> {
    let setup = ancestral_setup(Mutations::NucleotideSubstitution)?;
    let auspice = auspice_tree_of(&ancestral_graph(&setup), CommandKind::Ancestral)?;
    let child = auspice_child(&auspice, "A");
    assert_eq!(Some(UPDATED), auspice.data.meta.updated.as_deref());
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

    let mat = mat_tree_of(&ancestral_graph(&setup))?;
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
    let setup = ancestral_setup(Mutations::Indel)?;
    let auspice = auspice_tree_of(&ancestral_graph(&setup), CommandKind::Ancestral)?;
    let child = auspice_child(&auspice, "A");
    assert!(!child.branch_attrs.mutations.contains_key("nuc"));

    let setup = ancestral_setup(Mutations::AminoAcid)?;
    let auspice = auspice_tree_of(&ancestral_graph(&setup), CommandKind::Ancestral)?;
    let child = auspice_child(&auspice, "A");
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
    let mut setup = ancestral_setup(Mutations::AminoAcid)?;
    let a_key = node_key(&setup.topology, "A");
    setup
      .aa_node_data
      .as_mut()
      .expect("AA fixture must have amino-acid data")
      .node_aa_mutations
      .get_mut(&a_key)
      .expect("AA fixture must have mutations on A")
      .insert(o!("E"), vec![]);

    let auspice = auspice_tree_of(&ancestral_graph(&setup), CommandKind::Ancestral)?;

    let expected = btreemap! { o!("S") => vec![o!("A2T")] };
    assert_eq!(expected, auspice_child(&auspice, "A").branch_attrs.mutations);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nucleotide_substitution(Mutations::NucleotideSubstitution, true)]
  #[case::amino_acid_only(        Mutations::AminoAcid,              true)]
  #[case::nucleotide_indel_only(  Mutations::Indel,                  false)]
  #[case::no_mutation(            Mutations::None,                   false)]
  #[trace]
  fn test_tree_output_auspice_lists_genotype_coloring_only_when_a_branch_shows_a_mutation(
    #[case] mutations: Mutations,
    #[case] expected: bool,
  ) -> Result<(), Report> {
    let setup = ancestral_setup(mutations)?;

    let auspice = auspice_tree_of(&ancestral_graph(&setup), CommandKind::Ancestral)?;

    let actual = auspice.data.meta.colorings.iter().any(|coloring| coloring.key == "gt");
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_orders_colorings_date_excluded_trait_genotype() -> Result<(), Report> {
    let dated = dated_setup()?;
    let mugration = mugration_setup()?;
    let root_sequence = Seq::try_from_str("ACGT")?;
    let edge_mutations = mutation_on(&dated.topology, "A")?;
    let graph = AnnotatedGraph {
      sequences: Some(TreeSequences {
        alphabet: &NUC_ALPHABET,
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
      ..dated_graph(&dated, None)
    };

    let auspice = auspice_tree_of(&graph, CommandKind::Timetree)?;

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
    let dated = dated_setup()?;
    let mugration = mugration_setup()?;
    let root_sequence = Seq::try_from_str("ACGT")?;
    let edge_mutations = mutation_on(&dated.topology, "A")?;
    let dates = dated_graph(&dated, None).dates;
    let graph = AnnotatedGraph {
      sequences: with_sequences.then_some(TreeSequences {
        alphabet: &NUC_ALPHABET,
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
      ..annotated(&dated.topology)
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

    let document = json_value_read_file(&auspice_path)?;
    let errors = auspice_validator()?
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
    let mut setup = ancestral_setup(Mutations::NucleotideSubstitution)?;
    setup.root_sequence = Seq::try_from_str("NCG")?;
    let dir = TempDir::new().wrap_err("When creating a temporary directory")?;
    let nwk_path = dir.path().join("tree.nwk");
    let mat_path = dir.path().join("tree.mat.json");
    let outputs = btreemap! {
      TreeWriteKind::Nwk(NwkStyle::Plain) => nwk_path.clone(),
      TreeWriteKind::MatJson => mat_path.clone(),
    };
    let graph = ancestral_graph(&setup);

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
    let setup = ancestral_setup(Mutations::None)?;
    let dir = TempDir::new().wrap_err("When creating a temporary directory")?;
    let path = dir.path().join("graph.json");
    let outputs = btreemap! { TreeWriteKind::GraphJson => path.clone() };

    write_graph_outputs(&ancestral_graph(&setup), &outputs)?;

    let actual = json_value_read_file(&path)?;
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
    let topology = Topology {
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
      tree_outputs_by_partition: btreemap! {},
      non_tree_outputs_by_partition: btreemap! {},
    };
    let annotated = annotated(&topology);

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
    let topology = topology()?;
    let outputs = ResolvedOutputs {
      tree_outputs: btreemap! { TreeWriteKind::Dot => PathBuf::from("tree.dot") },
      non_tree_outputs: btreemap! {},
      tree_outputs_by_partition: btreemap! {},
      non_tree_outputs_by_partition: btreemap! {},
    };
    let annotated = annotated(&topology);

    let view = tree_view_for_outputs(&annotated, &outputs)?;

    assert!(view.is_none());
    Ok(())
  }

  #[test]
  fn test_tree_output_all_auspice_models_match_augur_v2_schema() -> Result<(), Report> {
    let documents = all_auspice_documents()?;
    let validator = auspice_validator()?;

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
    let setup = mugration_setup()?;
    let country = auspice_tree_of(&mugration_graph(&setup, "country"), CommandKind::Mugration)?;
    let confidence = auspice_tree_of(&mugration_graph(&setup, "confidence"), CommandKind::Mugration)?;

    let country_trait = auspice_child(&country, "A").node_attrs.other.get("country").cloned();
    let expected = Map::from_iter([("confidence", country_trait.expect("the country trait is written"))]);
    let actual = &auspice_child(&confidence, "A").node_attrs.other;
    assert_eq!(&expected, actual);
    Ok(())
  }

  #[test]
  fn test_tree_output_errors_name_the_command() -> Result<(), Report> {
    let mut topology = topology()?;
    set_branch_length(&mut topology, "A", Some(f64::NAN))?;
    let root_sequence = Seq::try_from_str("ACGT")?;
    let edge_mutations = no_mutations(&topology);
    let graph = optimize_graph(&topology, &root_sequence, &edge_mutations);
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
    let setup = dated_setup()?;
    let root_sequence = Seq::try_from_str("ACGT")?;
    let auspice = auspice_tree_of(&dated_graph(&setup, Some(&root_sequence)), CommandKind::Timetree)?;

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
        other: Map::new(),
      }),
      annotations.nuc
    );
    assert!(annotations.cdses.is_empty());
    let document = to_json_value(&auspice)?;
    let errors = auspice_validator()?
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
    let document = &all_auspice_documents()?[command];

    assert_eq!(json!(["tree"]), document["meta"]["panels"]);
    assert_eq!(None, document["meta"].get("genome_annotations"));

    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_rejects_node_without_divergence_or_date() -> Result<(), Report> {
    let mut topology = topology()?;
    set_branch_length(&mut topology, "A", None)?;
    let root_sequence = Seq::try_from_str("ACGT")?;
    let edge_mutations = no_mutations(&topology);
    assert_error!(
      auspice_tree_of(
        &optimize_graph(&topology, &root_sequence, &edge_mutations),
        CommandKind::Optimize
      ),
      "Auspice v2 node 'A' requires divergence or numerical date data"
    );
    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_divergence_sums_from_the_root() -> Result<(), Report> {
    let topology = topology_from("((A:0.125,B:0.25)AB:0.5,C:0.25)root;")?;

    let auspice = auspice_tree_of(&annotated(&topology), CommandKind::Ancestral)?;

    let ab = auspice_child(&auspice, "AB");
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
    let setup = ancestral_setup(Mutations::IndelAndAminoAcid)?;
    assert_error!(
      auspice_tree_of(&ancestral_graph(&setup), CommandKind::Ancestral),
      "Auspice v2 cannot represent amino-acid mutation track 'S/1:weird'"
    );
    Ok(())
  }

  #[test]
  fn test_tree_output_auspice_marks_excluded_nodes() -> Result<(), Report> {
    let mut setup = dated_setup()?;
    setup.excluded = btreeset! { node_key(&setup.topology, "B") };

    let auspice = auspice_tree_of(&dated_graph(&setup, None), CommandKind::Clock)?;

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
      Mutation::substitution(MutationTrack::Nucleotide, Sub::new(c(b'A'), 0_usize, c(b'T'))?),
      deletion(MutationTrack::Nucleotide, (1, 3), "CG")?,
      Mutation::substitution(
        MutationTrack::AminoAcid(o!("GENE")),
        Sub::new(c(b'K'), 4_usize, c(b'R'))?,
      ),
      deletion(MutationTrack::AminoAcid(o!("GENE")), (1, 3), "CG")?,
    ];

    let grouped = group_mutations(&mutations)?;

    let expected = btreemap! {
      o!("GENE") => vec![o!("K5R"), o!("C2-"), o!("G3-")],
      o!("nuc") => vec![o!("A1T")],
    };
    assert_eq!(expected, grouped);
    Ok(())
  }
}
