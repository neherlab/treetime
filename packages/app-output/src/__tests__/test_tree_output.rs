#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, and scratch collections"
)]

#[cfg(test)]
pub(super) mod tests {
  use crate::ancestral_result::AncestralNodeOut;
  use crate::ancestral_tree_output::{ancestral_to_auspice, ancestral_to_mat, write_ancestral_tree_outputs};
  use crate::clock_tree_output::{clock_to_auspice, clock_to_mat};
  use crate::mugration_tree_output::{mugration_to_auspice, mugration_to_mat};
  use crate::optimize_tree_output::{optimize_to_auspice, optimize_to_mat};
  use crate::prune_tree_output::{prune_to_auspice, prune_to_mat};
  use crate::timetree_result::TimetreeOutputMaps;
  use crate::timetree_tree_output::{timetree_to_auspice, timetree_to_mat};
  use crate::tree_output::{TraitValue, auspice_node, format_number, group_mutations, with_branch_support};
  use approx::assert_ulps_eq;
  use eyre::{Report, WrapErr};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::Value;
  use serde_json::json;
  use std::collections::BTreeMap;
  use tempfile::TempDir;
  use treetime::ancestral::aa::AaNodeData;
  use treetime::progress::NoopProgress;
  use treetime::seq::mutation::{AlignedMutation, Mutation, MutationEvent, MutationTrack, Sub};
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::auspice_types::AuspiceGenomeAnnotationNuc;
  use treetime_io::graph::TreeWriteKind;
  use treetime_io::nwk::{CommentProviders, NwkStyle, nwk_read_str};
  use treetime_primitives::{AsciiChar, Seq};
  use treetime_utils::io::json::{JsonPretty, json_read_file, json_read_str, json_write_str};
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_tree_output_ancestral_models_preserve_semantics() -> Result<(), Report> {
    let (graph, names, branch_lengths, maps, aa_node_data, aa_annotations) =
      helpers::ancestral_graph(helpers::Mutations::NucleotideSubstitution)?;

    let nodes = helpers::ancestral_nodes(&names, &graph, &helpers::ancestral_confidences(&names, &graph));
    let auspice = ancestral_to_auspice(
      &graph,
      &nodes,
      &branch_lengths,
      &maps,
      aa_node_data.as_ref(),
      &aa_annotations,
      "2026-07-19",
    )?;
    let child = helpers::auspice_child(&auspice, "A");
    assert_eq!(Some("2026-07-19"), auspice.data.meta.updated.as_deref());
    assert_eq!(vec!["tree".to_owned(), "entropy".to_owned()], auspice.data.meta.panels);
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
    assert_eq!(vec!["A1T".to_owned()], child.branch_attrs.mutations["nuc"]);

    let mat = ancestral_to_mat(&graph, &names, &branch_lengths, &maps)?.tree;
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
    let (graph, names, branch_lengths, maps, aa_node_data, aa_annotations) =
      helpers::ancestral_graph(helpers::Mutations::Indel)?;
    let auspice = ancestral_to_auspice(
      &graph,
      &helpers::ancestral_nodes(&names, &graph, &btreemap! {}),
      &branch_lengths,
      &maps,
      aa_node_data.as_ref(),
      &aa_annotations,
      "2026-07-19",
    )?;
    let child = helpers::auspice_child(&auspice, "A");
    assert!(!child.branch_attrs.mutations.contains_key("nuc"));

    let (graph, names, branch_lengths, maps, aa_node_data, aa_annotations) =
      helpers::ancestral_graph(helpers::Mutations::AminoAcid)?;
    let auspice = ancestral_to_auspice(
      &graph,
      &helpers::ancestral_nodes(&names, &graph, &btreemap! {}),
      &branch_lengths,
      &maps,
      aa_node_data.as_ref(),
      &aa_annotations,
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
  fn test_tree_output_conversion_failure_does_not_create_target_and_keeps_prior_file() -> Result<(), Report> {
    let (graph, names, branch_lengths, mut maps, aa_node_data, aa_annotations) =
      helpers::ancestral_graph(helpers::Mutations::NucleotideSubstitution)?;
    maps.root_sequence = Seq::try_from_str("NCG")?;
    let dir = TempDir::new().wrap_err("When creating a temporary directory")?;
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
      &maps,
      aa_node_data.as_ref(),
      &aa_annotations,
      &outputs,
      &CommentProviders::new(),
      &NoopProgress,
    )
    .expect_err("MAT conversion must fail");
    assert!(format!("{error:?}").contains("root reference nucleotide 'N'"));
    assert!(nwk_path.is_file());
    assert!(!mat_path.exists());

    Ok(())
  }

  #[test]
  fn test_tree_output_graph_json_dumps_topology() -> Result<(), Report> {
    let (graph, names, branch_lengths, maps, aa_node_data, aa_annotations) =
      helpers::ancestral_graph(helpers::Mutations::None)?;
    let dir = TempDir::new().wrap_err("When creating a temporary directory")?;
    let path = dir.path().join("graph.json");
    let outputs = btreemap! { TreeWriteKind::GraphJson => path.clone() };

    write_ancestral_tree_outputs(
      &graph,
      &helpers::ancestral_nodes(&names, &graph, &btreemap! {}),
      &branch_lengths,
      &maps,
      aa_node_data.as_ref(),
      &aa_annotations,
      &outputs,
      &CommentProviders::new(),
      &NoopProgress,
    )?;
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
  fn test_tree_output_all_auspice_models_match_augur_v2_schema() -> Result<(), Report> {
    let documents = helpers::all_auspice_documents(helpers::InputBranchSupport::Absent)?;
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

  #[rstest]
  #[case::ancestral("ancestral")]
  #[case::optimize("optimize")]
  #[case::prune("prune")]
  #[case::mugration("mugration")]
  #[trace]
  fn test_tree_output_auspice_writes_input_branch_support(#[case] command: &str) -> Result<(), Report> {
    let document = &helpers::all_auspice_documents(helpers::InputBranchSupport::OnLeafA)?[command];
    let support: Vec<&Value> = document["tree"]["children"]
      .as_array()
      .expect("fixture root must have children")
      .iter()
      .filter(|child| child["name"] == "A")
      .map(|child| &child["node_attrs"]["confidence"])
      .collect();
    assert_eq!(vec![&json!({ "value": 0.9 })], support);
    Ok(())
  }

  #[test]
  fn test_tree_output_branch_support_rejects_a_trait_named_confidence() {
    let traits = btreemap! {
      o!("confidence") => TraitValue {
        value: o!("high"),
        confidence: BTreeMap::new(),
        entropy: None,
      },
    };
    let node = auspice_node(o!("A"), Some(0.1), None, None, None, traits, BTreeMap::new(), None);
    assert_error!(
      with_branch_support(node, Some(0.9), "mugration"),
      "Node 'A' has a trait named 'confidence', which Auspice JSON also uses for the input branch support. Rename the metadata column of the trait."
    );
  }

  #[test]
  fn test_tree_output_branch_support_rejects_non_finite_value_with_command_name() {
    let node = auspice_node(
      o!("A"),
      Some(0.1),
      None,
      None,
      None,
      BTreeMap::new(),
      BTreeMap::new(),
      None,
    );
    assert_error!(
      with_branch_support(node, Some(f64::NAN), "optimize"),
      "optimize node 'A' has non-finite input branch support=NaN"
    );
  }

  #[test]
  fn test_tree_output_auspice_root_sequence_adds_nuc_annotation_and_entropy_panel() -> Result<(), Report> {
    let (graph, names, _branch_lengths) = helpers::timetree_graph()?;
    let maps = TimetreeOutputMaps {
      root_sequence: Some(Seq::try_from_str("ACGT")?),
      ..TimetreeOutputMaps::default()
    };
    let auspice = timetree_to_auspice(
      &graph,
      &helpers::timetree_nodes(&names, &graph, &btreemap! {}),
      &maps,
      None,
      None,
      "2026-07-19",
    )?;

    assert_eq!(vec!["tree".to_owned(), "entropy".to_owned()], auspice.data.meta.panels);
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
        strand: Some("+".to_owned()),
        r#type: Some("source".to_owned()),
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
    let document = &helpers::all_auspice_documents(helpers::InputBranchSupport::Absent)?[command];

    assert_eq!(json!(["tree"]), document["meta"]["panels"]);
    assert_eq!(None, document["meta"].get("genome_annotations"));

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
    let (graph, names, branch_lengths, maps, aa_node_data, aa_annotations) =
      helpers::ancestral_graph(helpers::Mutations::IndelAndAminoAcid)?;
    let error = ancestral_to_auspice(
      &graph,
      &helpers::ancestral_nodes(&names, &graph, &btreemap! {}),
      &branch_lengths,
      &maps,
      aa_node_data.as_ref(),
      &aa_annotations,
      "2026-07-19",
    )
    .expect_err("Auspice must reject an amino-acid track outside its schema grammar");
    assert!(error.to_string().contains("cannot represent amino-acid mutation track"));
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
        MutationTrack::AminoAcid("GENE".to_owned()),
        Sub::new(helpers::c(b'K'), 4_usize, helpers::c(b'R'))?,
      ),
      helpers::deletion(MutationTrack::AminoAcid("GENE".to_owned()), (1, 3), "CG")?,
    ];

    let grouped = group_mutations(mutations)?;

    let expected = btreemap! {
      "GENE".to_owned() => vec!["K5R".to_owned(), "C2-".to_owned(), "G3-".to_owned()],
      "nuc".to_owned() => vec!["A1T".to_owned()],
    };
    assert_eq!(expected, grouped);
    Ok(())
  }

  pub(crate) mod helpers {
    use super::*;
    use crate::ancestral_result::AncestralOutputMaps;
    use crate::clock_result::ClockNodeOut;
    use crate::mugration_result::MugrationNodeOut;
    use crate::optimize_result::{OptimizeNodeOut, OptimizeOutputMaps};
    use crate::prune_result::{PruneNodeOut, PruneOutputMaps};
    use crate::timetree_result::{TimetreeNodeOut, TimetreeOutputMaps};
    use jsonschema::{Retrieve, Uri, Validator};
    use ndarray::array;
    use serde::Serialize;
    use std::collections::BTreeMap;
    use std::error::Error;
    use std::io;
    use treetime::gtr::gtr::GTR;
    use treetime::mugration::pipeline::MugrationOutput;
    use treetime::partition::storage::discrete::DiscreteStates;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_io::auspice_types::{AuspiceTree, AuspiceTreeNode};
    use treetime_io::usher_mat::UsherTree;
    use treetime_utils::make_report;
    use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

    const AUSPICE_SCHEMA: &str = include_str!("schemas/auspice/schema-export-v2.json");
    const AUSPICE_CONFIG_SCHEMA: &str = include_str!("schemas/auspice/schema-auspice-config-v2.json");
    const ANNOTATIONS_SCHEMA: &str = include_str!("schemas/auspice/schema-annotations.json");
    const ROOT_SEQUENCE_SCHEMA: &str = include_str!("schemas/auspice/schema-export-root-sequence.json");
    const MODEL_TREE: &str = "(A:0.1,B:0,C:0.5)root;";

    #[derive(Clone, Copy)]
    pub(crate) enum Mutations {
      None,
      NucleotideSubstitution,
      NucleotideSubstitutionAndAminoAcid,
      Indel,
      AminoAcid,
      IndelAndAminoAcid,
    }

    fn optimize_maps(_graph: &Graph) -> Result<OptimizeOutputMaps, Report> {
      Ok(OptimizeOutputMaps {
        root_sequence: Seq::try_from_str("ACGT")?,
        edge_mutations: BTreeMap::new(),
        edge_mutation_counts: BTreeMap::new(),
      })
    }

    fn prune_maps(_graph: &Graph) -> PruneOutputMaps {
      PruneOutputMaps::default()
    }

    fn timetree_maps(_graph: &Graph) -> TimetreeOutputMaps {
      TimetreeOutputMaps::default()
    }

    type AncestralGraphSetup = (
      Graph,
      BTreeMap<GraphNodeKey, Option<String>>,
      BTreeMap<GraphEdgeKey, Option<f64>>,
      AncestralOutputMaps,
      Option<AaNodeData>,
      BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
    );

    pub(crate) fn ancestral_graph(mutations: Mutations) -> Result<AncestralGraphSetup, Report> {
      let nwk_parsed = nwk_read_str("(A:0.5,B:0)root;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;
      let a_key = node_key(&graph, &names, "A");
      let b_key = node_key(&graph, &names, "B");
      let a_edge = graph.node_parent(a_key)?.unwrap().1;
      let b_edge = graph.node_parent(b_key)?.unwrap().1;
      let root_sequence = Seq::try_from_str("ACG")?;

      let include_substitution = matches!(
        mutations,
        Mutations::NucleotideSubstitution | Mutations::NucleotideSubstitutionAndAminoAcid
      );
      let include_indel = matches!(mutations, Mutations::Indel | Mutations::IndelAndAminoAcid);
      let include_aa = matches!(
        mutations,
        Mutations::AminoAcid | Mutations::NucleotideSubstitutionAndAminoAcid | Mutations::IndelAndAminoAcid
      );
      let mut a_mutations = Vec::new();
      if include_substitution {
        a_mutations.push(Mutation::substitution(
          MutationTrack::Nucleotide,
          Sub::new(c(b'A'), 0_usize, c(b'T'))?,
        ));
      }
      if include_indel {
        a_mutations.push(deletion(MutationTrack::Nucleotide, (1, 3), "CG")?);
      }
      let maps = AncestralOutputMaps {
        root_sequence,
        edge_mutations: btreemap! {
          a_edge => a_mutations,
          b_edge => vec![],
        },
      };

      let aa_node_data = include_aa.then(|| {
        let mut aa = AaNodeData::default();
        let track = if matches!(mutations, Mutations::IndelAndAminoAcid) {
          "S/1:weird"
        } else {
          "S"
        };
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
      let aa_annotations = if include_aa {
        btreemap! {
          "S".to_owned() => AugurNodeDataJsonAnnotationEntry {
            start: Some(1),
            end: Some(3),
            strand: Some("+".to_owned()),
            entry_type: Some("CDS".to_owned()),
            ..AugurNodeDataJsonAnnotationEntry::default()
          },
        }
      } else {
        btreemap! {}
      };
      Ok((graph, names, branch_lengths, maps, aa_node_data, aa_annotations))
    }

    pub(crate) fn ancestral_nodes(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    ) -> BTreeMap<GraphNodeKey, AncestralNodeOut> {
      graph
        .get_nodes()
        .map(|node| {
          let key = node.key();
          (
            key,
            AncestralNodeOut {
              name: names.get(&node.key()).cloned().flatten(),
              branch_support: confidences.get(&key).copied().flatten(),
            },
          )
        })
        .collect()
    }

    pub(crate) fn ancestral_confidences(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
    ) -> BTreeMap<GraphNodeKey, Option<f64>> {
      graph
        .get_nodes()
        .filter(|&node| names.get(&node.key()).and_then(|x| x.as_deref()) == Some("A"))
        .map(|node| (node.key(), Some(0.9)))
        .collect()
    }

    #[derive(Clone, Copy)]
    pub(crate) enum InputBranchSupport {
      Absent,
      OnLeafA,
    }

    fn branch_support(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      support: InputBranchSupport,
    ) -> BTreeMap<GraphNodeKey, Option<f64>> {
      match support {
        InputBranchSupport::Absent => BTreeMap::new(),
        InputBranchSupport::OnLeafA => ancestral_confidences(names, graph),
      }
    }

    pub(crate) fn ancestral_topology() -> Result<
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

    pub(crate) fn all_auspice_documents(support: InputBranchSupport) -> Result<BTreeMap<&'static str, Value>, Report> {
      let (ancestral_graph, ancestral_names, ancestral_bl, ancestral_maps, ancestral_aa, ancestral_aa_annotations) =
        ancestral_graph(Mutations::NucleotideSubstitution)?;
      let ancestral = ancestral_to_auspice(
        &ancestral_graph,
        &ancestral_nodes(
          &ancestral_names,
          &ancestral_graph,
          &branch_support(&ancestral_names, &ancestral_graph, support),
        ),
        &ancestral_bl,
        &ancestral_maps,
        ancestral_aa.as_ref(),
        &ancestral_aa_annotations,
        "2026-07-19",
      )?;
      let (optimize_graph, optimize_names, optimize_bl) = optimize_graph()?;
      let optimize = optimize_to_auspice(
        &optimize_graph,
        &optimize_nodes(
          &optimize_names,
          &optimize_graph,
          &branch_support(&optimize_names, &optimize_graph, support),
        ),
        &optimize_bl,
        &optimize_maps(&optimize_graph)?,
        "2026-07-19",
      )?;
      let (prune_graph, prune_names, prune_bl) = prune_graph()?;
      let prune = prune_to_auspice(
        &prune_graph,
        &prune_nodes(
          &prune_names,
          &prune_graph,
          &branch_support(&prune_names, &prune_graph, support),
        ),
        &prune_bl,
        &prune_maps(&prune_graph),
        "2026-07-19",
      )?;
      let (clock_graph, clock_names, _clock_bl) = clock_graph()?;
      let clock = clock_to_auspice(&clock_graph, &clock_nodes(&clock_names, &clock_graph), "2026-07-19")?;
      let (mugration_output, mugration_names, mugration_bl) = mugration_graph()?;
      let mugration = mugration_to_auspice(
        &mugration_output.graph,
        &mugration_nodes(
          &mugration_names,
          &mugration_output.graph,
          &branch_support(&mugration_names, &mugration_output.graph, support),
        ),
        &mugration_bl,
        &mugration_output,
        "country",
        "2026-07-19",
      )?;
      let (timetree_graph, timetree_names, _timetree_bl) = timetree_graph()?;
      let timetree = timetree_to_auspice(
        &timetree_graph,
        &timetree_nodes(
          &timetree_names,
          &timetree_graph,
          &branch_support(&timetree_names, &timetree_graph, support),
        ),
        &timetree_maps(&timetree_graph),
        None,
        None,
        "2026-07-19",
      )?;

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
      let (ancestral, ancestral_names, mut ancestral_bl) = ancestral_topology()?;
      set_mat_branch_lengths(&ancestral, &ancestral_names, &mut ancestral_bl)?;
      let (optimize, optimize_names, mut optimize_bl) = optimize_graph()?;
      set_mat_branch_lengths(&optimize, &optimize_names, &mut optimize_bl)?;
      let (prune, prune_names, mut prune_bl) = prune_graph()?;
      set_mat_branch_lengths(&prune, &prune_names, &mut prune_bl)?;
      let (clock, clock_names, mut clock_bl) = clock_graph()?;
      set_mat_branch_lengths(&clock, &clock_names, &mut clock_bl)?;
      let (mugration_output, mugration_names, mut mugration_bl) = mugration_graph()?;
      let mugration = &mugration_output.graph;
      set_mat_branch_lengths(mugration, &mugration_names, &mut mugration_bl)?;
      let (timetree, timetree_names, _timetree_bl) = timetree_graph()?;
      let timetree_weights = timetree_mat_nwk_weights(&timetree, &timetree_names)?;

      Ok(vec![
        ancestral_to_mat(
          &ancestral,
          &ancestral_names,
          &ancestral_bl,
          &AncestralOutputMaps {
            root_sequence: Seq::try_from_str("ACG")?,
            edge_mutations: BTreeMap::new(),
          },
        )?
        .tree,
        optimize_to_mat(&optimize, &optimize_names, &optimize_bl, &optimize_maps(&optimize)?)?.tree,
        prune_to_mat(&prune, &prune_names, &prune_bl, &prune_maps(&prune))?.tree,
        clock_to_mat(&clock, &clock_names, &clock_bl)?.tree,
        mugration_to_mat(mugration, &mugration_names, &mugration_bl)?.tree,
        timetree_to_mat(&timetree, &timetree_names, &timetree_weights, &timetree_maps(&timetree))?.tree,
      ])
    }

    pub(crate) fn auspice_validator() -> Result<Validator, Report> {
      let schema = json_read_str(AUSPICE_SCHEMA)?;
      Ok(
        jsonschema::draft6::options()
          .with_retriever(AuspiceSchemaRetriever::new()?)
          .build(&schema)?,
      )
    }

    pub(crate) fn optimize_auspice_without_required_node_data() -> Result<AuspiceTree, Report> {
      let (graph, names, mut branch_lengths) = optimize_graph()?;
      set_branch_length(&graph, &names, &mut branch_lengths, "A", None)?;
      optimize_to_auspice(
        &graph,
        &optimize_nodes(&names, &graph, &btreemap! {}),
        &branch_lengths,
        &optimize_maps(&graph)?,
        "2026-07-19",
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

    pub(crate) fn branch_length(
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
        .find_map(|node| (names.get(&node.key()).and_then(|x| x.as_deref()) == Some(name)).then(|| node.key()))
        .expect("fixture node must exist")
    }

    pub(crate) fn c(value: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(value)
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

    fn optimize_nodes(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    ) -> BTreeMap<GraphNodeKey, OptimizeNodeOut> {
      graph
        .get_nodes()
        .map(|node| {
          let key = node.key();
          (
            key,
            OptimizeNodeOut {
              name: names.get(&node.key()).cloned().flatten(),
              branch_support: confidences.get(&key).copied().flatten(),
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

    fn prune_nodes(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    ) -> BTreeMap<GraphNodeKey, PruneNodeOut> {
      graph
        .get_nodes()
        .map(|node| {
          let key = node.key();
          (
            key,
            PruneNodeOut {
              name: names.get(&node.key()).cloned().flatten(),
              branch_support: confidences.get(&key).copied().flatten(),
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

    fn clock_nodes(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
    ) -> BTreeMap<GraphNodeKey, ClockNodeOut> {
      graph
        .get_nodes()
        .enumerate()
        .map(|(index, node)| {
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

    fn mugration_nodes(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    ) -> BTreeMap<GraphNodeKey, MugrationNodeOut> {
      graph
        .get_nodes()
        .map(|node| {
          let key = node.key();
          (
            key,
            MugrationNodeOut {
              name: names.get(&node.key()).cloned().flatten(),
              branch_support: confidences.get(&key).copied().flatten(),
            },
          )
        })
        .collect()
    }

    fn mugration_graph() -> Result<
      (
        MugrationOutput,
        BTreeMap<GraphNodeKey, Option<String>>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str(MODEL_TREE)?;
      let names = nwk_parsed.names();
      let graph: Graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let states = DiscreteStates::from_values(["CH", "US"].into_iter(), "?");
      let gtr = GTR::builder().n_states(2).mu(1.0).pi(array![0.5, 0.5]).build()?;
      let (reconstructed_traits, confidences) = graph
        .get_nodes()
        .enumerate()
        .map(|(index, node)| {
          let key = node.key();
          let (state, profile) = if index % 2 == 0 {
            ("CH", array![1.0, 0.0])
          } else {
            ("US", array![0.0, 1.0])
          };
          ((key, Some(state.to_owned())), (key, Some(profile)))
        })
        .unzip();
      let output = MugrationOutput {
        n_states: 2,
        states,
        gtr,
        graph,
        reconstructed_traits,
        confidences,
      };
      Ok((output, names, branch_lengths))
    }

    pub(crate) fn timetree_graph() -> Result<
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

    pub(crate) fn timetree_nodes(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
    ) -> BTreeMap<GraphNodeKey, TimetreeNodeOut> {
      graph
        .get_nodes()
        .enumerate()
        .map(|(index, node)| {
          let key = node.key();
          (
            key,
            TimetreeNodeOut {
              name: names.get(&node.key()).cloned().flatten(),
              branch_support: confidences.get(&key).copied().flatten(),
              time: Some(2020.0 + index as f64),
              div: index as f64 / 2.0,
              is_outlier: false,
              bad_branch: false,
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
        .find_map(|node| (names.get(&node.key()).and_then(|x| x.as_deref()) == Some(name)).then(|| node.key()))
        .expect("fixture node must exist");
      let edge_key = graph.node_parent(key)?.expect("fixture node must have a parent").1;
      branch_lengths.insert(edge_key, length);
      Ok(())
    }

    fn timetree_mat_nwk_weights(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
    ) -> Result<BTreeMap<GraphEdgeKey, Option<f64>>, Report> {
      let mut weights: BTreeMap<GraphEdgeKey, Option<f64>> = graph.get_edges().map(|edge| (edge.key(), None)).collect();
      for (name, length) in [("A", None), ("B", Some(0.0)), ("C", Some(0.5))] {
        let key = graph
          .get_nodes()
          .into_iter()
          .find_map(|node| (names.get(&node.key()).and_then(|x| x.as_deref()) == Some(name)).then(|| node.key()))
          .expect("fixture node must exist");
        let edge_key = graph.node_parent(key)?.expect("fixture node must have a parent").1;
        weights.insert(edge_key, length);
      }
      Ok(weights)
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
