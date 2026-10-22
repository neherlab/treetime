#[cfg(test)]
pub(crate) mod tests {
  use crate::ambiguity_outputs::AmbiguityOutputs;
  use crate::annotated_graph::{
    AnnotatedGraph, AnnotatedTreeView, Divergence, TreeAminoAcids, TreeDates, TreeSequences, TreeTraits,
  };
  use crate::auspice::auspice_tree;
  use crate::command_outputs::CommandKind;
  use crate::output_plan::OutputSelection;
  use crate::usher_mat::mat_tree;
  use eyre::Report;
  use jsonschema::{Retrieve, Uri, Validator};
  use maplit::btreemap;
  use ndarray::{Array1, array};
  use serde_json::Value;
  use std::collections::{BTreeMap, BTreeSet};
  use std::error::Error;
  use std::io;
  use std::sync::LazyLock;
  use treetime::alphabet::alphabet::Alphabet;
  use treetime::alphabet::state_alphabet::StateAlphabet;
  use treetime::ancestral::aa::AaNodeData;
  use treetime::partition::marginal::discrete::state_posterior::StatePosterior;
  use treetime::partition::marginal::posterior_states::AmbiguityThreshold;
  use treetime::seq::mutation::{AlignedMutation, Mutation, MutationEvent, MutationTrack, PositionBase, Sub};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_graph::node_names::NodeNames;
  use treetime_io::auspice_types::{AuspiceTree, AuspiceTreeNode};
  use treetime_io::nwk::nwk_read;
  use treetime_io::usher_mat::UsherTree;
  use treetime_primitives::{AsciiChar, Seq};
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::io::json::{json_value_read_str, to_json_value};
  use treetime_utils::{make_report, o};
  use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

  pub(crate) static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

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
    pub(crate) names: NodeNames,
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
    pub(crate) states: StateAlphabet,
    pub(crate) posteriors: BTreeMap<GraphNodeKey, StatePosterior>,
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
      names: parsed.names.clone(),
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

  pub(crate) fn auspice_tree_of(graph: &AnnotatedGraph<'_>, command: CommandKind) -> Result<AuspiceTree, Report> {
    auspice_tree(&AnnotatedTreeView::new(graph)?, command, UPDATED)
  }

  pub(crate) fn mat_tree_of(graph: &AnnotatedGraph<'_>) -> Result<UsherTree, Report> {
    Ok(mat_tree(&AnnotatedTreeView::new(graph)?)?.tree)
  }

  pub(crate) fn no_mutations(topology: &Topology) -> BTreeMap<GraphEdgeKey, Vec<Mutation>> {
    topology.graph.get_edges().map(|edge| (edge.key(), vec![])).collect()
  }

  pub(crate) fn mutation_on(topology: &Topology, name: &str) -> Result<BTreeMap<GraphEdgeKey, Vec<Mutation>>, Report> {
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
        alphabet: &NUC_ALPHABET,
        root_sequence: &setup.root_sequence,
        edge_mutations: &setup.edge_mutations,
        amino_acids: setup.aa_node_data.as_ref().map(|node_data| TreeAminoAcids {
          node_data,
          cdses: &setup.aa_annotations,
        }),
        comment_positions: PositionBase::One,
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
        alphabet: &NUC_ALPHABET,
        root_sequence,
        edge_mutations,
        amino_acids: None,
        comment_positions: PositionBase::One,
      }),
      ..annotated(topology)
    }
  }

  pub(crate) fn mugration_setup() -> Result<MugrationSetup, Report> {
    let topology = topology()?;
    let values = alternating_trait_values(&topology);
    mugration_setup_from(topology, &["CH", "US"], &values)
  }

  pub(crate) fn alternating_trait_values(topology: &Topology) -> BTreeMap<GraphNodeKey, &'static str> {
    topology
      .graph
      .get_nodes()
      .enumerate()
      .map(|(index, node)| (node.key(), if index % 2 == 0 { "CH" } else { "US" }))
      .collect()
  }

  pub(crate) fn mugration_setup_from(
    topology: Topology,
    state_names: &[&str],
    values: &BTreeMap<GraphNodeKey, &str>,
  ) -> Result<MugrationSetup, Report> {
    let states = StateAlphabet::new(state_names.iter().copied())?;
    let posteriors = values
      .iter()
      .map(|(&key, value)| {
        let index = states
          .index_of(value)
          .ok_or_else(|| make_report!("'{value}' is not one of the trait states"))?;
        let profile = Array1::from_shape_fn(states.len(), |state| if state == index { 1.0 } else { 0.0 });
        Ok((key, StatePosterior::new(profile)?))
      })
      .collect::<Result<_, Report>>()?;
    Ok(MugrationSetup {
      topology,
      states,
      posteriors,
    })
  }

  pub(crate) fn mugration_graph<'a>(setup: &'a MugrationSetup, attribute: &'a str) -> AnnotatedGraph<'a> {
    AnnotatedGraph {
      traits: Some(TreeTraits {
        attribute,
        states: &setup.states,
        posteriors: &setup.posteriors,
        ambiguity: None,
      }),
      ..annotated(&setup.topology)
    }
  }

  pub(crate) fn mugration_graph_with_ambiguity<'a>(
    setup: &'a MugrationSetup,
    attribute: &'a str,
    ambiguity: &'a AmbiguityOutputs,
  ) -> AnnotatedGraph<'a> {
    let graph = mugration_graph(setup, attribute);
    AnnotatedGraph {
      traits: graph.traits.map(|traits| TreeTraits {
        ambiguity: Some(ambiguity),
        ..traits
      }),
      ..graph
    }
  }

  pub(crate) fn ambiguity_for(outputs: &[OutputSelection], threshold: f64) -> Result<AmbiguityOutputs, Report> {
    Ok(AmbiguityOutputs::new(
      AmbiguityThreshold::new(threshold)?,
      &BTreeSet::new(),
      outputs.iter().copied(),
    ))
  }

  pub(crate) fn ambiguous_mugration_setup() -> Result<MugrationSetup, Report> {
    let mut setup = mugration_setup()?;
    let root = node_key(&setup.topology, "root");
    setup.posteriors.insert(root, StatePosterior::new(array![0.45, 0.55])?);
    Ok(setup)
  }

  pub(crate) fn dated_setup() -> Result<DatedSetup, Report> {
    let topology = topology()?;
    let div = topology
      .graph
      .get_nodes()
      .enumerate()
      .map(|(index, node)| (node.key(), index.to_f64_exact() / 2.0))
      .collect();
    let num_date = topology
      .graph
      .get_nodes()
      .enumerate()
      .map(|(index, node)| (node.key(), Some(2020.0 + index.to_f64_exact())))
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
        alphabet: &NUC_ALPHABET,
        root_sequence,
        edge_mutations: &setup.edge_mutations,
        amino_acids: None,
        comment_positions: PositionBase::One,
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
    let ancestral = auspice_tree_of(&ancestral_graph(&ancestral), CommandKind::Ancestral)?;

    let optimize = topology()?;
    let optimize_reference = Seq::try_from_str("ACGT")?;
    let optimize_mutations = no_mutations(&optimize);
    let optimize = auspice_tree_of(
      &optimize_graph(&optimize, &optimize_reference, &optimize_mutations),
      CommandKind::Optimize,
    )?;

    let prune = topology()?;
    let prune = auspice_tree_of(&annotated(&prune), CommandKind::Prune)?;

    let clock = dated_setup()?;
    let clock = auspice_tree_of(&dated_graph(&clock, None), CommandKind::Clock)?;

    let mugration = mugration_setup()?;
    let mugration = auspice_tree_of(&mugration_graph(&mugration, "country"), CommandKind::Mugration)?;

    let ambiguous = ambiguous_mugration_setup()?;
    let ambiguity = ambiguity_for(&[OutputSelection::Auspice], 0.125)?;
    let ambiguous = auspice_tree_of(
      &mugration_graph_with_ambiguity(&ambiguous, "country", &ambiguity),
      CommandKind::Mugration,
    )?;

    let timetree = dated_setup()?;
    let timetree = auspice_tree_of(&dated_graph(&timetree, None), CommandKind::Timetree)?;

    Ok(btreemap! {
      "ancestral" => to_json_value(&ancestral)?,
      "optimize" => to_json_value(&optimize)?,
      "prune" => to_json_value(&prune)?,
      "clock" => to_json_value(&clock)?,
      "mugration" => to_json_value(&mugration)?,
      "mugration-ambiguity" => to_json_value(&ambiguous)?,
      "timetree" => to_json_value(&timetree)?,
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
      mat_tree_of(&AnnotatedGraph {
        sequences: root_sequence.map(|root_sequence| TreeSequences {
          alphabet: &NUC_ALPHABET,
          root_sequence,
          edge_mutations: &no_mutations,
          amino_acids: None,
          comment_positions: PositionBase::One,
        }),
        ..annotated(&topology)
      })
    })
    .collect()
  }

  pub(crate) fn auspice_validator() -> Result<Validator, Report> {
    let schema = json_value_read_str(AUSPICE_SCHEMA)?;
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
      .graph
      .get_nodes()
      .map(|node| node.key())
      .find(|key| &topology.names[*key] == name)
      .expect("fixture node must exist")
  }

  pub(crate) fn rename_node(topology: &mut Topology, key: GraphNodeKey, name: &str) -> Result<(), Report> {
    let input = topology
      .graph
      .get_nodes()
      .map(|node| node.key())
      .map(|other| {
        let other_name = if other == key { name } else { &topology.names[other] };
        (other, Some(other_name.to_owned()))
      })
      .collect();
    topology.names = NodeNames::from_input(input, &topology.graph)?.names;
    Ok(())
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

  struct AuspiceSchemaRetriever {
    schemas: BTreeMap<String, Value>,
  }

  impl AuspiceSchemaRetriever {
    fn new() -> Result<Self, Report> {
      Ok(Self {
        schemas: [AUSPICE_CONFIG_SCHEMA, ANNOTATIONS_SCHEMA, ROOT_SEQUENCE_SCHEMA]
          .into_iter()
          .map(|schema| {
            let schema = json_value_read_str(schema)?;
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
