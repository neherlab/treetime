use crate::commands::ancestral::aa_node_data::AaNodeData;
use crate::commands::ancestral::result::{AncestralNodeOut, AncestralOutputMaps};
use crate::commands::shared::tree_output::{
  NUC_TRACK, auspice_data, auspice_from_graph, cumulative_branch_length_from, generation_date, mat_from_graph,
  node_name_value, phyloxml_from_graph, sequence_auspice_node, sequence_phyloxml_clade, write_tree_outputs,
};
use crate::seq::mutation::{Mutation, MutationTrack};
use eyre::{Report, WrapErr};
use maplit::btreemap;
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::auspice_types::{
  AuspiceGenomeAnnotationCds, AuspiceGenomeAnnotationNuc, AuspiceGenomeAnnotations, AuspiceTree, Segments, StartEnd,
};
use treetime_io::graph::TreeWriteKind;
use treetime_io::nwk::CommentProviders;
use treetime_io::phyloxml::Phyloxml;
use treetime_io::usher_mat::UsherTree;
use treetime_utils::{make_error, make_report};
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

pub fn write_ancestral_tree_outputs(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, AncestralNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &AncestralOutputMaps,
  aa_node_data: Option<&AaNodeData>,
  outputs: &BTreeMap<TreeWriteKind, PathBuf>,
  providers: &CommentProviders,
) -> Result<(), Report> {
  let updated = generation_date();
  let names: BTreeMap<GraphNodeKey, Option<String>> =
    nodes.iter().map(|(key, node)| (*key, node.name.clone())).collect();
  write_tree_outputs(
    graph,
    &names,
    branch_lengths,
    branch_lengths,
    outputs,
    providers,
    "ancestral",
    || ancestral_to_auspice(graph, nodes, branch_lengths, maps, aa_node_data, &updated),
    || ancestral_to_phyloxml(graph, nodes, branch_lengths, maps, aa_node_data),
    || ancestral_to_mat(graph, &names, branch_lengths, maps, aa_node_data),
  )
}

pub(crate) fn ancestral_to_auspice(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, AncestralNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &AncestralOutputMaps,
  aa_node_data: Option<&AaNodeData>,
  updated: &str,
) -> Result<AuspiceTree, Report> {
  let root_sequences = ancestral_root_sequences(maps, aa_node_data);
  let genome_annotations = ancestral_genome_annotations(&root_sequences, aa_node_data)?;
  let data = auspice_data(
    "TreeTime ancestral analysis",
    updated,
    vec![],
    vec![],
    None,
    genome_annotations,
    Some(root_sequences),
    !ancestral_all_mutations(graph, maps, aa_node_data).is_empty(),
  );
  auspice_from_graph(graph, data, |context| {
    let out = &nodes[&context.node_key];
    let name = node_name_value(context.node_key, out.name.as_deref());
    let div = cumulative_branch_length_from(graph, branch_lengths, context.node_key)?;
    let mutations = ancestral_node_mutations(graph, maps, context.node_key, context.edge_key, aa_node_data);
    sequence_auspice_node(&name, div, out.confidence, mutations, None, None)
  })
}

pub(crate) fn ancestral_to_phyloxml(
  graph: &Graph,
  nodes: &BTreeMap<GraphNodeKey, AncestralNodeOut>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &AncestralOutputMaps,
  aa_node_data: Option<&AaNodeData>,
) -> Result<Phyloxml, Report> {
  phyloxml_from_graph(graph, "TreeTime ancestral analysis", |context| {
    let out = &nodes[&context.node_key];
    let name = node_name_value(context.node_key, out.name.as_deref());
    let div = cumulative_branch_length_from(graph, branch_lengths, context.node_key)?;
    let branch_length = context.edge_key.and_then(|edge_key| branch_lengths[&edge_key]);
    let mutations = ancestral_node_mutations(graph, maps, context.node_key, context.edge_key, aa_node_data);
    sequence_phyloxml_clade(
      &name,
      out.name.clone(),
      div,
      branch_length,
      out.confidence,
      mutations,
      &ancestral_node_sequences(graph, maps, context.node_key, aa_node_data),
    )
  })
}

pub(crate) fn ancestral_to_mat(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  nwk_weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
  maps: &AncestralOutputMaps,
  aa_node_data: Option<&AaNodeData>,
) -> Result<UsherTree, Report> {
  let reference = ancestral_root_sequences(maps, aa_node_data).remove(NUC_TRACK);
  mat_from_graph(graph, names, nwk_weights, reference.as_deref(), |node_key, edge_key| {
    Ok(ancestral_node_mutations(
      graph,
      maps,
      node_key,
      Some(edge_key),
      aa_node_data,
    ))
  })
}

fn ancestral_root_sequences(maps: &AncestralOutputMaps, aa_node_data: Option<&AaNodeData>) -> BTreeMap<String, String> {
  let mut sequences = BTreeMap::new();
  if let Some(sequence) = maps.root_sequence.as_ref() {
    sequences.insert(NUC_TRACK.to_owned(), sequence.to_string());
  }
  if let Some(aa) = aa_node_data {
    sequences.extend(aa.root_aa_sequences.clone());
  }
  sequences
}

fn ancestral_node_sequences(
  graph: &Graph,
  maps: &AncestralOutputMaps,
  node_key: GraphNodeKey,
  aa_node_data: Option<&AaNodeData>,
) -> BTreeMap<String, String> {
  let mut sequences = maps
    .node_sequences
    .get(&node_key)
    .map(|sequence| btreemap! { NUC_TRACK.to_owned() => sequence.to_string() })
    .unwrap_or_default();
  if graph.is_root(node_key)
    && let Some(aa) = aa_node_data
  {
    sequences.extend(aa.root_aa_sequences.clone());
  }
  sequences
}

fn ancestral_genome_annotations(
  root_sequences: &BTreeMap<String, String>,
  aa_node_data: Option<&AaNodeData>,
) -> Result<Option<AuspiceGenomeAnnotations>, Report> {
  let nuc = root_sequences
    .get(NUC_TRACK)
    .map(|sequence| -> Result<_, Report> {
      Ok(AuspiceGenomeAnnotationNuc {
        start: 1,
        end: isize::try_from(sequence.len()).wrap_err("Nucleotide sequence length does not fit Auspice coordinates")?,
        strand: Some("+".to_owned()),
        r#type: Some("source".to_owned()),
        other: Value::default(),
      })
    })
    .transpose()?;
  let cdses = aa_node_data
    .map(|data| {
      data
        .annotations
        .iter()
        .map(|(name, annotation)| Ok((name.clone(), auspice_cds_annotation(name, annotation)?)))
        .collect::<Result<BTreeMap<_, _>, Report>>()
    })
    .transpose()?
    .unwrap_or_default();
  if nuc.is_none() && cdses.is_empty() {
    return Ok(None);
  }
  Ok(Some(AuspiceGenomeAnnotations {
    nuc,
    cdses,
    other: Value::default(),
  }))
}

fn auspice_cds_annotation(
  name: &str,
  annotation: &AugurNodeDataJsonAnnotationEntry,
) -> Result<AuspiceGenomeAnnotationCds, Report> {
  let strand = annotation
    .strand
    .clone()
    .ok_or_else(|| make_report!("CDS annotation '{name}' has no strand for Auspice output"))?;
  let other = Value::Object(annotation.other.clone().into_iter().collect());
  let segments = if let Some(segments) = annotation.segments.as_ref() {
    if segments.is_empty() {
      return make_error!("CDS annotation '{name}' has no segments for Auspice output");
    }
    Segments::MultipleSegments {
      segments: segments
        .iter()
        .map(|segment| {
          Ok(StartEnd {
            start: auspice_coordinate(segment.start, name, "segment start")?,
            end: auspice_coordinate(segment.end, name, "segment end")?,
            other: Value::Object(segment.other.clone().into_iter().collect()),
          })
        })
        .collect::<Result<Vec<_>, Report>>()?,
      other,
    }
  } else {
    Segments::OneSegment(StartEnd {
      start: auspice_coordinate(
        annotation
          .start
          .ok_or_else(|| make_report!("CDS annotation '{name}' has no start for Auspice output"))?,
        name,
        "start",
      )?,
      end: auspice_coordinate(
        annotation
          .end
          .ok_or_else(|| make_report!("CDS annotation '{name}' has no end for Auspice output"))?,
        name,
        "end",
      )?,
      other,
    })
  };
  Ok(AuspiceGenomeAnnotationCds {
    r#type: annotation.entry_type.clone(),
    gene: None,
    color: None,
    display_name: None,
    description: None,
    strand: Some(strand),
    segments,
  })
}

fn auspice_coordinate(coordinate: i64, name: &str, field: &str) -> Result<isize, Report> {
  isize::try_from(coordinate)
    .wrap_err_with(|| format!("CDS annotation '{name}' {field} does not fit Auspice coordinates"))
}

fn ancestral_node_mutations(
  graph: &Graph,
  maps: &AncestralOutputMaps,
  node_key: GraphNodeKey,
  edge_key: Option<GraphEdgeKey>,
  aa_node_data: Option<&AaNodeData>,
) -> Vec<Mutation> {
  let mut mutations = edge_key
    .and_then(|edge_key| maps.edge_mutations.get(&edge_key).cloned())
    .unwrap_or_default();
  if let Some(aa) = aa_node_data
    && let Some(tracks) = aa.node_aa_mutations.get(&node_key)
  {
    mutations.extend(tracks.iter().flat_map(|(track, events)| {
      events.iter().cloned().map(|event| Mutation {
        track: MutationTrack::AminoAcid(track.clone()),
        event,
      })
    }));
  }
  mutations
}

fn ancestral_all_mutations(
  graph: &Graph,
  maps: &AncestralOutputMaps,
  aa_node_data: Option<&AaNodeData>,
) -> Vec<Mutation> {
  graph
    .get_nodes()
    .into_iter()
    .flat_map(|node| {
      let node = node.read_arc();
      ancestral_node_mutations(graph, maps, node.key(), node.inbound().first().copied(), aa_node_data)
    })
    .collect()
}
