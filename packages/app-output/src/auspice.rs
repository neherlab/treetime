use crate::annotated_graph::{AnnotatedTreeView, Divergence, TreeSequences, TreeTraits};
use crate::output_plan::CommandKind;
use crate::trait_profile::{build_confidence_map, compute_entropy};
use eyre::{Report, WrapErr};
use itertools::izip;
use maplit::btreemap;
use serde_json::{Map, Value, json};
use std::collections::BTreeMap;
use treetime::clock::divergence::root_to_node_divergences_where_known;
use treetime::seq::mutation::{Mutation, MutationEvent, MutationTrack, mutation_event_strings};
use treetime_graph::assign_node_names::node_name_or_key;
use treetime_graph::node::GraphNodeKey;
use treetime_io::auspice_types::{
  AuspiceColoring, AuspiceDisplayDefaults, AuspiceGenomeAnnotationCds, AuspiceGenomeAnnotationNuc,
  AuspiceGenomeAnnotations, AuspiceNumDate, AuspiceTree, AuspiceTreeBranchAttrs, AuspiceTreeBranchAttrsLabels,
  AuspiceTreeData, AuspiceTreeMeta, AuspiceTreeNode, AuspiceTreeNodeAttr, AuspiceTreeNodeAttrs, Segments, StartEnd,
};
use treetime_utils::{make_error, make_internal_report, make_report};
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

const COLORING_BAD_BRANCH: &str = "bad_branch";
const COLORING_GENOTYPE: &str = "gt";
const COLORING_NUM_DATE: &str = "num_date";
const NUC_TRACK: &str = "nuc";

pub(crate) fn auspice_tree(
  tree: &AnnotatedTreeView<'_>,
  command: CommandKind,
  updated: &str,
) -> Result<AuspiceTree, Report> {
  let divergences = match tree.graph().divergence {
    Divergence::CumulativeBranchLength => {
      root_to_node_divergences_where_known(tree.tree(), tree.graph().divergence_branch_lengths)
    },
    Divergence::Values(values) => values.iter().map(|(&key, &div)| (key, Some(div))).collect(),
  };
  let mut nodes: BTreeMap<GraphNodeKey, AuspiceTreeNode> = BTreeMap::new();
  let mut has_mutations = false;
  for &key in tree.tree().preorder().iter().rev() {
    let mut node = auspice_node(tree, key, divergences[&key])?;
    has_mutations |= !node.branch_attrs.mutations.is_empty();
    node.children = tree
      .tree()
      .children(key)
      .iter()
      .map(|(child_key, _)| {
        nodes
          .remove(child_key)
          .ok_or_else(|| make_internal_report!("Auspice child node {child_key} was not converted"))
      })
      .collect::<Result<_, Report>>()?;
    nodes.insert(key, node);
  }
  let root_key = tree.tree().root();
  let root = nodes
    .remove(&root_key)
    .ok_or_else(|| make_internal_report!("Auspice root node {root_key} was not converted"))?;
  Ok(AuspiceTree {
    data: auspice_data(tree, command, updated, has_mutations)?,
    tree: root,
  })
}

pub(crate) fn ensure_finite(value: f64, node_name: &str, field: &str) -> Result<(), Report> {
  if value.is_finite() {
    Ok(())
  } else {
    make_error!("Node '{node_name}' has non-finite {field}={value}")
  }
}

pub(crate) fn group_mutations(mutations: &[Mutation]) -> Result<BTreeMap<String, Vec<String>>, Report> {
  let mut grouped: BTreeMap<String, Vec<String>> = BTreeMap::new();
  for mutation in mutations {
    let track = match &mutation.track {
      MutationTrack::Nucleotide if !matches!(mutation.event, MutationEvent::Substitution(_)) => continue,
      MutationTrack::Nucleotide => NUC_TRACK,
      MutationTrack::AminoAcid(track) => auspice_track_name(track)?,
    };
    grouped
      .entry(track.to_owned())
      .or_default()
      .extend(mutation_event_strings(&mutation.event)?);
  }
  Ok(grouped)
}

#[expect(
  clippy::as_conversions,
  clippy::expect_used,
  reason = "count/index numeric cast is exact for the domain range; expect on a value an upstream invariant guarantees is present"
)]
pub(crate) fn format_number(number: f64, precision: i32) -> f64 {
  if number == 0.0 || !number.is_finite() {
    return number;
  }
  let integral = number.abs().trunc();
  let significand = if integral >= 1.0 {
    integral.log10().floor() as i32 + 1
  } else {
    0
  };
  let significant_figures = (significand + precision).max(1) as usize;
  format!("{number:.*e}", significant_figures - 1)
    .parse()
    .expect("a float formatted in scientific notation must parse back")
}

fn auspice_data(
  tree: &AnnotatedTreeView<'_>,
  command: CommandKind,
  updated: &str,
  has_mutations: bool,
) -> Result<AuspiceTreeData, Report> {
  let graph = tree.graph();
  let mut colorings = vec![];
  let mut filters = vec![];
  if graph.dates.is_some() {
    colorings.push(coloring(COLORING_NUM_DATE, "Date", "continuous"));
    colorings.push(coloring(COLORING_BAD_BRANCH, "Excluded", "categorical"));
    filters.push(COLORING_BAD_BRANCH.to_owned());
  }
  if let Some(traits) = &graph.traits {
    colorings.push(coloring(traits.attribute, traits.attribute, "categorical"));
    filters.push(traits.attribute.to_owned());
  }
  if has_mutations {
    colorings.push(coloring(COLORING_GENOTYPE, "Genotype", "categorical"));
  }
  let color_by = filters.last().cloned();

  let root_sequences = graph.sequences.as_ref().map(root_sequences).unwrap_or_default();
  let cdses = graph
    .sequences
    .as_ref()
    .and_then(|sequences| sequences.amino_acids.as_ref())
    .map(|amino_acids| auspice_cds_annotations(amino_acids.cdses))
    .transpose()?
    .unwrap_or_default();
  let genome_annotations = genome_annotations(&root_sequences, cdses)?;
  let mut panels = vec!["tree".to_owned()];
  if genome_annotations.is_some() {
    panels.push("entropy".to_owned());
  }

  Ok(AuspiceTreeData {
    version: Some("v2".to_owned()),
    meta: AuspiceTreeMeta {
      title: Some(format!("TreeTime {} analysis", command.stem())),
      updated: Some(updated.to_owned()),
      panels,
      genome_annotations,
      colorings,
      filters,
      display_defaults: AuspiceDisplayDefaults {
        color_by,
        ..AuspiceDisplayDefaults::default()
      },
      ..AuspiceTreeMeta::default()
    },
    root_sequence: (!root_sequences.is_empty()).then_some(root_sequences),
    other: Value::default(),
  })
}

fn root_sequences(sequences: &TreeSequences<'_>) -> BTreeMap<String, String> {
  let mut root_sequences = btreemap! { NUC_TRACK.to_owned() => sequences.root_sequence.to_string() };
  if let Some(amino_acids) = &sequences.amino_acids {
    root_sequences.extend(amino_acids.node_data.root_aa_sequences.clone());
  }
  root_sequences
}

fn coloring(key: &str, title: &str, type_: &str) -> AuspiceColoring {
  AuspiceColoring {
    key: key.to_owned(),
    title: title.to_owned(),
    type_: type_.to_owned(),
    ..AuspiceColoring::default()
  }
}

fn genome_annotations(
  root_sequences: &BTreeMap<String, String>,
  cdses: BTreeMap<String, AuspiceGenomeAnnotationCds>,
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
  if nuc.is_none() && cdses.is_empty() {
    return Ok(None);
  }
  Ok(Some(AuspiceGenomeAnnotations {
    nuc,
    cdses,
    other: Value::default(),
  }))
}

fn auspice_node(tree: &AnnotatedTreeView<'_>, key: GraphNodeKey, div: Option<f64>) -> Result<AuspiceTreeNode, Report> {
  let graph = tree.graph();
  let name = node_name_or_key(key, graph.names[&key].as_deref());

  let div = finite_number(div, 6, &name, "div")?;
  let num_date = match &graph.dates {
    Some(dates) => finite_number(dates.num_date[&key], 3, &name, "date")?,
    None => None,
  };
  if div.is_none() && num_date.is_none() {
    return make_error!("Auspice v2 node '{name}' requires divergence or numerical date data");
  }
  let num_date_confidence = graph
    .dates
    .as_ref()
    .and_then(|dates| dates.confidence?.get(&key))
    .map(|&interval| date_confidence(interval, &name))
    .transpose()?;
  let (traits, labels) = match &graph.traits {
    Some(traits) => (
      node_trait(traits, key, &name)?,
      node_transition(tree, traits, key).map(|transition| AuspiceTreeBranchAttrsLabels {
        aa: None,
        clade: None,
        other: json!({ traits.attribute: transition }),
      }),
    ),
    None => (None, None),
  };
  let branch_support = graph
    .branch_support
    .and_then(|support| support.get(&key).copied().flatten());

  Ok(AuspiceTreeNode {
    branch_attrs: AuspiceTreeBranchAttrs {
      mutations: node_mutations(tree, key)?,
      labels,
      other: Value::default(),
    },
    node_attrs: AuspiceTreeNodeAttrs {
      div,
      num_date: num_date.map(|value| AuspiceNumDate {
        value,
        confidence: num_date_confidence,
      }),
      bad_branch: graph
        .dates
        .as_ref()
        .map(|dates| AuspiceTreeNodeAttr::new(if dates.excluded.contains(&key) { "Yes" } else { "No" })),
      clade_membership: None,
      region: None,
      country: None,
      division: None,
      other: node_attrs_other(traits, branch_support, &name)?,
    },
    name,
    children: vec![],
    other: Value::default(),
  })
}

fn node_mutations(tree: &AnnotatedTreeView<'_>, key: GraphNodeKey) -> Result<BTreeMap<String, Vec<String>>, Report> {
  let Some(sequences) = &tree.graph().sequences else {
    return Ok(BTreeMap::new());
  };
  let mut grouped = match tree.tree().parent(key) {
    Some((_, edge_key)) => group_mutations(&sequences.edge_mutations[&edge_key])?,
    None => BTreeMap::new(),
  };
  if let Some(tracks) = sequences
    .amino_acids
    .as_ref()
    .and_then(|amino_acids| amino_acids.node_data.node_aa_mutations.get(&key))
  {
    for (track, events) in tracks.iter().filter(|(_, events)| !events.is_empty()) {
      let strings = grouped.entry(auspice_track_name(track)?.to_owned()).or_default();
      for event in events {
        strings.extend(mutation_event_strings(event)?);
      }
    }
  }
  Ok(grouped)
}

fn auspice_track_name(track: &str) -> Result<&str, Report> {
  let valid = !track.is_empty()
    && track
      .bytes()
      .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'*' | b'_' | b'.' | b'(' | b')' | b'-'));
  if valid {
    Ok(track)
  } else {
    make_error!("Auspice v2 cannot represent amino-acid mutation track '{track}'")
  }
}

fn node_trait(traits: &TreeTraits<'_>, key: GraphNodeKey, node_name: &str) -> Result<Option<(String, Value)>, Report> {
  let Some(value) = &traits.values[&key] else {
    return Ok(None);
  };
  let mut fields = Map::new();
  fields.insert("value".to_owned(), json!(value));
  if let Some(profile) = &traits.profiles[&key] {
    for (state, probability) in izip!(traits.states.iter(), profile) {
      ensure_finite(*probability, node_name, &format!("trait state '{state}' probability"))?;
    }
    let confidence: BTreeMap<String, f64> = build_confidence_map(traits.states, profile)
      .into_iter()
      .map(|(state, probability)| (state, format_number(probability, 3)))
      .collect();
    let entropy = compute_entropy(profile);
    ensure_finite(entropy, node_name, "trait entropy")?;
    if !confidence.is_empty() {
      fields.insert("confidence".to_owned(), json!(confidence));
    }
    fields.insert("entropy".to_owned(), json!(format_number(entropy, 3)));
  }
  Ok(Some((traits.attribute.to_owned(), Value::Object(fields))))
}

fn node_transition(tree: &AnnotatedTreeView<'_>, traits: &TreeTraits<'_>, key: GraphNodeKey) -> Option<String> {
  let (parent_key, _) = tree.tree().parent(key)?;
  match (&traits.values[&parent_key], &traits.values[&key]) {
    (Some(parent), Some(child)) if parent != child => Some(format!("{parent} \u{2192} {child}")),
    _ => None,
  }
}

fn node_attrs_other(
  traits: Option<(String, Value)>,
  branch_support: Option<f64>,
  node_name: &str,
) -> Result<Value, Report> {
  let mut attrs = Map::new();
  if let Some((attribute, value)) = traits {
    attrs.insert(attribute, value);
  }
  if let Some(branch_support) = branch_support {
    ensure_finite(branch_support, node_name, "input branch support")?;
    if attrs.contains_key("confidence") {
      return make_error!(
        "Node '{node_name}' has a trait named 'confidence', which Auspice JSON also uses for the input branch support. \
         Rename the metadata column of the trait."
      );
    }
    attrs.insert("confidence".to_owned(), json!({ "value": branch_support }));
  }
  Ok(Value::Object(attrs))
}

fn finite_number(value: Option<f64>, precision: i32, node_name: &str, field: &str) -> Result<Option<f64>, Report> {
  value
    .map(|value| {
      ensure_finite(value, node_name, field)?;
      Ok(format_number(value, precision))
    })
    .transpose()
}

fn date_confidence([lower, upper]: [f64; 2], node_name: &str) -> Result<[f64; 2], Report> {
  ensure_finite(lower, node_name, "date confidence lower bound")?;
  ensure_finite(upper, node_name, "date confidence upper bound")?;
  Ok([format_number(lower, 3), format_number(upper, 3)])
}

fn auspice_cds_annotations(
  aa_annotations: &BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
) -> Result<BTreeMap<String, AuspiceGenomeAnnotationCds>, Report> {
  aa_annotations
    .iter()
    .map(|(name, annotation)| Ok((name.clone(), auspice_cds_annotation(name, annotation)?)))
    .collect()
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
    let start = annotation
      .start
      .ok_or_else(|| make_report!("CDS annotation '{name}' has no start for Auspice output"))?;
    let end = annotation
      .end
      .ok_or_else(|| make_report!("CDS annotation '{name}' has no end for Auspice output"))?;
    Segments::OneSegment(StartEnd {
      start: auspice_coordinate(start, name, "start")?,
      end: auspice_coordinate(end, name, "end")?,
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
