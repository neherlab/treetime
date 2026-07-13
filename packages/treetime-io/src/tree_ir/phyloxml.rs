//! PhyloXML adapter for the format-neutral TreeIR graph.
//!
//! PhyloXML carries name, branch length, and date natively. TreeTime-specific
//! data (divergence, discrete traits, mutations, indels, branch-rate multiplier,
//! bad-branch flag) has no native PhyloXML element, so it is projected into
//! generic `<property>` elements with descriptive `ref` attributes and the XSD
//! `datatype`/`applies_to` qualifiers. The `treetime:` ref prefix is a placeholder
//! pending a dedicated property namespace (see
//! `kb/proposals/phyloxml-treetime-property-namespace.md`).

use crate::phyloxml::{
  Phyloxml, PhyloxmlClade, PhyloxmlContext, PhyloxmlDataFromGraphData, PhyloxmlDataToGraphData, PhyloxmlDate,
  PhyloxmlFromGraph, PhyloxmlGraphContext, PhyloxmlPhylogeny, PhyloxmlProperty, PhyloxmlToGraph,
};
use crate::tree_ir::mutation::{IndelKind, TreeIrIndel, TreeIrSub};
use crate::tree_ir::types::{TreeIrData, TreeIrEdge, TreeIrNode, TreeIrTrait};
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use treetime_primitives::AsciiChar;
use treetime_utils::{make_error, make_report};

const BRANCH_LENGTH_UNIT: &str = "subs/site";

const REF_DIV: &str = "treetime:divergence";
const REF_BAD_BRANCH: &str = "treetime:bad_branch";
const REF_DATE_INFERRED: &str = "treetime:date_inferred";
const REF_GAMMA: &str = "treetime:gamma";
const REF_MUTATION: &str = "treetime:mutation";
const REF_INDEL: &str = "treetime:indel";
const REF_TRAIT_PREFIX: &str = "treetime:trait:";

const APPLIES_NODE: &str = "node";
const APPLIES_BRANCH: &str = "parent_branch";

const DT_DOUBLE: &str = "xsd:double";
const DT_STRING: &str = "xsd:string";
const DT_BOOLEAN: &str = "xsd:boolean";

impl PhyloxmlDataFromGraphData for TreeIrData {
  fn phyloxml_data_from_graph_data(&self) -> Result<Phyloxml, Report> {
    let phylogeny = PhyloxmlPhylogeny {
      rooted: true,
      rerootable: None,
      branch_length_unit: Some(BRANCH_LENGTH_UNIT.to_owned()),
      phylogeny_type: None,
      name: self.title.clone(),
      id: None,
      description: self.description.clone(),
      date: None,
      confidence: vec![],
      clade: None,
      clade_relation: vec![],
      sequence_relation: vec![],
      property: vec![],
      other: BTreeMap::new(),
    };
    Ok(Phyloxml {
      phylogeny: vec![phylogeny],
      other: BTreeMap::new(),
    })
  }
}

impl PhyloxmlDataToGraphData for TreeIrData {
  fn phyloxml_data_to_graph_data(pxml: &Phyloxml) -> Result<Self, Report> {
    let phylogeny = pxml
      .phylogeny
      .first()
      .ok_or_else(|| make_report!("PhyloXML has no phylogeny"))?;
    Ok(Self {
      title: phylogeny.name.clone(),
      description: phylogeny.description.clone(),
      ..Self::default()
    })
  }
}

impl PhyloxmlFromGraph<TreeIrNode, TreeIrEdge, TreeIrData> for () {
  fn phyloxml_node_from_graph_components(
    context: &PhyloxmlGraphContext<TreeIrNode, TreeIrEdge, TreeIrData>,
  ) -> Result<PhyloxmlClade, Report> {
    let node = context.node;
    let edge = context.edge;

    let mut property = vec![];

    if let Some(div) = node.div {
      property.push(make_property(REF_DIV, DT_DOUBLE, APPLIES_NODE, &div.to_string()));
    }
    if node.bad_branch {
      property.push(make_property(REF_BAD_BRANCH, DT_BOOLEAN, APPLIES_NODE, "true"));
    }
    if node.date_is_inferred {
      property.push(make_property(REF_DATE_INFERRED, DT_BOOLEAN, APPLIES_NODE, "true"));
    }
    for (attr, value) in &node.traits {
      let ref_ = format!("{REF_TRAIT_PREFIX}{attr}");
      property.push(make_property(&ref_, DT_STRING, APPLIES_NODE, &value.value));
    }

    if let Some(edge) = edge {
      if let Some(gamma) = edge.gamma {
        property.push(make_property(REF_GAMMA, DT_DOUBLE, APPLIES_BRANCH, &gamma.to_string()));
      }
      for sub in &edge.mutations {
        property.push(make_property(
          REF_MUTATION,
          DT_STRING,
          APPLIES_BRANCH,
          &encode_mutation(sub),
        ));
      }
      for indel in &edge.indels {
        property.push(make_property(REF_INDEL, DT_STRING, APPLIES_BRANCH, &encode_indel(indel)));
      }
    }

    let date = node.date.map(|value| {
      let (minimum, maximum) = match node.date_confidence {
        Some([lo, hi]) => (Some(lo), Some(hi)),
        None => (None, None),
      };
      PhyloxmlDate {
        desc: None,
        value: Some(value),
        minimum,
        maximum,
        unit: Some("year".to_owned()),
      }
    });

    Ok(PhyloxmlClade {
      name: node.name.clone(),
      branch_length_elem: edge.and_then(|edge| edge.branch_length),
      branch_length_attr: None,
      confidence: vec![],
      width: None,
      color: None,
      node_id: None,
      taxonomy: vec![],
      sequence: vec![],
      events: None,
      binary_characters: None,
      distribution: vec![],
      date,
      reference: vec![],
      property,
      clade: vec![],
      other: BTreeMap::new(),
    })
  }
}

impl PhyloxmlToGraph<TreeIrNode, TreeIrEdge, TreeIrData> for () {
  fn phyloxml_node_to_graph_components(context: &PhyloxmlContext) -> Result<(TreeIrNode, TreeIrEdge), Report> {
    let clade = context.clade;

    let mut node = TreeIrNode {
      name: clade.name.clone(),
      ..TreeIrNode::default()
    };
    let mut edge = TreeIrEdge {
      branch_length: clade.branch_length_elem.or(clade.branch_length_attr),
      ..TreeIrEdge::default()
    };

    if let Some(date) = &clade.date {
      node.date = date.value;
      node.date_confidence = match (date.minimum, date.maximum) {
        (Some(lo), Some(hi)) => Some([lo, hi]),
        _ => None,
      };
    }

    for prop in &clade.property {
      match prop.ref_.as_str() {
        REF_DIV => node.div = Some(parse_f64(&prop.value, REF_DIV)?),
        REF_BAD_BRANCH => node.bad_branch = prop.value == "true",
        REF_DATE_INFERRED => node.date_is_inferred = prop.value == "true",
        REF_GAMMA => edge.gamma = Some(parse_f64(&prop.value, REF_GAMMA)?),
        REF_MUTATION => edge.mutations.push(decode_mutation(&prop.value)?),
        REF_INDEL => edge.indels.push(decode_indel(&prop.value)?),
        other if other.starts_with(REF_TRAIT_PREFIX) => {
          let attr = other.trim_start_matches(REF_TRAIT_PREFIX).to_owned();
          node.traits.insert(
            attr,
            TreeIrTrait {
              value: prop.value.clone(),
              ..TreeIrTrait::default()
            },
          );
        },
        _ => {},
      }
    }

    Ok((node, edge))
  }
}

fn make_property(ref_: &str, datatype: &str, applies_to: &str, value: &str) -> PhyloxmlProperty {
  PhyloxmlProperty {
    value: value.to_owned(),
    ref_: ref_.to_owned(),
    unit: None,
    datatype: datatype.to_owned(),
    applies_to: applies_to.to_owned(),
    id_ref: None,
  }
}

/// Encode a substitution as `gene:A123T` so the gene track survives the round-trip.
fn encode_mutation(sub: &TreeIrSub) -> String {
  format!("{}:{}", sub.gene, sub.to_auspice_string())
}

fn decode_mutation(value: &str) -> Result<TreeIrSub, Report> {
  let (gene, mut_str) = value
    .split_once(':')
    .ok_or_else(|| make_report!("Invalid mutation property '{value}': expected 'gene:A123T'"))?;
  TreeIrSub::from_auspice_string(gene, mut_str)
}

/// Encode an indel as `gene:ins|del:start:SEQ`.
fn encode_indel(indel: &TreeIrIndel) -> String {
  let kind = match indel.kind {
    IndelKind::Insertion => "ins",
    IndelKind::Deletion => "del",
  };
  let seq: String = indel.seq.iter().map(|c| char::from(*c)).collect();
  format!("{}:{}:{}:{}", indel.gene, kind, indel.start, seq)
}

fn decode_indel(value: &str) -> Result<TreeIrIndel, Report> {
  let parts: Vec<&str> = value.splitn(4, ':').collect();
  if parts.len() != 4 {
    return make_error!("Invalid indel property '{value}': expected 'gene:ins|del:start:SEQ'");
  }
  let kind = match parts[1] {
    "ins" => IndelKind::Insertion,
    "del" => IndelKind::Deletion,
    other => return make_error!("Invalid indel kind '{other}' in '{value}'"),
  };
  let start: usize = parts[2]
    .parse()
    .wrap_err_with(|| format!("Invalid indel start in '{value}'"))?;
  let seq = parts[3]
    .bytes()
    .map(AsciiChar::try_new)
    .collect::<Result<Vec<_>, _>>()?;
  Ok(TreeIrIndel {
    gene: parts[0].to_owned(),
    kind,
    start,
    seq,
  })
}

fn parse_f64(value: &str, ref_: &str) -> Result<f64, Report> {
  value
    .parse()
    .wrap_err_with(|| format!("Invalid numeric value '{value}' for property '{ref_}'"))
}
