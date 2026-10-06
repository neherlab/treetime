use crate::results::year_date::YearDate;
use eyre::Report;
use itertools::{Itertools, izip};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_with::skip_serializing_none;
use std::collections::BTreeSet;
use treetime::timetree::confidence::CI_FRACTION;
use treetime_io::auspice_types::{AuspiceColoring, AuspiceTree, AuspiceTreeNode};
use treetime_utils::datetime::year_fraction::year_fraction_days_between;

const CATEGORICAL: &str = "categorical";

const BAD_BRANCH: &str = "bad_branch";

const GENOTYPE: &str = "gt";

const MUTED_PALETTE: [&str; 9] = [
  "#332288", "#88ccee", "#44aa99", "#117733", "#999933", "#ddcc77", "#cc6677", "#882255", "#aa4499",
];

const BAD_BRANCH_YES: &str = "Yes";

const NUCLEOTIDE_MUTATIONS: &str = "nuc";

pub fn preorder(root: &AuspiceTreeNode) -> Vec<(&AuspiceTreeNode, Option<usize>)> {
  let mut nodes = vec![];
  let mut stack = vec![(root, None)];
  while let Some((node, parent)) = stack.pop() {
    let index = nodes.len();
    nodes.push((node, parent));
    stack.extend(node.children.iter().rev().map(|child| (child, Some(index))));
  }
  nodes
}

/// A tree a run wrote, read from its Auspice file.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
#[deser(skip_serializing_optionals)]
pub struct ResultTree {
  /// Nodes in preorder; the root comes first and every parent precedes its children.
  pub nodes: Vec<ResultNode>,
  /// Colorings the Auspice file offers.
  pub colorings: Vec<ResultColoring>,
  /// Coloring the Auspice file selects by default.
  pub default_color_by: Option<String>,
}

impl ResultTree {
  pub fn from_auspice(tree: &AuspiceTree) -> Result<Self, Report> {
    let order = preorder(&tree.tree);
    let mut nodes = order
      .iter()
      .map(|(node, parent)| result_node(node, *parent))
      .collect::<Result<Vec<_>, Report>>()?;
    for index in 1..nodes.len() {
      if let Some(parent) = nodes[index].parent {
        nodes[parent].children.push(index);
      }
    }
    for index in (0..nodes.len()).rev() {
      let tips = if nodes[index].children.is_empty() {
        1
      } else {
        nodes[index].children.iter().map(|&child| nodes[child].tips).sum()
      };
      nodes[index].tips = tips;
    }
    let colorings = result_colorings(tree)?;
    Ok(Self {
      nodes,
      colorings,
      default_color_by: tree.data.meta.display_defaults.color_by.clone(),
    })
  }

  pub fn root(&self) -> &ResultNode {
    &self.nodes[0]
  }

  pub fn tips(&self) -> impl Iterator<Item = &ResultNode> {
    self.nodes.iter().filter(|node| node.is_tip())
  }

  pub fn find(&self, name: &str) -> Option<usize> {
    self.nodes.iter().position(|node| node.name == name)
  }
}

/// One node of a result tree.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
#[deser(skip_serializing_optionals)]
pub struct ResultNode {
  /// Name of the node.
  pub name: String,
  /// Index of the parent node; absent for the root.
  pub parent: Option<usize>,
  /// Indices of the child nodes.
  pub children: Vec<usize>,
  /// Number of samples below the node, 1 for a sample.
  pub tips: usize,
  /// Divergence from the root.
  pub div: Option<f64>,
  /// Date of the node.
  pub date: Option<YearDate>,
  /// Confidence interval of the date; absent when the run computed none or the interval is empty.
  pub date_interval: Option<DateInterval>,
  /// Whether the clock model left the sample out, because it had no usable date or was a clock outlier.
  pub excluded: Option<bool>,
  /// Nucleotide mutations on the branch above the node.
  pub mutations: Vec<String>,
}

impl ResultNode {
  pub fn is_tip(&self) -> bool {
    self.children.is_empty()
  }
}

/// Confidence interval of a date.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
pub struct DateInterval {
  /// Earliest date.
  pub lower: YearDate,
  /// Latest date.
  pub upper: YearDate,
  /// Width of the interval in days.
  pub days: f64,
  /// Probability that the date lies inside the interval, for example `0.9`.
  pub level: f64,
}

impl DateInterval {
  pub fn new(lower: f64, upper: f64) -> Option<Self> {
    (upper > lower).then(|| Self {
      days: year_fraction_days_between(lower, upper),
      lower: YearDate::new(lower),
      upper: YearDate::new(upper),
      level: CI_FRACTION,
    })
  }
}

/// A coloring of the Auspice tree.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
pub struct ResultColoring {
  /// Node attribute the coloring reads.
  pub key: String,
  /// Title of the coloring.
  pub title: String,
  /// Kind of scale, for example `categorical` or `continuous`.
  pub kind: String,
  /// Distinct states of a categorical coloring, sorted; empty for other kinds.
  pub states: Vec<String>,
  /// Colors of the states in the displayed Auspice tree; empty when Auspice chooses them.
  pub scale: Vec<StateColor>,
}

/// Color of one state of a categorical coloring.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema, deser::Serialize, deser::Deserialize)]
pub struct StateColor {
  /// Value of the node attribute.
  pub state: String,
  /// CSS hex color.
  pub color: String,
}

pub fn result_colorings(tree: &AuspiceTree) -> Result<Vec<ResultColoring>, Report> {
  let order = preorder(&tree.tree);
  tree
    .data
    .meta
    .colorings
    .iter()
    .map(|coloring| {
      let states = if coloring.type_ == CATEGORICAL {
        order
          .iter()
          .map(|(node, _)| node.node_attrs.attr(&coloring.key))
          .filter_map_ok(|attr| attr.map(|attr| attr.value().to_owned()))
          .collect::<Result<BTreeSet<_>, Report>>()?
          .into_iter()
          .collect()
      } else {
        vec![]
      };
      let scale = coloring_scale(coloring, &states);
      Ok(ResultColoring {
        key: coloring.key.clone(),
        title: coloring.title.clone(),
        kind: coloring.type_.clone(),
        states,
        scale,
      })
    })
    .collect()
}

fn coloring_scale(coloring: &AuspiceColoring, states: &[String]) -> Vec<StateColor> {
  let muted = coloring.type_ == CATEGORICAL
    && coloring.key != GENOTYPE
    && coloring.key != BAD_BRANCH
    && (1..=MUTED_PALETTE.len()).contains(&states.len());
  if muted {
    izip!(states, MUTED_PALETTE)
      .map(|(state, color)| StateColor {
        state: state.clone(),
        color: color.to_owned(),
      })
      .collect()
  } else {
    coloring
      .scale
      .iter()
      .map(|[state, color]| StateColor {
        state: state.clone(),
        color: color.clone(),
      })
      .collect()
  }
}

fn result_node(node: &AuspiceTreeNode, parent: Option<usize>) -> Result<ResultNode, Report> {
  let attrs = &node.node_attrs;
  let num_date = attrs.num_date.as_ref();
  let excluded = attrs.attr(BAD_BRANCH)?.map(|attr| attr.value() == BAD_BRANCH_YES);
  Ok(ResultNode {
    name: node.name.clone(),
    parent,
    children: vec![],
    tips: 0,
    div: attrs.div,
    date: num_date.map(|date| YearDate::new(date.value)),
    date_interval: num_date
      .and_then(|date| date.confidence)
      .and_then(|[lower, upper]| DateInterval::new(lower, upper)),
    excluded,
    mutations: node
      .branch_attrs
      .mutations
      .get(NUCLEOTIDE_MUTATIONS)
      .cloned()
      .unwrap_or_default(),
  })
}
