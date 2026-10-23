use crate::alphabet::alphabet::Alphabet;
use crate::gtr::get_gtr::{GtrModel, GtrModelName};
use crate::partition::resolve::{PartitionLayout, PartitionPlan, ResolvedPartitions, resolve_partitions};
use crate::seq::alignment::{NodeSeqInput, SequenceSource};
use eyre::Report;
use itertools::{Itertools, izip};
use std::collections::BTreeMap;
use std::iter::repeat_with;
use std::path::PathBuf;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::fmt::string::counted;
use treetime_utils::{make_internal_error, make_internal_report};

#[derive(Debug)]
pub struct SequenceInput {
  pub alignments: Vec<SequenceSource>,
  pub partitions: Vec<PartitionPlan>,
}

impl SequenceInput {
  pub fn whole_alignments(alignments: Vec<SequenceSource>, model: &GtrModel) -> Result<Self, Report> {
    let lengths = alignments.iter().map(|alignment| alignment.length).collect_vec();
    let alphabet = alignments.first().map(|alignment| alignment.alphabet.name());
    let ResolvedPartitions { plans, .. } = resolve_partitions(&[], &lengths, model, alphabet)?;
    Ok(Self {
      alignments,
      partitions: plans,
    })
  }

  pub fn is_empty(&self) -> bool {
    self.alignments.is_empty()
  }

  pub fn total_length(&self) -> usize {
    self.alignments.iter().map(|alignment| alignment.length).sum()
  }

  pub(crate) fn has_inferred_model(&self) -> bool {
    self.partitions.iter().any(|plan| plan.model.is_infer())
  }

  pub(crate) fn split(self) -> Result<SplitInput, Report> {
    let (layouts, models): (Vec<_>, Vec<_>) = self
      .partitions
      .into_iter()
      .map(|PartitionPlan { layout, model }| (layout, model))
      .unzip();
    let SplitAlignments { alignments, partitions } = split_alignments(self.alignments, layouts)?;
    let partitions = izip!(partitions, models)
      .map(|(data, model)| SplitPartition { data, model })
      .collect();
    Ok(SplitInput { alignments, partitions })
  }
}

pub(crate) fn split_alignments(
  alignments: Vec<SequenceSource>,
  layouts: Vec<PartitionLayout>,
) -> Result<SplitAlignments, Report> {
  let n_partitions = layouts.len();
  let mut groups = alignments.iter().map(|_| vec![]).collect_vec();
  for (index, layout) in layouts.into_iter().enumerate() {
    let Some(group) = layout.alignment.checked_sub(1).and_then(|i| groups.get_mut(i)) else {
      return make_internal_error!(
        "{} belongs to alignment {}, but the run has {}",
        layout.label(),
        layout.alignment,
        counted(alignments.len(), "alignment")
      );
    };
    group.push((index, layout));
  }
  let mut partitions = repeat_with(|| None).take(n_partitions).collect_vec();
  let alignments = izip!(1.., alignments, groups)
    .map(|(number, source, group)| {
      let SequenceSource {
        path,
        alphabet,
        length,
        nodes,
        mask,
      } = source;
      for (index, data) in split_alignment(&alphabet, length, nodes, group)? {
        partitions[index] = Some(data);
      }
      Ok(AlignmentInfo {
        number,
        path,
        alphabet,
        length,
        mask,
      })
    })
    .collect::<Result<Vec<_>, Report>>()?;
  let partitions = partitions
    .into_iter()
    .map(|data| data.ok_or_else(|| make_internal_report!("A partition layout was not split from its alignment")))
    .try_collect()?;
  Ok(SplitAlignments { alignments, partitions })
}

pub(crate) struct SplitInput {
  pub(crate) alignments: Vec<AlignmentInfo>,
  pub(crate) partitions: Vec<SplitPartition>,
}

impl SplitInput {
  pub(crate) fn model_names(&self) -> Vec<GtrModelName> {
    self.partitions.iter().map(|partition| partition.model.name()).collect()
  }
}

pub(crate) struct SplitPartition {
  pub(crate) data: PartitionData,
  pub(crate) model: GtrModel,
}

pub(crate) struct SplitAlignments {
  pub(crate) alignments: Vec<AlignmentInfo>,
  pub(crate) partitions: Vec<PartitionData>,
}

#[derive(Clone, Debug)]
pub(crate) struct AlignmentInfo {
  pub(crate) number: usize,
  pub(crate) path: PathBuf,
  pub(crate) alphabet: Alphabet,
  pub(crate) length: usize,
  pub(crate) mask: Vec<bool>,
}

fn split_alignment(
  alphabet: &Alphabet,
  length: usize,
  nodes: BTreeMap<GraphNodeKey, NodeSeqInput>,
  mut group: Vec<(usize, PartitionLayout)>,
) -> Result<Vec<(usize, PartitionData)>, Report> {
  if let Some(layout) = group
    .iter()
    .map(|(_, layout)| layout)
    .find(|layout| layout.columns.iter().any(|&column| column >= length))
  {
    return make_internal_error!(
      "{} has a column outside its alignment of {}",
      layout.label(),
      counted(length, "site")
    );
  }
  if let [(_, layout)] = group.as_slice()
    && layout.covers(length)
  {
    let (index, layout) = group.remove(0);
    let alphabet = alphabet.clone();
    return Ok(vec![(
      index,
      PartitionData {
        layout,
        alphabet,
        nodes,
      },
    )]);
  }
  let mut parts = group
    .into_iter()
    .map(|(index, layout)| (index, layout, BTreeMap::new()))
    .collect_vec();
  for (key, NodeSeqInput { name, seq }) in nodes {
    if let Some(seq) = &seq
      && seq.len() != length
    {
      return make_internal_error!(
        "The sequence of '{name}' has {}, but its alignment has {}",
        counted(seq.len(), "site"),
        counted(length, "site")
      );
    }
    for (_, layout, part_nodes) in &mut parts {
      let seq = seq
        .as_ref()
        .map(|seq| layout.columns.iter().map(|&column| seq[column]).collect());
      part_nodes.insert(
        key,
        NodeSeqInput {
          name: name.clone(),
          seq,
        },
      );
    }
  }
  Ok(
    parts
      .into_iter()
      .map(|(index, layout, nodes)| {
        let alphabet = alphabet.clone();
        (
          index,
          PartitionData {
            layout,
            alphabet,
            nodes,
          },
        )
      })
      .collect(),
  )
}

pub(crate) struct PartitionData {
  pub(crate) layout: PartitionLayout,
  pub(crate) alphabet: Alphabet,
  pub(crate) nodes: BTreeMap<GraphNodeKey, NodeSeqInput>,
}
