use crate::partition::resolve::{PartitionLayout, partition_label};
use crate::seq::indel::{InDel, InDelKind, SplitFrom, compose_indels, sort_indels};
use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
use crate::seq::sink::AmbiguousSites;
use crate::{make_internal_error, make_internal_report};
use eyre::Report;
use itertools::{Itertools, izip};
use std::collections::BTreeMap;
use treetime_primitives::{AsciiChar, Seq};

pub(crate) fn alignment_groups<'a>(
  partitions: impl IntoIterator<Item = (Option<&'a PartitionLayout>, usize)>,
) -> Result<Vec<AlignmentGroup>, Report> {
  izip!(0.., 1.., partitions)
    .into_group_map_by(|&(_, number, (layout, _))| layout.map_or(number, |layout| layout.alignment))
    .into_iter()
    .sorted_by_key(|&(alignment, _)| alignment)
    .map(|(alignment, members)| {
      let layouts: Option<Vec<&PartitionLayout>> = members.iter().map(|&(_, _, (layout, _))| layout).collect();
      let layout = match (layouts, members.as_slice()) {
        (Some(layouts), _) => {
          let length = layouts.iter().map(|layout| layout.columns.len()).sum();
          AlignmentLayout::new(length, layouts)?
        },
        (None, &[(_, number, (None, length))]) => {
          AlignmentLayout::whole(alignment, partition_label(number, None), length)
        },
        (None, _) => {
          return make_internal_error!(
            "Alignment {alignment} has a partition without a layout next to other partitions"
          );
        },
      };
      Ok(AlignmentGroup {
        layout,
        partitions: members.iter().map(|&(index, ..)| index).collect(),
      })
    })
    .collect()
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AlignmentGroup {
  pub(crate) layout: AlignmentLayout,
  pub(crate) partitions: Vec<usize>,
}

impl AlignmentGroup {
  pub(crate) fn alignment(&self) -> usize {
    self.layout.alignment
  }

  pub(crate) fn members<'a, T>(&self, items: &'a [T]) -> Vec<&'a T> {
    self.partitions.iter().map(|&index| &items[index]).collect()
  }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlignmentLayout {
  alignment: usize,
  columns: AlignmentColumns,
}

impl AlignmentLayout {
  pub fn new<'a>(length: usize, layouts: impl IntoIterator<Item = &'a PartitionLayout>) -> Result<Self, Report> {
    let layouts = layouts.into_iter().collect_vec();
    let Some(alignment) = layouts.first().map(|layout| layout.alignment) else {
      return make_internal_error!("An alignment layout needs at least one partition");
    };
    if let [layout] = layouts.as_slice()
      && layout.covers(length)
    {
      return Ok(Self::whole(alignment, layout.label(), length));
    }
    let mut owners: Vec<Option<Site>> = vec![None; length];
    for (partition, layout) in layouts.iter().enumerate() {
      if layout.alignment != alignment {
        return make_internal_error!(
          "{} belongs to alignment {}, but the layout is of alignment {alignment}",
          layout.label(),
          layout.alignment
        );
      }
      if let Some((a, b)) = layout.columns.iter().tuple_windows().find(|(a, b)| a >= b) {
        return make_internal_error!("Columns of {} are not ascending: {a} before {b}", layout.label());
      }
      for (local, &column) in layout.columns.iter().enumerate() {
        let owner = owners.get_mut(column).ok_or_else(|| {
          make_internal_report!(
            "{} takes column {column}, but alignment {alignment} has {length} columns",
            layout.label()
          )
        })?;
        if let Some(other) = owner {
          return make_internal_error!(
            "{} and {} both take column {column} of alignment {alignment}",
            layouts[other.partition].label(),
            layout.label()
          );
        }
        *owner = Some(Site { partition, local });
      }
    }
    let sites = izip!(0.., owners)
      .map(|(column, owner)| {
        owner.ok_or_else(|| make_internal_report!("No partition takes column {column} of alignment {alignment}"))
      })
      .try_collect()?;
    let partitions = layouts
      .iter()
      .map(|layout| PartitionColumns {
        label: layout.label(),
        columns: layout.columns.clone(),
      })
      .collect();
    Ok(Self {
      alignment,
      columns: AlignmentColumns::Split(SplitColumns { partitions, sites }),
    })
  }

  pub(crate) fn whole(alignment: usize, label: String, length: usize) -> Self {
    Self {
      alignment,
      columns: AlignmentColumns::Whole(WholeColumns { label, length }),
    }
  }

  pub fn is_whole(&self) -> bool {
    matches!(self.columns, AlignmentColumns::Whole(_))
  }

  pub fn assemble_sequence(&self, sequences: &[Seq]) -> Result<Seq, Report> {
    self.check_count(sequences.len(), "sequences")?;
    for (partition, seq) in sequences.iter().enumerate() {
      let (label, length) = self.partition_span(partition);
      if seq.len() != length {
        return make_internal_error!(
          "The sequence of {label} has length {}, but the partition has {length} columns",
          seq.len(),
        );
      }
    }
    Ok(match &self.columns {
      AlignmentColumns::Whole(_) => sequences[0].clone(),
      AlignmentColumns::Split(split) => split
        .sites
        .iter()
        .map(|site| sequences[site.partition][site.local])
        .collect(),
    })
  }

  pub fn assemble_ambiguity(&self, ambiguous_sites: &[AmbiguousSites]) -> Result<AmbiguousSites, Report> {
    self.check_count(ambiguous_sites.len(), "ambiguity reports")?;
    izip!(0.., ambiguous_sites)
      .flat_map(|(partition, sites)| {
        sites
          .iter()
          .map(move |(&local, &code)| Ok((self.column(partition, local)?, code)))
      })
      .collect()
  }

  pub fn assemble_substitutions(&self, substitutions: &[Vec<Sub>]) -> Result<Vec<Sub>, Report> {
    self.check_count(substitutions.len(), "substitution lists")?;
    let mut assembled: Vec<Sub> = izip!(0.., substitutions)
      .flat_map(|(partition, subs)| {
        subs
          .iter()
          .map(move |sub| Sub::new(sub.reff(), self.column(partition, sub.pos())?, sub.qry()))
      })
      .try_collect()?;
    assembled.sort();
    Ok(assembled)
  }

  pub fn assemble_indels(&self, indels: &[Vec<InDel>]) -> Result<Vec<InDel>, Report> {
    self.check_count(indels.len(), "indel lists")?;
    for (partition, pieces) in indels.iter().enumerate() {
      for piece in pieces {
        self.check_piece(partition, piece)?;
      }
    }
    let split = match &self.columns {
      AlignmentColumns::Whole(_) => return Ok(indels[0].clone()),
      AlignmentColumns::Split(split) => split,
    };
    let pieces = izip!(0.., indels)
      .flat_map(|(partition, pieces)| pieces.iter().map(move |piece| (partition, piece)))
      .collect_vec();
    match pieces
      .iter()
      .map(|(_, piece)| piece.split_from.is_some())
      .all_equal_value()
    {
      Ok(true) => split.join_split_pieces(self.alignment, &pieces),
      Ok(false) => Ok(join_runs(
        pieces
          .iter()
          .flat_map(|&(partition, piece)| split.partitions[partition].runs(partition, piece))
          .collect(),
      )),
      Err(None) => Ok(vec![]),
      Err(Some(_)) => make_internal_error!(
        "An edge of alignment {} holds indel pieces of the marginal passes next to pieces split from composed indels",
        self.alignment
      ),
    }
  }

  pub fn assemble_edge_mutations(
    &self,
    substitutions: &[Vec<Sub>],
    indels: &[Vec<InDel>],
    track: &MutationTrack,
  ) -> Result<Vec<Mutation>, Report> {
    combine_edge_mutations(
      self.assemble_substitutions(substitutions)?,
      &self.assemble_indels(indels)?,
      track,
    )
  }

  pub fn compose_indels(&self, parent: &[Vec<InDel>], child: &[Vec<InDel>]) -> Result<Vec<Vec<InDel>>, Report> {
    let mut parent = self.assemble_indels(parent)?;
    let mut child = self.assemble_indels(child)?;
    sort_indels(&mut parent);
    sort_indels(&mut child);
    self.split_indels(&compose_indels(&parent, &child))
  }

  pub fn split_indels(&self, indels: &[InDel]) -> Result<Vec<Vec<InDel>>, Report> {
    for indel in indels {
      self.check_alignment_indel(indel)?;
    }
    let split = match &self.columns {
      AlignmentColumns::Whole(_) => return Ok(vec![indels.to_vec()]),
      AlignmentColumns::Split(split) => split,
    };
    let mut pieces = vec![vec![]; split.partitions.len()];
    for (index, indel) in indels.iter().enumerate() {
      let split_from = Some(SplitFrom {
        index,
        range: indel.range,
      });
      let mut indel_pieces: BTreeMap<usize, InDel> = BTreeMap::new();
      for (site, &state) in izip!(&split.sites[indel.range.0..indel.range.1], indel.seq.as_slice()) {
        let piece = indel_pieces.entry(site.partition).or_insert_with(|| InDel {
          range: (site.local, site.local),
          seq: Seq::new(),
          kind: indel.kind,
          split_from,
        });
        piece.range.1 = site.local + 1;
        piece.seq.push(state);
      }
      for (partition, piece) in indel_pieces {
        pieces[partition].push(piece);
      }
    }
    Ok(pieces)
  }

  fn partition_count(&self) -> usize {
    match &self.columns {
      AlignmentColumns::Whole(_) => 1,
      AlignmentColumns::Split(split) => split.partitions.len(),
    }
  }

  fn length(&self) -> usize {
    match &self.columns {
      AlignmentColumns::Whole(whole) => whole.length,
      AlignmentColumns::Split(split) => split.sites.len(),
    }
  }

  fn partition_span(&self, partition: usize) -> (&str, usize) {
    match &self.columns {
      AlignmentColumns::Whole(whole) => (&whole.label, whole.length),
      AlignmentColumns::Split(split) => {
        let columns = &split.partitions[partition];
        (&columns.label, columns.columns.len())
      },
    }
  }

  fn column(&self, partition: usize, local: usize) -> Result<usize, Report> {
    match &self.columns {
      AlignmentColumns::Whole(whole) if local < whole.length => Ok(local),
      AlignmentColumns::Whole(whole) => Err(make_internal_report!(
        "Position {local} of {} is outside its {} columns",
        whole.label,
        whole.length
      )),
      AlignmentColumns::Split(split) => split.partitions[partition].column(local),
    }
  }

  fn check_piece(&self, partition: usize, piece: &InDel) -> Result<(), Report> {
    let (label, length) = self.partition_span(partition);
    let (lo, hi) = piece.range;
    if lo >= hi || hi > length || piece.seq.len() != hi - lo {
      return make_internal_error!("Indel {piece} does not fit {label}, which has {length} columns");
    }
    Ok(())
  }

  fn check_count(&self, count: usize, values: &str) -> Result<(), Report> {
    if count == self.partition_count() {
      Ok(())
    } else {
      make_internal_error!(
        "Alignment {} has {} partitions, but assembly received {count} {values}",
        self.alignment,
        self.partition_count()
      )
    }
  }

  fn check_alignment_indel(&self, indel: &InDel) -> Result<(), Report> {
    let (lo, hi) = indel.range;
    if let Some(split_from) = indel.split_from {
      return make_internal_error!(
        "Indel {indel} of alignment {} is a piece of the indel at {}--{}",
        self.alignment,
        split_from.range.0,
        split_from.range.1
      );
    }
    if lo >= hi || hi > self.length() || indel.seq.len() != hi - lo {
      return make_internal_error!(
        "Indel {indel} does not fit alignment {}, which has {} columns",
        self.alignment,
        self.length()
      );
    }
    Ok(())
  }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum AlignmentColumns {
  Whole(WholeColumns),
  Split(SplitColumns),
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct WholeColumns {
  label: String,
  length: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct SplitColumns {
  partitions: Vec<PartitionColumns>,
  sites: Vec<Site>,
}

impl SplitColumns {
  fn join_split_pieces(&self, alignment: usize, pieces: &[(usize, &InDel)]) -> Result<Vec<InDel>, Report> {
    let mut groups: BTreeMap<usize, (SplitFrom, InDelKind, Vec<(usize, AsciiChar)>)> = BTreeMap::new();
    for &(partition, piece) in pieces {
      let split_from = piece
        .split_from
        .ok_or_else(|| make_internal_report!("Indel {piece} is not a split piece"))?;
      let (group_split_from, group_kind, positions) = groups
        .entry(split_from.index)
        .or_insert_with(|| (split_from, piece.kind, vec![]));
      if (*group_split_from, *group_kind) != (split_from, piece.kind) {
        return make_internal_error!(
          "Indel pieces of alignment {alignment} disagree on indel {} of the split list: {:?} and {:?}",
          split_from.index,
          (*group_split_from, *group_kind),
          (split_from, piece.kind)
        );
      }
      positions.extend(self.partitions[partition].positions(piece));
    }
    groups
      .into_values()
      .map(|(split_from, kind, positions)| join_split_piece(alignment, split_from.range, kind, positions))
      .collect()
  }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct PartitionColumns {
  label: String,
  columns: Vec<usize>,
}

impl PartitionColumns {
  fn column(&self, local: usize) -> Result<usize, Report> {
    self.columns.get(local).copied().ok_or_else(|| {
      make_internal_report!(
        "Position {local} of {} is outside its {} columns",
        self.label,
        self.columns.len()
      )
    })
  }

  fn positions<'a>(&'a self, piece: &'a InDel) -> impl Iterator<Item = (usize, AsciiChar)> + 'a {
    izip!(&self.columns[piece.range.0..piece.range.1], piece.seq.as_slice()).map(|(&column, &state)| (column, state))
  }

  fn runs(&self, partition: usize, piece: &InDel) -> Vec<Run> {
    self.columns[piece.range.0..piece.range.1]
      .chunk_by(|column, next| column + 1 == *next)
      .scan(0, |offset, block| {
        let start = *offset;
        *offset += block.len();
        Some(Run {
          kind: piece.kind,
          range: (block[0], block[0] + block.len()),
          seq: Seq::from(&piece.seq.as_slice()[start..*offset]),
          partition,
        })
      })
      .collect()
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Site {
  partition: usize,
  local: usize,
}

fn join_split_piece(
  alignment: usize,
  range: (usize, usize),
  kind: InDelKind,
  mut positions: Vec<(usize, AsciiChar)>,
) -> Result<InDel, Report> {
  positions.sort_by_key(|&(position, _)| position);
  let (lo, hi) = range;
  if !positions.iter().map(|&(position, _)| position).eq(lo..hi) {
    return make_internal_error!(
      "The pieces of indel {lo}--{hi} ({kind:?}) of alignment {alignment} cover the positions {}",
      positions.iter().map(|&(position, _)| position).join(", ")
    );
  }
  Ok(InDel {
    range,
    seq: positions.into_iter().map(|(_, state)| state).collect(),
    kind,
    split_from: None,
  })
}

fn join_runs(runs: Vec<Run>) -> Vec<InDel> {
  runs
    .into_iter()
    .sorted_by_key(|run| (run.kind == InDelKind::Insertion, run.range))
    .coalesce(|mut run, next| {
      #[expect(
        clippy::suspicious_operation_groupings,
        reason = "runs join when the end of one equals the start of the next"
      )]
      let joins = run.kind == next.kind && run.range.1 == next.range.0 && run.partition != next.partition;
      if joins {
        run.range.1 = next.range.1;
        run.seq.extend(next.seq);
        run.partition = next.partition;
        Ok(run)
      } else {
        Err((run, next))
      }
    })
    .map(|run| InDel {
      range: run.range,
      seq: run.seq,
      kind: run.kind,
      split_from: None,
    })
    .collect()
}

struct Run {
  kind: InDelKind,
  range: (usize, usize),
  seq: Seq,
  partition: usize,
}
