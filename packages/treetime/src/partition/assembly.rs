use crate::partition::resolve::PartitionLayout;
use crate::seq::indel::{InDel, InDelKind, SplitFrom, compose_indels, sort_indels};
use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
use crate::seq::sink::AmbiguousSites;
use crate::{make_internal_error, make_internal_report};
use eyre::Report;
use itertools::{Itertools, izip};
use std::collections::BTreeMap;
use treetime_primitives::{AsciiChar, Seq};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlignmentLayout {
  alignment: usize,
  partitions: Vec<PartitionColumns>,
  sites: Vec<Site>,
}

impl AlignmentLayout {
  pub fn new<'a>(length: usize, layouts: impl IntoIterator<Item = &'a PartitionLayout>) -> Result<Self, Report> {
    let layouts = layouts.into_iter().collect_vec();
    let Some(alignment) = layouts.first().map(|layout| layout.alignment) else {
      return make_internal_error!("An alignment layout needs at least one partition");
    };
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
      partitions,
      sites,
    })
  }

  pub fn is_whole(&self) -> bool {
    self.partitions.len() == 1
  }

  pub fn assemble_sequence(&self, sequences: &[Seq]) -> Result<Seq, Report> {
    self.check_count(sequences.len(), "sequences")?;
    for (partition, seq) in izip!(&self.partitions, sequences) {
      if seq.len() != partition.columns.len() {
        return make_internal_error!(
          "The sequence of {} has length {}, but the partition has {} columns",
          partition.label,
          seq.len(),
          partition.columns.len()
        );
      }
    }
    Ok(
      self
        .sites
        .iter()
        .map(|site| sequences[site.partition][site.local])
        .collect(),
    )
  }

  pub fn assemble_ambiguity(&self, ambiguous_sites: &[AmbiguousSites]) -> Result<AmbiguousSites, Report> {
    self.check_count(ambiguous_sites.len(), "ambiguity reports")?;
    izip!(&self.partitions, ambiguous_sites)
      .flat_map(|(partition, sites)| {
        sites
          .iter()
          .map(move |(&local, &code)| Ok((partition.column(local)?, code)))
      })
      .collect()
  }

  pub fn assemble_substitutions(&self, substitutions: &[Vec<Sub>]) -> Result<Vec<Sub>, Report> {
    self.check_count(substitutions.len(), "substitution lists")?;
    let mut assembled: Vec<Sub> = izip!(&self.partitions, substitutions)
      .flat_map(|(partition, subs)| {
        subs
          .iter()
          .map(move |sub| Sub::new(sub.reff(), partition.column(sub.pos())?, sub.qry()))
      })
      .try_collect()?;
    assembled.sort();
    Ok(assembled)
  }

  pub fn assemble_indels(&self, indels: &[Vec<InDel>]) -> Result<Vec<InDel>, Report> {
    self.check_count(indels.len(), "indel lists")?;
    for (columns, pieces) in izip!(&self.partitions, indels) {
      for piece in pieces {
        columns.check_piece(piece)?;
      }
    }
    if self.is_whole() {
      return Ok(indels[0].clone());
    }
    let pieces = izip!(0.., indels)
      .flat_map(|(partition, pieces)| pieces.iter().map(move |piece| (partition, piece)))
      .collect_vec();
    match pieces
      .iter()
      .map(|(_, piece)| piece.split_from.is_some())
      .all_equal_value()
    {
      Ok(true) => self.join_split_pieces(&pieces),
      Ok(false) => Ok(join_runs(
        pieces
          .iter()
          .flat_map(|&(partition, piece)| self.partitions[partition].runs(partition, piece))
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
    if self.is_whole() {
      return Ok(vec![indels.to_vec()]);
    }
    let mut pieces = vec![vec![]; self.partitions.len()];
    for (index, indel) in indels.iter().enumerate() {
      let split_from = Some(SplitFrom {
        index,
        range: indel.range,
      });
      let mut indel_pieces: BTreeMap<usize, InDel> = BTreeMap::new();
      for (site, &state) in izip!(&self.sites[indel.range.0..indel.range.1], indel.seq.as_slice()) {
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

  fn check_count(&self, count: usize, values: &str) -> Result<(), Report> {
    if count == self.partitions.len() {
      Ok(())
    } else {
      make_internal_error!(
        "Alignment {} has {} partitions, but assembly received {count} {values}",
        self.alignment,
        self.partitions.len()
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
    if lo >= hi || hi > self.sites.len() || indel.seq.len() != hi - lo {
      return make_internal_error!(
        "Indel {indel} does not fit alignment {}, which has {} columns",
        self.alignment,
        self.sites.len()
      );
    }
    Ok(())
  }

  fn join_split_pieces(&self, pieces: &[(usize, &InDel)]) -> Result<Vec<InDel>, Report> {
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
          "Indel pieces of alignment {} disagree on indel {} of the split list: {:?} and {:?}",
          self.alignment,
          split_from.index,
          (*group_split_from, *group_kind),
          (split_from, piece.kind)
        );
      }
      positions.extend(self.partitions[partition].positions(piece));
    }
    groups
      .into_values()
      .map(|(split_from, kind, positions)| self.join_split_piece(split_from.range, kind, positions))
      .collect()
  }

  fn join_split_piece(
    &self,
    range: (usize, usize),
    kind: InDelKind,
    mut positions: Vec<(usize, AsciiChar)>,
  ) -> Result<InDel, Report> {
    positions.sort_by_key(|&(position, _)| position);
    let (lo, hi) = range;
    if !positions.iter().map(|&(position, _)| position).eq(lo..hi) {
      return make_internal_error!(
        "The pieces of indel {lo}--{hi} ({kind:?}) of alignment {} cover the positions {}",
        self.alignment,
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

  fn check_piece(&self, piece: &InDel) -> Result<(), Report> {
    let (lo, hi) = piece.range;
    if lo >= hi || hi > self.columns.len() || piece.seq.len() != hi - lo {
      return make_internal_error!(
        "Indel {piece} does not fit {}, which has {} columns",
        self.label,
        self.columns.len()
      );
    }
    Ok(())
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
