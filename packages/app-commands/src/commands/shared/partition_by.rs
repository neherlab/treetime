use crate::commands::shared::argument_error::argument_error;
use crate::commands::shared::dependent_settings::check_applies;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use smart_default::SmartDefault;
use std::cmp::Reverse;
use std::collections::BTreeSet;
use std::path::PathBuf;
use treetime::alphabet::alphabet::AlphabetName;
use treetime::gtr::get_gtr::GtrModel;
use treetime::partition::resolve::PartitionEntry;
use treetime::partition::sites::site_ranges_from_columns;
use treetime::progress::{LogSink, RunWarning, RunWarningKind};
use treetime::progress_info;
use treetime_io::gff::{GffCdsAnnotation, GffCdsFeature, GffStrand};
use treetime_utils::cast::ToF64Exact;
use treetime_utils::fmt::string::{counted, name_list};
use treetime_utils::o;

const NONCODING: &str = "noncoding";

#[derive(Copy, Debug, Clone, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum PartitionScheme {
  /// Codon positions 1, 2 and 3 of every CDS of `--annotation` as the partitions `cp1`, `cp2` and `cp3`, and the
  /// other sites as `noncoding`
  CodonPositions,
  /// One partition per CDS of `--annotation`, named by the CDS, and the other sites as `noncoding`
  Genes,
  /// The partitions `slow`, `medium` and `fast` from the number of changes at each site in a Fitch reconstruction on
  /// the input tree, by the rule of `--rate-category-rule`
  RateCategories,
}

impl PartitionScheme {
  pub const fn uses_annotation(self) -> bool {
    matches!(self, Self::CodonPositions | Self::Genes)
  }

  pub const fn value_name(self) -> &'static str {
    match self {
      Self::CodonPositions => "codon-positions",
      Self::Genes => "genes",
      Self::RateCategories => "rate-categories",
    }
  }
}

#[derive(Copy, Debug, Clone, PartialEq, Eq, SmartDefault, JsonSchema, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum RateCategoryRule {
  /// `slow` for sites without changes, `medium` for one change, `fast` for two or more
  #[default]
  ChangeCounts,
  /// `slow` for sites without changes; changed sites at or below the median count of the changed sites go to
  /// `medium`, sites above it to `fast`
  MedianSplit,
  /// The three groups of change counts with the smallest sum of squared deviations from their group means (exact
  /// one-dimensional k-means), named by increasing mean count
  KMeans,
}

pub struct SchemeArguments<'a> {
  pub partition_by: Option<PartitionScheme>,
  pub rate_category_rule: Option<RateCategoryRule>,
  pub annotation: &'a [PathBuf],
  pub alignments: usize,
  pub alphabet: Option<AlphabetName>,
  pub translations: Option<bool>,
}

pub fn check_scheme_arguments(args: &SchemeArguments<'_>) -> Result<(), Report> {
  let scheme = args.partition_by;
  check_applies(
    "rate_category_rule",
    args.rate_category_rule.is_some(),
    "partition_by",
    scheme == Some(PartitionScheme::RateCategories),
    || o!("applies only with --partition-by rate-categories"),
  )?;
  let uses_annotation = scheme.is_some_and(PartitionScheme::uses_annotation);
  let annotation_given = !args.annotation.is_empty();
  match args.translations {
    None => check_applies("annotation", annotation_given, "partition_by", uses_annotation, || {
      o!("applies only with --partition-by codon-positions or genes, which read the CDSs of each alignment from it")
    })?,
    Some(translations) => check_applies(
      "annotation",
      annotation_given,
      "partition_by",
      uses_annotation || translations,
      || {
        o!(
          "applies only with --partition-by codon-positions or genes, which read the CDSs of each alignment from it, \
           or with --translations"
        )
      },
    )?,
  }

  let Some(scheme) = scheme else {
    return Ok(());
  };
  if args.alignments == 0 {
    return Err(argument_error(
      &["partition_by", "alignment"],
      "--partition-by needs at least one --alignment file, whose sites it puts into partitions",
    ));
  }
  if let Some(alphabet) = args.alphabet.filter(|alphabet| *alphabet != AlphabetName::Nuc) {
    return Err(argument_error(
      &["partition_by", "alphabet"],
      format!(
        "--partition-by needs the alphabet nuc, because CDS positions are nucleotide positions and the rate \
         categories count nucleotide changes, but the alphabet is {alphabet}"
      ),
    ));
  }
  if scheme.uses_annotation() && args.annotation.len() != args.alignments {
    return Err(argument_error(
      &["annotation", "partition_by", "alignment"],
      format!(
        "--partition-by {} needs one --annotation file per --alignment file, matched in list order, but the run has \
         {} and {}",
        scheme.value_name(),
        counted(args.alignments, "--alignment file"),
        counted(args.annotation.len(), "--annotation file")
      ),
    ));
  }
  Ok(())
}

pub struct AnnotatedAlignment<'a> {
  pub length: usize,
  pub annotation: &'a GffCdsAnnotation,
}

pub fn codon_position_entries(
  alignments: &[AnnotatedAlignment<'_>],
  model: &GtrModel,
  log: &dyn LogSink,
) -> Result<Vec<PartitionEntry>, Report> {
  let n_alignments = alignments.len();
  let entries = alignments
    .iter()
    .zip(1..)
    .map(|(alignment, number)| {
      let owners = owned_sites(alignment, number)?;
      warn_cdses_without_sites(&owners, number, n_alignments, log);
      let groups = (1..=3)
        .map(|position| {
          let columns = owners
            .sites
            .iter()
            .enumerate()
            .filter(|(_, owner)| owner.is_some_and(|owner| owner.codon_position == position))
            .map(|(column, _)| column)
            .collect::<BTreeSet<usize>>();
          (format!("cp{position}"), columns)
        })
        .collect_vec();
      Ok(scheme_entries(
        groups,
        owners.has_noncoding(),
        number,
        n_alignments,
        model,
      ))
    })
    .collect::<Result<Vec<_>, Report>>()?
    .concat();
  check_scheme_entries(&entries)?;
  Ok(entries)
}

pub fn gene_entries(
  alignments: &[AnnotatedAlignment<'_>],
  model: &GtrModel,
  log: &dyn LogSink,
) -> Result<Vec<PartitionEntry>, Report> {
  let n_alignments = alignments.len();
  let entries = alignments
    .iter()
    .zip(1..)
    .map(|(alignment, number)| {
      let owners = owned_sites(alignment, number)?;
      warn_cdses_without_sites(&owners, number, n_alignments, log);
      let groups = owners
        .sites
        .iter()
        .enumerate()
        .filter_map(|(column, owner)| owner.map(|owner| (owner.cds, column)))
        .into_group_map()
        .into_iter()
        .map(|(cds, columns)| {
          (
            owners.ranked[cds].label.clone(),
            columns.into_iter().collect::<BTreeSet<usize>>(),
          )
        })
        .sorted_by_key(|(_, columns)| columns.first().copied())
        .collect_vec();
      Ok(scheme_entries(
        groups,
        owners.has_noncoding(),
        number,
        n_alignments,
        model,
      ))
    })
    .collect::<Result<Vec<_>, Report>>()?
    .concat();
  check_scheme_entries(&entries)?;
  Ok(entries)
}

pub fn rate_category_entries(
  rule: RateCategoryRule,
  counts: &[Vec<usize>],
  model: &GtrModel,
  log: &dyn LogSink,
) -> Result<Vec<PartitionEntry>, Report> {
  let n_alignments = counts.len();
  let entries = counts
    .iter()
    .zip(1..)
    .flat_map(|(counts, number)| {
      let buckets = rate_category_buckets(rule, counts);
      for bucket in &buckets {
        progress_info!(
          log,
          "Partition {}: {} sites with {} to {} changes",
          entry_name(&bucket.name, number, n_alignments),
          bucket.columns.len(),
          bucket.min_changes,
          bucket.max_changes
        );
      }
      let groups = buckets
        .into_iter()
        .map(|bucket| (bucket.name, bucket.columns))
        .collect_vec();
      scheme_entries(groups, false, number, n_alignments, model)
    })
    .collect_vec();
  check_scheme_entries(&entries)?;
  Ok(entries)
}

pub fn rate_category_buckets(rule: RateCategoryRule, counts: &[usize]) -> Vec<RateBucket> {
  let thresholds = match rule {
    RateCategoryRule::ChangeCounts => vec![0, 1],
    RateCategoryRule::MedianSplit => median_split_thresholds(counts),
    RateCategoryRule::KMeans => k_means_thresholds(counts),
  };
  let bucket_of = |count: usize| thresholds.iter().take_while(|&&threshold| count > threshold).count();
  let by_bucket = counts
    .iter()
    .enumerate()
    .map(|(column, &count)| (bucket_of(count), (column, count)))
    .into_group_map();
  ["slow", "medium", "fast"]
    .into_iter()
    .enumerate()
    .filter_map(|(bucket, name)| {
      let sites = by_bucket.get(&bucket)?;
      Some(RateBucket {
        name: name.to_owned(),
        columns: sites.iter().map(|&(column, _)| column).collect(),
        min_changes: sites.iter().map(|&(_, count)| count).min()?,
        max_changes: sites.iter().map(|&(_, count)| count).max()?,
      })
    })
    .collect()
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RateBucket {
  pub name: String,
  pub columns: BTreeSet<usize>,
  pub min_changes: usize,
  pub max_changes: usize,
}

pub fn k_means_thresholds(counts: &[usize]) -> Vec<usize> {
  let distinct = counts.iter().copied().counts().into_iter().sorted().collect_vec();
  if distinct.len() <= 3 {
    return distinct.iter().map(|&(value, _)| value).take(2).collect();
  }
  let prefix = distinct
    .iter()
    .scan(PrefixSums::default(), |sums, &(value, weight)| {
      *sums = sums.add(value, weight);
      Some(*sums)
    })
    .collect_vec();
  let cost = |first: usize, last: usize| {
    let before = first
      .checked_sub(1)
      .map_or_else(PrefixSums::default, |index| prefix[index]);
    prefix[last].sub(before).squared_deviations()
  };
  let n = distinct.len();
  let (first_end, second_end) = (0..n - 2)
    .flat_map(|first_end| (first_end + 1..n - 1).map(move |second_end| (first_end, second_end)))
    .min_by(|&(a1, a2), &(b1, b2)| {
      let total = |end1: usize, end2: usize| cost(0, end1) + cost(end1 + 1, end2) + cost(end2 + 1, n - 1);
      total(a1, a2).total_cmp(&total(b1, b2))
    })
    .unwrap_or((0, 1));
  vec![distinct[first_end].0, distinct[second_end].0]
}

fn median_split_thresholds(counts: &[usize]) -> Vec<usize> {
  let changed = counts.iter().copied().filter(|&count| count > 0).sorted().collect_vec();
  let Some(last) = changed.len().checked_sub(1) else {
    return vec![0];
  };
  let lower_middle = changed[last.div_euclid(2)];
  let upper_middle = changed[changed.len().div_euclid(2)];
  vec![0, usize::midpoint(lower_middle, upper_middle)]
}

#[derive(Clone, Copy, Debug, Default)]
struct PrefixSums {
  weight: f64,
  sum: f64,
  sum_squares: f64,
}

impl PrefixSums {
  fn add(self, value: usize, weight: usize) -> Self {
    let (value, weight) = (value.to_f64_exact(), weight.to_f64_exact());
    Self {
      weight: self.weight + weight,
      sum: weight.mul_add(value, self.sum),
      sum_squares: (weight * value).mul_add(value, self.sum_squares),
    }
  }

  fn sub(self, other: Self) -> Self {
    Self {
      weight: self.weight - other.weight,
      sum: self.sum - other.sum,
      sum_squares: self.sum_squares - other.sum_squares,
    }
  }

  fn squared_deviations(self) -> f64 {
    self.sum_squares - self.sum * self.sum / self.weight
  }
}

fn scheme_entries(
  groups: Vec<(String, BTreeSet<usize>)>,
  noncoding: bool,
  alignment: usize,
  n_alignments: usize,
  model: &GtrModel,
) -> Vec<PartitionEntry> {
  let explicit = groups
    .into_iter()
    .filter(|(_, columns)| !columns.is_empty())
    .map(|(name, columns)| PartitionEntry {
      name: Some(entry_name(&name, alignment, n_alignments)),
      alignment: Some(alignment),
      sites: Some(site_ranges_from_columns(&columns)),
      model: model.clone(),
    });
  let rest = noncoding.then(|| PartitionEntry {
    name: Some(entry_name(NONCODING, alignment, n_alignments)),
    alignment: Some(alignment),
    sites: None,
    model: model.clone(),
  });
  explicit.chain(rest).collect()
}

fn entry_name(name: &str, alignment: usize, n_alignments: usize) -> String {
  if n_alignments > 1 {
    format!("{alignment}/{name}")
  } else {
    name.to_owned()
  }
}

fn check_scheme_entries(entries: &[PartitionEntry]) -> Result<(), Report> {
  if entries.len() >= 2 {
    return Ok(());
  }
  let names = entries
    .iter()
    .filter_map(|entry| entry.name.as_deref())
    .map(|name| format!(" ('{name}')"))
    .join("");
  Err(argument_error(
    &["partition_by"],
    format!(
      "--partition-by makes only {}{names} over all --alignment files, but a partition scheme needs at least two \
       partitions",
      counted(entries.len(), "partition")
    ),
  ))
}

fn owned_sites<'a>(alignment: &AnnotatedAlignment<'a>, number: usize) -> Result<SiteOwners<'a>, Report> {
  check_annotation(alignment, number)?;
  let ranked = alignment
    .annotation
    .features
    .iter()
    .sorted_by_key(|cds| (Reverse(cds.coding_length()), first_site(cds), cds.label.clone()))
    .collect_vec();
  let sites = ranked.iter().enumerate().fold(
    vec![None; alignment.length],
    |mut sites: Vec<Option<SiteOwner>>, (cds, feature)| {
      for segment in &feature.segments {
        for site in segment.start..=segment.end {
          let owner = &mut sites[site - 1];
          if owner.is_none() {
            *owner = Some(SiteOwner {
              cds,
              codon_position: codon_position(feature.strand, segment.start, segment.end, segment.phase.offset(), site),
            });
          }
        }
      }
      sites
    },
  );
  Ok(SiteOwners { ranked, sites })
}

pub const fn codon_position(strand: GffStrand, start: usize, end: usize, phase: usize, site: usize) -> usize {
  let from_five_prime = match strand {
    GffStrand::Forward => site - start,
    GffStrand::Reverse => end - site,
  };
  (from_five_prime + 3 - phase) % 3 + 1
}

fn first_site(cds: &GffCdsFeature) -> usize {
  cds
    .segments
    .iter()
    .map(|segment| segment.start)
    .min()
    .unwrap_or(usize::MAX)
}

fn check_annotation(alignment: &AnnotatedAlignment<'_>, number: usize) -> Result<(), Report> {
  let settings = ["annotation", "alignment"];
  let length = alignment.length;
  let seqids = alignment
    .annotation
    .features
    .iter()
    .map(|cds| cds.seqid.as_str())
    .collect::<BTreeSet<&str>>();
  let seqid = match seqids.iter().exactly_one() {
    Ok(seqid) => *seqid,
    Err(_) => {
      return Err(argument_error(
        &settings,
        format!(
          "The annotation of alignment {number} places its CDSs on the sequences {}, but one annotation must describe \
           the positions of one alignment, so all its CDSs need the same sequence ID",
          seqids.iter().map(|seqid| format!("'{seqid}'")).join(", ")
        ),
      ));
    },
  };

  let regions = alignment
    .annotation
    .sequence_regions
    .iter()
    .filter(|region| region.seqid == seqid)
    .collect_vec();
  if regions.is_empty() {
    return Err(argument_error(
      &settings,
      format!(
        "The annotation of alignment {number} has no ##sequence-region line for the sequence '{seqid}': TreeTime \
         needs '##sequence-region {seqid} 1 {length}' to check that the annotation positions are the {length} columns \
         of the alignment"
      ),
    ));
  }
  if let Some(region) = regions.iter().find(|region| region.start != 1 || region.end != length) {
    return Err(argument_error(
      &settings,
      format!(
        "The annotation of alignment {number} gives the sequence '{seqid}' the region {}-{}, but the alignment has \
         {length} columns: annotation positions are alignment columns only for an alignment on that reference, \
         without insertion columns, so the region must run from 1 to {length}",
        region.start, region.end
      ),
    ));
  }
  if let Some((cds, end)) = alignment
    .annotation
    .features
    .iter()
    .flat_map(|cds| cds.segments.iter().map(move |segment| (cds, segment.end)))
    .find(|&(_, end)| end > length)
  {
    return Err(argument_error(
      &settings,
      format!(
        "The annotation of alignment {number} has CDS '{}' with a segment that ends at position {end}, after the \
         last column {length} of the alignment",
        cds.label
      ),
    ));
  }
  Ok(())
}

fn warn_cdses_without_sites(owners: &SiteOwners<'_>, number: usize, n_alignments: usize, log: &dyn LogSink) {
  let with_sites = owners
    .sites
    .iter()
    .flatten()
    .map(|owner| owner.cds)
    .collect::<BTreeSet<usize>>();
  let names = owners
    .ranked
    .iter()
    .enumerate()
    .filter(|(cds, _)| !with_sites.contains(cds))
    .map(|(_, cds)| entry_name(&cds.label, number, n_alignments))
    .sorted()
    .collect_vec();
  if names.is_empty() {
    return;
  }
  log.warning(&RunWarning {
    kind: RunWarningKind::CdsWithoutPartitionSites,
    message: format!(
      "--partition-by makes no partition for CDS {}: each of their sites is also covered by a longer CDS, and the \
       longest CDS that covers a site owns it.",
      name_list(&names)
    ),
    names,
  });
}

struct SiteOwners<'a> {
  ranked: Vec<&'a GffCdsFeature>,
  sites: Vec<Option<SiteOwner>>,
}

impl SiteOwners<'_> {
  fn has_noncoding(&self) -> bool {
    self.sites.iter().any(Option::is_none)
  }
}

#[derive(Clone, Copy, Debug)]
struct SiteOwner {
  cds: usize,
  codon_position: usize,
}
