use crate::commands::homoplasy::result::{
  DrmAnnotation, IndelResult, MultiplicityRow, RankedMutation, SiteBranchesRow, SiteHitsRow, SubstitutionResult,
  TaxonResult,
};
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_with::skip_serializing_none;
use std::collections::BTreeMap;

pub const AMBIGUOUS_SITES_SHOWN: usize = 100;

pub fn homoplasy_results(stats: Option<&HomoplasyStatsFile>) -> HomoplasyResults {
  HomoplasyResults {
    statistics: stats.map(homoplasy_statistics),
  }
}

pub fn homoplasy_statistics(stats: &HomoplasyStatsFile) -> HomoplasyStatistics {
  let substitutions = &stats.substitutions;
  let tree_offset = usize::from(stats.zero_based);
  let recurrent = recurrent_rows(&substitutions.all.ranked, &substitutions.terminal.ranked, tree_offset);
  let recurrent_indels = recurrent_rows(&stats.indels.all.ranked, &stats.indels.terminal.ranked, tree_offset);
  let sites = recurrent_sites(&substitutions.all.ranked, tree_offset);
  HomoplasyStatistics {
    drm_annotated: stats.drm_annotated,
    zero_based: stats.zero_based,
    genome_length: substitutions.genome_length,
    total_branch_length: substitutions.total_branch_length,
    substitutions: substitutions.all.mutations,
    distinct_substitutions: substitutions.all.ranked.len(),
    recurrent_substitutions: recurrent.len(),
    sites_hit_more_than_once: substitutions
      .site_hits
      .iter()
      .filter(|row| row.hits >= 2)
      .map(|row| row.sites)
      .sum(),
    expected_sites_hit_more_than_once: expected_sites_hit_more_than_once(
      substitutions.genome_length,
      &substitutions.site_hits,
    ),
    log_likelihood_difference: substitutions.log_likelihood_difference,
    samples_with_homoplasies: stats
      .taxa
      .iter()
      .filter(|taxon| !taxon.homoplasic_mutations.is_empty())
      .count(),
    recurrent_drm_substitutions: stats
      .drm_annotated
      .then(|| recurrent.iter().filter(|row| row.drm.is_some()).count()),
    ambiguous_changes: stats.ambiguous.all.mutations,
    indels: stats.indels.all.mutations,
    site_hits: substitutions.site_hits.clone(),
    multiplicities: substitutions.all.multiplicities.clone(),
    recurrent,
    sites,
    recurrent_indels,
    taxa: stats.taxa.clone(),
    ambiguous_sites: stats
      .ambiguous
      .sites
      .iter()
      .take(AMBIGUOUS_SITES_SHOWN)
      .map(|site| AmbiguousSite {
        position: site.position + tree_offset,
        display_position: site.position,
        branches: site.branches,
      })
      .collect(),
    ambiguous_site_count: stats.ambiguous.sites.len(),
  }
}

#[derive(Clone, Debug, PartialEq, Deserialize)]
pub struct HomoplasyStatsFile {
  pub zero_based: bool,
  pub drm_annotated: bool,
  pub substitutions: SubstitutionResult,
  pub ambiguous: AmbiguousStats,
  pub indels: IndelResult,
  pub taxa: Vec<TaxonResult>,
}

#[derive(Clone, Debug, PartialEq, Eq, Deserialize)]
pub struct AmbiguousStats {
  pub all: AmbiguousCount,
  pub sites: Vec<SiteBranchesRow>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Deserialize)]
pub struct AmbiguousCount {
  pub mutations: usize,
}

/// Results of a `homoplasy` run.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct HomoplasyResults {
  /// Statistics of the run; absent when its statistics file is missing or unreadable.
  pub statistics: Option<HomoplasyStatistics>,
}

/// Statistics of a `homoplasy` run: mutations that occur on more than one branch of the tree.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct HomoplasyStatistics {
  /// Whether `--drms` annotated the substitutions with drug resistance mutations.
  pub drm_annotated: bool,
  /// Whether the run counts positions from 0 (`--zero-based`). Display positions follow this
  /// setting; tree positions always count from 1.
  pub zero_based: bool,
  /// Number of sites of the genome: the alignment length plus the constant sites.
  pub genome_length: usize,
  /// Sum of all branch lengths.
  pub total_branch_length: f64,
  /// Substitutions between determined states on all branches, counted once per branch.
  pub substitutions: usize,
  /// Number of distinct substitutions.
  pub distinct_substitutions: usize,
  /// Number of distinct substitutions that occur on two or more branches.
  pub recurrent_substitutions: usize,
  /// Number of sites with substitutions on two or more branches, of any alleles.
  pub sites_hit_more_than_once: usize,
  /// Expected number of sites with two or more substitutions under a Poisson distribution with the
  /// same mean.
  pub expected_sites_hit_more_than_once: f64,
  /// Log-likelihood of the site counts under the Poisson distribution minus its expected value.
  /// Negative values mean that substitutions cluster at fewer sites than the Poisson distribution
  /// predicts.
  pub log_likelihood_difference: f64,
  /// Number of samples whose terminal branch has a substitution at a site hit more than once.
  pub samples_with_homoplasies: usize,
  /// Number of recurrent substitutions at positions listed in `--drms`; absent without `--drms`.
  pub recurrent_drm_substitutions: Option<usize>,
  /// Changes from or to an ambiguous character on all branches, counted once per branch.
  pub ambiguous_changes: usize,
  /// Insertions and deletions on all branches, counted once per branch.
  pub indels: usize,
  /// Number of sites by the number of substitutions at the site, with the Poisson expectation.
  pub site_hits: Vec<SiteHitsRow>,
  /// Number of distinct substitutions by the number of branches they occur on.
  pub multiplicities: Vec<MultiplicityRow>,
  /// Substitutions that occur on two or more branches, most branches first.
  pub recurrent: Vec<RecurrentMutation>,
  /// Sites with substitutions on two or more branches, by position.
  pub sites: Vec<HomoplasySite>,
  /// Insertions and deletions that occur on two or more branches, most branches first.
  pub recurrent_indels: Vec<RecurrentMutation>,
  /// Samples whose terminal branch carries homoplasic, ambiguous, or recurrent indel changes.
  pub taxa: Vec<TaxonResult>,
  /// Sites with changes involving ambiguous characters, most branches first, at most 100.
  pub ambiguous_sites: Vec<AmbiguousSite>,
  /// Number of sites with changes involving ambiguous characters.
  pub ambiguous_site_count: usize,
}

/// A mutation that occurs on two or more branches.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct RecurrentMutation {
  /// The mutation, such as `G9343A` or `del:100-102:ACG`.
  pub mutation: String,
  /// Position counted from 1, as in the tree: the site of a substitution, the first column of an
  /// insertion or deletion.
  pub position: usize,
  /// Position counted as the run reports positions.
  pub display_position: usize,
  /// Number of branches the mutation occurs on.
  pub branches: usize,
  /// Number of terminal branches the mutation occurs on.
  pub terminal_branches: usize,
  /// Names of the nodes below the branches.
  pub branch_names: Vec<String>,
  /// Drug resistance annotation from `--drms`.
  pub drm: Option<DrmAnnotation>,
}

/// A site with substitutions on two or more branches.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct HomoplasySite {
  /// Position counted from 1, as in the tree.
  pub position: usize,
  /// Position counted as the run reports positions.
  pub display_position: usize,
  /// Number of branches with a substitution at the site.
  pub branches: usize,
  /// Every substitution at the site, most branches first.
  pub substitutions: Vec<SiteSubstitution>,
}

/// A substitution at a site hit more than once.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct SiteSubstitution {
  /// The substitution, such as `G9343A`.
  pub mutation: String,
  /// Number of branches the substitution occurs on.
  pub branches: usize,
  /// Names of the nodes below the branches.
  pub branch_names: Vec<String>,
  /// Drug resistance annotation from `--drms`.
  pub drm: Option<DrmAnnotation>,
}

/// A site with changes involving ambiguous characters.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct AmbiguousSite {
  /// Position counted from 1, as in the tree.
  pub position: usize,
  /// Position counted as the run reports positions.
  pub display_position: usize,
  /// Number of branches with such a change at the site.
  pub branches: usize,
}

fn recurrent_rows(all: &[RankedMutation], terminal: &[RankedMutation], tree_offset: usize) -> Vec<RecurrentMutation> {
  let terminal: BTreeMap<&str, usize> = terminal
    .iter()
    .map(|row| (row.mutation.as_str(), row.multiplicity))
    .collect();
  all
    .iter()
    .filter(|row| row.multiplicity >= 2)
    .map(|row| RecurrentMutation {
      mutation: row.mutation.clone(),
      position: row.position + tree_offset,
      display_position: row.position,
      branches: row.multiplicity,
      terminal_branches: terminal.get(row.mutation.as_str()).copied().unwrap_or(0),
      branch_names: row.branches.clone(),
      drm: row.drm.clone(),
    })
    .collect()
}

fn recurrent_sites(all: &[RankedMutation], tree_offset: usize) -> Vec<HomoplasySite> {
  all
    .iter()
    .into_group_map_by(|row| row.position)
    .into_iter()
    .map(|(position, rows)| HomoplasySite {
      position: position + tree_offset,
      display_position: position,
      branches: rows.iter().map(|row| row.multiplicity).sum(),
      substitutions: rows
        .into_iter()
        .map(|row| SiteSubstitution {
          mutation: row.mutation.clone(),
          branches: row.multiplicity,
          branch_names: row.branches.clone(),
          drm: row.drm.clone(),
        })
        .collect(),
    })
    .filter(|site| site.branches >= 2)
    .sorted_by_key(|site| site.position)
    .collect()
}

#[expect(
  clippy::as_conversions,
  reason = "site counts are far below 2^53, so the conversion to f64 is exact"
)]
fn expected_sites_hit_more_than_once(genome_length: usize, site_hits: &[SiteHitsRow]) -> f64 {
  let expected_below_two: f64 = site_hits
    .iter()
    .filter(|row| row.hits < 2)
    .map(|row| row.expected)
    .sum();
  (genome_length as f64 - expected_below_two).max(0.0)
}
