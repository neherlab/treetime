use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_with::skip_serializing_none;

/// Recurrent mutations of a `homoplasy` run: mutations that occur on more than one branch of the tree.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct HomoplasyResult {
  /// Whether positions count from 0 (`--zero-based`) instead of from 1.
  pub zero_based: bool,
  /// Whether `--drms` annotated the substitutions with drug resistance mutations.
  pub drm_annotated: bool,
  /// Substitutions between determined states (`A`, `C`, `G`, `T` for nucleotides).
  pub substitutions: SubstitutionResult,
  /// Changes from or to an ambiguous character, such as `N` or `R`.
  pub ambiguous: AmbiguousResult,
  /// Insertions and deletions.
  pub indels: IndelResult,
  /// Samples whose terminal branch carries recurrent or ambiguous changes, most homoplasic
  /// substitutions first.
  pub taxa: Vec<TaxonResult>,
}

/// Statistics of substitutions between determined states.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct SubstitutionResult {
  /// Number of sites of the genome: the alignment length plus `--const`.
  pub genome_length: usize,
  /// Sum of all branch lengths.
  pub total_branch_length: f64,
  /// Sum over terminal branches of the probability of an odd number of substitutions on a branch of
  /// length t, (1 - exp(-2t)) / 2.
  pub terminal_branch_length: f64,
  /// Substitutions on all branches.
  pub all: MutationTable,
  /// Substitutions on terminal branches.
  pub terminal: MutationTable,
  /// Number of sites by the number of substitutions at the site, with the expected number under a
  /// Poisson distribution with the same mean.
  pub site_hits: Vec<SiteHitsRow>,
  /// Log-likelihood of the site counts under the Poisson distribution minus its expected value under
  /// the same distribution. Negative values mean that substitutions cluster at fewer sites than the
  /// Poisson distribution predicts.
  pub log_likelihood_difference: f64,
}

/// Statistics of changes involving ambiguous characters.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct AmbiguousResult {
  /// Changes on all branches.
  pub all: MutationTable,
  /// Sites with changes, most branches first.
  pub sites: Vec<SiteBranchesRow>,
}

/// Statistics of insertions and deletions.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct IndelResult {
  /// Insertions and deletions on all branches.
  pub all: MutationTable,
  /// Insertions and deletions on terminal branches.
  pub terminal: MutationTable,
}

/// Mutations of one class with the number of branches each occurs on.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct MutationTable {
  /// Number of mutations, counted once per branch.
  pub mutations: usize,
  /// Number of distinct mutations by the number of branches they occur on.
  pub multiplicities: Vec<MultiplicityRow>,
  /// Distinct mutations, most branches first, then by position.
  pub ranked: Vec<RankedMutation>,
}

/// Number of distinct mutations that occur on the same number of branches.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct MultiplicityRow {
  /// Number of branches.
  pub branches: usize,
  /// Number of distinct mutations that occur on this many branches.
  pub mutations: usize,
}

/// A distinct mutation and the branches it occurs on.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct RankedMutation {
  /// The mutation: `G9343A` for a substitution, `del:100-102:ACG` or `ins:100-102:ACG` for a
  /// deletion or insertion of alignment columns 100 to 102.
  pub mutation: String,
  /// Number of branches the mutation occurs on.
  pub multiplicity: usize,
  /// Names of the nodes below the branches.
  pub branches: Vec<String>,
  /// Drug resistance annotation from `--drms`, for substitutions at a listed position.
  pub drm: Option<DrmAnnotation>,
}

/// Drug resistance annotation of a substitution.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct DrmAnnotation {
  /// Gene of the position.
  pub gene: String,
  /// Drug of the position.
  pub drug: String,
  /// Amino-acid substitution of the derived base, when the table lists the base.
  pub substitution: Option<String>,
}

/// Number of sites hit by the same number of substitutions.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct SiteHitsRow {
  /// Number of substitutions at a site.
  pub hits: usize,
  /// Number of sites with this many substitutions.
  pub sites: usize,
  /// Expected number of such sites under the Poisson distribution.
  pub expected: f64,
}

/// A site and the number of branches with a change at it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct SiteBranchesRow {
  /// Position of the site.
  pub position: usize,
  /// Number of branches with a change at the site.
  pub branches: usize,
}

/// Recurrent and ambiguous changes on the terminal branch of one sample.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct TaxonResult {
  /// Name of the sample.
  pub name: String,
  /// Substitutions on the terminal branch at sites with substitutions on two or more branches.
  pub homoplasic_mutations: Vec<String>,
  /// Number of the homoplasic substitutions at positions listed in `--drms`.
  pub drm_mutations: Option<usize>,
  /// Number of changes involving ambiguous characters on the terminal branch.
  pub ambiguous_changes: usize,
  /// Number of insertions and deletions on the terminal branch that occur on two or more branches.
  pub recurrent_indels: usize,
}
