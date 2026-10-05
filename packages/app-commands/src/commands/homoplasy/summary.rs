use crate::commands::homoplasy::drms::DrmTable;
use crate::commands::homoplasy::result::{
  AmbiguousResult, DrmAnnotation, HomoplasyResult, IndelResult, MultiplicityRow, MutationTable, RankedMutation,
  SiteBranchesRow, SiteHitsRow, SubstitutionResult, TaxonResult,
};
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use treetime::homoplasy::pipeline::{HomoplasyOutput, IndelKey};
use treetime::homoplasy::recurrence::RecurrenceTable;
use treetime::seq::indel::InDelKind;
use treetime::seq::mutation::Sub;
use treetime_graph::assign_node_names::node_name_or_key;
use treetime_graph::node::GraphNodeKey;

pub struct ResultContext<'a> {
  pub names: &'a BTreeMap<GraphNodeKey, Option<String>>,
  pub drms: Option<&'a DrmTable>,
  pub zero_based: bool,
}

pub fn homoplasy_result(output: &HomoplasyOutput, context: &ResultContext<'_>) -> HomoplasyResult {
  let substitutions = &output.substitutions;
  let sub_name = |sub: &Sub| substitution_name(sub, context.zero_based);
  let sub_drm = |sub: &Sub| substitution_drm(sub, context.drms);
  let indel_name = |key: &IndelKey| indel_name(key, context.zero_based);
  HomoplasyResult {
    zero_based: context.zero_based,
    drm_annotated: context.drms.is_some(),
    substitutions: SubstitutionResult {
      genome_length: substitutions.sites.genome_length,
      total_branch_length: substitutions.total_branch_length,
      terminal_branch_length: substitutions.terminal_branch_length,
      all: mutation_table(&substitutions.all, context, sub_name, sub_drm),
      terminal: mutation_table(&substitutions.terminal, context, sub_name, sub_drm),
      site_hits: substitutions
        .sites
        .rows
        .iter()
        .map(|row| SiteHitsRow {
          hits: row.hits,
          sites: row.sites,
          expected: row.expected,
        })
        .collect(),
      log_likelihood_difference: substitutions.sites.log_likelihood_difference,
    },
    ambiguous: AmbiguousResult {
      all: mutation_table(&output.ambiguous.all, context, sub_name, |_| None),
      sites: output
        .ambiguous
        .sites
        .iter()
        .map(|site| SiteBranchesRow {
          position: site.position + offset(context.zero_based),
          branches: site.branches,
        })
        .collect(),
    },
    indels: IndelResult {
      all: mutation_table(&output.indels.all, context, indel_name, |_| None),
      terminal: mutation_table(&output.indels.terminal, context, indel_name, |_| None),
    },
    taxa: taxa(output, context),
  }
}

fn offset(zero_based: bool) -> usize {
  if zero_based { 0 } else { 1 }
}

fn substitution_name(sub: &Sub, zero_based: bool) -> String {
  format!("{}{}{}", sub.reff(), sub.pos() + offset(zero_based), sub.qry())
}

fn substitution_drm(sub: &Sub, drms: Option<&DrmTable>) -> Option<DrmAnnotation> {
  drms?.annotate(sub.pos(), &sub.qry().to_string())
}

fn indel_name(key: &IndelKey, zero_based: bool) -> String {
  let kind = match key.kind {
    InDelKind::Deletion => "del",
    InDelKind::Insertion => "ins",
  };
  let (start, end) = key.range;
  let offset = offset(zero_based);
  format!("{kind}:{}-{}:{}", start + offset, end - 1 + offset, key.sequence)
}

fn mutation_table<K>(
  table: &RecurrenceTable<K>,
  context: &ResultContext<'_>,
  name: impl Fn(&K) -> String,
  drm: impl Fn(&K) -> Option<DrmAnnotation>,
) -> MutationTable {
  MutationTable {
    mutations: table.count,
    multiplicities: table
      .histogram
      .iter()
      .map(|(&branches, &mutations)| MultiplicityRow { branches, mutations })
      .collect(),
    ranked: table
      .ranked
      .iter()
      .map(|recurrence| RankedMutation {
        mutation: name(&recurrence.mutation),
        multiplicity: recurrence.multiplicity(),
        branches: recurrence
          .branches
          .iter()
          .map(|&key| node_name_or_key(key, context.names[&key].as_deref()))
          .collect(),
        drm: drm(&recurrence.mutation),
      })
      .collect(),
  }
}

fn taxa(output: &HomoplasyOutput, context: &ResultContext<'_>) -> Vec<TaxonResult> {
  let leaves: BTreeSet<GraphNodeKey> = output
    .substitutions
    .leaves
    .keys()
    .chain(output.ambiguous.leaves.keys())
    .chain(output.indels.leaves.keys())
    .copied()
    .collect();
  leaves
    .into_iter()
    .map(|key| {
      let homoplasic = output
        .substitutions
        .leaves
        .get(&key)
        .map(Vec::as_slice)
        .unwrap_or_default();
      TaxonResult {
        name: node_name_or_key(key, context.names[&key].as_deref()),
        homoplasic_mutations: homoplasic
          .iter()
          .map(|sub| substitution_name(sub, context.zero_based))
          .collect(),
        drm_mutations: context
          .drms
          .map(|drms| homoplasic.iter().filter(|sub| drms.contains(sub.pos())).count()),
        ambiguous_changes: output.ambiguous.leaves.get(&key).copied().unwrap_or_default(),
        recurrent_indels: output.indels.leaves.get(&key).copied().unwrap_or_default(),
      }
    })
    .sorted_by(|a, b| {
      b.homoplasic_mutations
        .len()
        .cmp(&a.homoplasic_mutations.len())
        .then_with(|| a.name.cmp(&b.name))
    })
    .collect()
}
