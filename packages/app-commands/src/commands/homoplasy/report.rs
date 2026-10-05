use crate::commands::homoplasy::result::{DrmAnnotation, HomoplasyResult, MutationTable, TaxonResult};
use itertools::Itertools;
use treetime_utils::fmt::float::float_to_exponential;

pub fn render_homoplasy_report(result: &HomoplasyResult, n: usize, detailed: bool) -> String {
  let substitutions = &result.substitutions;
  let mut blocks = vec![multiplicity_block(
    &format!(
      "The TOTAL tree length is {} and {} mutations were observed.",
      float_to_exponential(substitutions.total_branch_length, 3),
      substitutions.all.mutations
    ),
    "mutations",
    &substitutions.all,
  )];
  if detailed {
    blocks.push(multiplicity_block(
      &format!(
        "The TERMINAL branch length is {} and {} mutations were observed.",
        float_to_exponential(substitutions.terminal_branch_length, 3),
        substitutions.terminal.mutations
      ),
      "mutations",
      &substitutions.terminal,
    ));
  }
  blocks.push(site_hits_block(result));
  blocks.push(format!(
    "log-likelihood difference to Poisson distribution with same mean: {}",
    float_to_exponential(substitutions.log_likelihood_difference, 3)
  ));
  blocks.push(top_mutations_block(
    &format!("\nThe {n} most homoplasic mutations are:"),
    "mut",
    &substitutions.all,
    n,
    result.drm_annotated,
  ));
  if detailed {
    blocks.push(top_mutations_block(
      &format!("\nThe {n} most homoplasic mutations on terminal branches are:"),
      "mut",
      &substitutions.terminal,
      n,
      result.drm_annotated,
    ));
    blocks.push(taxa_block(&result.taxa, n, result.drm_annotated));
  }
  blocks.push(multiplicity_block(
    &format!(
      "\nChanges involving ambiguous characters: {} were observed.",
      result.ambiguous.all.mutations
    ),
    "changes",
    &result.ambiguous.all,
  ));
  blocks.push(top_mutations_block(
    &format!("The {n} most frequent changes involving ambiguous characters are:"),
    "mut",
    &result.ambiguous.all,
    n,
    false,
  ));
  blocks.push(multiplicity_block(
    &format!(
      "\nInsertions and deletions: {} were observed.",
      result.indels.all.mutations
    ),
    "insertions and deletions",
    &result.indels.all,
  ));
  blocks.push(top_mutations_block(
    &format!("The {n} most frequent insertions and deletions are:"),
    "indel",
    &result.indels.all,
    n,
    false,
  ));
  blocks.iter().map(|block| format!("{block}\n")).join("\n")
}

fn multiplicity_block(header: &str, noun: &str, table: &MutationTable) -> String {
  let rows = table
    .multiplicities
    .iter()
    .map(|row| format!("\n\t - {} occur {} times", row.mutations, row.branches))
    .join("");
  format!("{header}\nOf these {} {noun},{rows}", table.mutations)
}

fn site_hits_block(result: &HomoplasyResult) -> String {
  let substitutions = &result.substitutions;
  let rows = substitutions
    .site_hits
    .iter()
    .filter(|row| row.sites > 0)
    .map(|row| {
      format!(
        "\n\t - {} were hit {} times (expected {:.2})",
        row.sites, row.hits, row.expected
      )
    })
    .join("");
  format!(
    "Of the {} positions in the genome,{rows}",
    substitutions.genome_length
  )
}

fn top_mutations_block(header: &str, label: &str, table: &MutationTable, n: usize, drms: bool) -> String {
  let drm_header = if drms { "\tDRM details (gene drug AAmut)" } else { "" };
  let rows = table
    .ranked
    .iter()
    .take(n)
    .take_while(|mutation| mutation.multiplicity > 1)
    .map(|mutation| {
      let drm = if drms {
        format!("\t{}", mutation.drm.as_ref().map(drm_details).unwrap_or_default())
      } else {
        String::new()
      };
      format!("\n\t{}\t{}{drm}", mutation.mutation, mutation.multiplicity)
    })
    .join("");
  format!("{header}\n\t{label}\tmultiplicity{drm_header}{rows}")
}

fn drm_details(drm: &DrmAnnotation) -> String {
  [Some(drm.gene.as_str()), Some(drm.drug.as_str()), drm.substitution.as_deref()]
    .into_iter()
    .flatten()
    .join(" ")
}

fn taxa_block(taxa: &[TaxonResult], n: usize, drms: bool) -> String {
  let drm_header = if drms { "\t# DRM" } else { "" };
  let rows = taxa
    .iter()
    .take(n)
    .map(|taxon| {
      let drm = taxon
        .drm_mutations
        .filter(|_| drms)
        .map(|count| format!("\t{count}"))
        .unwrap_or_default();
      format!(
        "\n\t{}\t{}{drm}\t{}\t{}",
        taxon.name,
        taxon.homoplasic_mutations.len(),
        taxon.ambiguous_changes,
        taxon.recurrent_indels
      )
    })
    .join("");
  format!(
    "\nTaxons that carry positions that mutated elsewhere in the tree:\n\ttaxon name\t#of homoplasic mutations{drm_header}\t# ambiguous changes\t# recurrent indels{rows}"
  )
}
