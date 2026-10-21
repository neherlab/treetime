use crate::ancestral::__tests__::prop_generators::alignment::{arb_alignment, arb_alignment_no_gaps};
use crate::ancestral::__tests__::prop_generators::tree::{arb_tree_topology, taxa_names};
use crate::gtr::gtr::GTR;
use crate::test_utils::arb_gtr_nuc;
use proptest::prelude::*;
use std::collections::BTreeSet;
use treetime_graph::graph::Graph;
use treetime_io::nwk::nwk_read;
use treetime_primitives::AlignmentRecord;

#[derive(Debug, Clone)]
pub(crate) struct MarginalTestInput {
  pub(crate) newick: String,
  pub(crate) alignment: Vec<AlignmentRecord>,
  pub(crate) gtr: GTR,
  pub(crate) n_taxa: usize,
  pub(crate) seq_len: usize,
}

fn arb_marginal_input_with_params(n_taxa: usize, seq_len: usize) -> impl Strategy<Value = MarginalTestInput> {
  let taxa = taxa_names(n_taxa);
  let taxa_for_aln = taxa.clone();

  (
    arb_tree_topology(taxa),
    arb_alignment(taxa_for_aln, seq_len),
    arb_gtr_nuc(),
  )
    .prop_map(move |(tree, alignment, gtr)| {
      let newick = format!("({tree})root:0.001;");
      MarginalTestInput {
        newick,
        alignment,
        gtr,
        n_taxa,
        seq_len,
      }
    })
}

pub(crate) fn arb_marginal_input_no_gaps(n_taxa: usize, seq_len: usize) -> impl Strategy<Value = MarginalTestInput> {
  let taxa = taxa_names(n_taxa);
  let taxa_for_aln = taxa.clone();

  (
    arb_tree_topology(taxa),
    arb_alignment_no_gaps(taxa_for_aln, seq_len),
    arb_gtr_nuc(),
  )
    .prop_map(move |(tree, alignment, gtr)| {
      let newick = format!("({tree})root:0.001;");
      MarginalTestInput {
        newick,
        alignment,
        gtr,
        n_taxa,
        seq_len,
      }
    })
}

pub(crate) fn arb_marginal_input() -> impl Strategy<Value = MarginalTestInput> {
  (3_usize..=6, 5_usize..=20).prop_flat_map(|(n_taxa, seq_len)| arb_marginal_input_with_params(n_taxa, seq_len))
}

pub(crate) fn arb_marginal_input_small() -> impl Strategy<Value = MarginalTestInput> {
  (3_usize..=4, 3_usize..=10).prop_flat_map(|(n_taxa, seq_len)| arb_marginal_input_with_params(n_taxa, seq_len))
}

#[cfg(test)]
mod tests {
  use super::*;
  use proptest::proptest;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]

    #[test]
    fn test_prop_input_arb_marginal_input_valid(input in arb_marginal_input_small()) {
      prop_assert!(input.newick.ends_with(';'), "Invalid Newick: {}", input.newick);

      prop_assert_eq!(
        input.alignment.len(),
        input.n_taxa,
        "Alignment size mismatch: {} vs {}",
        input.alignment.len(),
        input.n_taxa
      );

      for record in &input.alignment {
        prop_assert_eq!(
          record.seq.len(),
          input.seq_len,
          "Sequence length mismatch for {}: {} vs {}",
          record.name,
          record.seq.len(),
          input.seq_len
        );
      }

      let pi_sum = input.gtr.pi().sum();
      prop_assert!(
        (pi_sum - 1.0).abs() < 1e-10,
        "GTR pi should sum to 1: {pi_sum}"
      );
    }

    #[test]
    fn test_prop_input_arb_marginal_input_taxa_match(input in arb_marginal_input_small()) {
      for record in &input.alignment {
        prop_assert!(
          input.newick.contains(&record.name),
          "Taxon {} not found in Newick: {}",
          record.name,
          input.newick
        );
      }
    }

    #[test]
    fn test_prop_input_arb_marginal_input_parseable_and_taxa_exact(input in arb_marginal_input_small()) {
      let nwk_parsed = nwk_read(input.newick.as_bytes()).unwrap();
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let graph: Graph = graph;

      let leaf_names: Vec<String> = graph.get_leaves().map(|leaf| names[leaf.key()].to_owned()).collect();

      let leaf_name_count = leaf_names.len();
      let leaf_names: BTreeSet<String> = leaf_names.into_iter().collect();
      let aln_names: BTreeSet<String> = input.alignment.iter().map(|record| record.name.clone()).collect();

      prop_assert_eq!(
        leaf_name_count,
        input.n_taxa,
        "Parsed tree leaf count should match n_taxa for {}",
        input.newick
      );
      prop_assert_eq!(
        leaf_names.len(),
        input.n_taxa,
        "Parsed tree leaves should be unique for {}",
        input.newick
      );
      prop_assert_eq!(
        leaf_names,
        aln_names,
        "Tree leaves and alignment taxa should match exactly for {}",
        input.newick
      );
    }
  }
}
