use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::__tests__::prop_generators::alignment::{arb_alignment, arb_alignment_no_gaps};
use crate::ancestral::__tests__::prop_generators::tree::{arb_tree_topology, taxa_names};
use crate::gtr::gtr::GTR;
use ndarray::{Array1, Array2};
use proptest::prelude::*;
use std::collections::BTreeSet;
use treetime_graph::graph::Graph;
use treetime_io::nwk::nwk_read_str;
use treetime_primitives::AlignmentRecord;

fn arb_pi_nuc() -> impl Strategy<Value = Array1<f64>> {
  prop::collection::vec(0.0_f64..0.5, 4).prop_map(|offsets| {
    let bases = [0.5, 0.75, 1.0, 1.25];
    let raw: Vec<f64> = bases.iter().zip(offsets.iter()).map(|(b, o)| b + o).collect();
    let sum: f64 = raw.iter().sum();
    Array1::from_vec(raw.iter().map(|x| x / sum).collect())
  })
}

fn arb_w_nuc() -> impl Strategy<Value = Array2<f64>> {
  prop::collection::vec(0.0_f64..1.0, 6).prop_map(|offsets| {
    let bases = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5];
    let mut w = Array2::zeros((4, 4));
    let mut idx = 0;
    for i in 0..4 {
      for j in (i + 1)..4 {
        w[[i, j]] = bases[idx] + offsets[idx];
        w[[j, i]] = w[[i, j]];
        idx += 1;
      }
    }
    w
  })
}

fn arb_gtr_nuc() -> impl Strategy<Value = GTR> {
  (arb_pi_nuc(), arb_w_nuc(), 0.1_f64..5.0).prop_map(|(pi, w, mu)| {
    let alphabet = Alphabet::new(AlphabetName::Nuc).expect("Nuc alphabet should be valid");
    let n_states = alphabet.n_canonical();
    GTR::builder()
      .n_states(n_states)
      .mu(mu)
      .W(w)
      .pi(pi)
      .build()
      .expect("GTR construction should succeed with valid parameters")
  })
}

#[derive(Debug, Clone)]
pub(crate) struct MarginalTestInput {
  pub newick: String,
  pub alignment: Vec<AlignmentRecord>,
  pub gtr: GTR,
  pub n_taxa: usize,
  pub seq_len: usize,
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

      let pi_sum = input.gtr.pi.sum();
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
      let nwk_parsed = nwk_read_str(&input.newick).unwrap();
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let graph: Graph = graph;

      let mut leaf_names = Vec::new();
      for leaf in graph.get_leaves() {
        let maybe_name = names.get(&leaf.key()).cloned().flatten();
        prop_assert!(
          maybe_name.is_some(),
          "Leaf node is missing name in generated Newick: {}",
          input.newick
        );
        if let Some(name) = maybe_name {
          leaf_names.push(name);
        }
      }

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
