use crate::alphabet::find_mutations::Mutation;
use crate::commands::ancestral::anc_args::TreetimeAncestralArgs;
use crate::commands::ancestral::anc_graph::AncestralGraph;
use crate::commands::ancestral::run_ancestral_reconstruction::run_ancestral_reconstruction;
use crate::commands::homoplasy::homoplasy_args::TreetimeHomoplasyArgs;
use crate::graph::node::GraphNodeKey;
use crate::utils::container::count_occurences;
use eyre::Report;
use itertools::Itertools;
use multimap::MultiMap;
use ordered_float::OrderedFloat;
use rstat::univariate::poisson::Poisson;
use rstat::DiscreteDistribution;
use std::collections::BTreeMap;
use std::hash::Hash;

pub fn run_homoplasy(homoplasy_args: TreetimeHomoplasyArgs) -> Result<(), Report> {
  let TreetimeHomoplasyArgs {
    ancestral_args,
    constant_sites,
    rescale,
    detailed,
    drms,
    num_mut,
  } = homoplasy_args;

  let TreetimeAncestralArgs {
    input_fastas,
    aln,
    vcf_reference,
    tree,
    alphabet,
    gtr,
    gtr_params,
    method_anc,
    aa,
    keep_overhangs,
    zero_based,
    reconstruct_tip_states,
    report_ambiguous,
    outdir,
    seed,
  } = &ancestral_args;

  if drms.is_some() {
    unimplemented!("DRMs are not yet implemented")
  }

  let (gtr, sequence_data, graph) = run_ancestral_reconstruction(&ancestral_args)?;

  let aln_len_full = sequence_data.len_full();

  let MappedMutations {
    positions,
    mutations,
    terminal_mutations,
  } = map_mutations(&graph);

  let mutations_by_strain = find_homoplasic_mutations_by_strain(&graph, &positions);

  let total_branch_length = calculate_total_branch_length(&graph);
  let corrected_branch_length = calculate_corrected_branch_length(&graph);
  let corrected_terminal_branch_length = calculate_corrected_terminal_branch_length(&graph);

  let expected_mutations = aln_len_full as f64 * corrected_branch_length;
  let expected_terminal_mutations = aln_len_full as f64 * corrected_terminal_branch_length;

  let (multiplicities, total_mutations) = calculate_histogram(&mutations);
  let (multiplicities_terminal, terminal_mutation_count) = calculate_histogram(&terminal_mutations);
  let (multiplicities_positions, positions_count) = calculate_histogram(&positions);

  println!("TOTAL tree length: {total_branch_length:.3}. Mutations observed: {total_mutations}.");
  println!("Of these {total_mutations} mutations,");
  for (n_occurences, n_muts) in &multiplicities {
    println!("  - {n_muts} occur {n_occurences} times");
  }

  println!();
  println!(
    "TERMINAL branch length: {corrected_terminal_branch_length:.3}. Mutations observed: {terminal_mutation_count}."
  );
  println!("Of these {total_mutations} mutations,");
  for (n_occurences, n_muts) in &multiplicities_terminal {
    println!("  - {n_muts} occur {n_occurences} times");
  }

  println!();
  println!("Of the {aln_len_full} positions in the genome,");
  for (n_occurences, n_muts) in &multiplicities_positions {
    let expected = aln_len_full as f64 * poisson_pmf(*n_occurences, total_mutations as f64 / aln_len_full as f64);
    println!("  - {n_muts} were hit {n_occurences} times (expected: {expected:.1})");
  }

  // TODO: calculate log-likelihood difference to Poisson distribution
  // p = poisson.pmf(np.arange(3*len(multiplicities_positions)),1.0*total_mutations/L);
  // diff = - L*np.sum(p*np.log(p+1e-100)) + np.sum(multiplicities_positions*np.log(p[_:len(multiplicities_positions)]+1e-100)))
  // println!("Log-likelihood difference to Poisson distribution with same mean: {diff:.3}");

  let muts = sort_mutations(&mutations, aln_len_full).take(num_mut);
  println!();
  println!("The {num_mut} most homoplasic mutations are:\n\tmut\tmultiplicity");
  for (m, count) in muts {
    println!("\t{m:<10}\t{count}");
  }

  println!();
  println!("The {num_mut} most homoplasic mutations on terminal branches are:\n\tmut\tmultiplicity");
  let muts = sort_mutations(&terminal_mutations, aln_len_full).take(num_mut);
  for (m, count) in muts {
    println!("\t{m:<10}\t{count}");
  }

  println!(
    "\n\nTaxons that carry positions that mutated elsewhere in the tree:\n\ttaxon name\t#of homoplasic mutations"
  );

  let muts = mutations_by_strain
    .into_iter()
    .map(|(node_key, mut_counter)| (node_key, mut_counter.len()))
    .sorted_by_key(|(node_key, count)| *count)
    .rev()
    .take(num_mut);

  for (node_key, count) in muts {
    let node = graph.get_node(node_key).expect("Node expected to exist, but not found");
    let node = node.read().payload();
    let node = node.read();
    let name = node.name.as_str();
    println!("\t{name}\t{count}");
  }

  Ok(())
}

fn sort_mutations(
  map: &MultiMap<Mutation, GraphNodeKey>,
  aln_len_full: usize,
) -> std::iter::Rev<std::vec::IntoIter<(&Mutation, usize)>> {
  map
    .iter_all()
    .map(|(mutation, node_keys)| (mutation, node_keys.len()))
    .sorted_by_key(|(mutation, count)| {
      OrderedFloat(*count as f64 - 0.1 * (mutation.pos as f64) / (aln_len_full as f64))
    })
    .rev()
}

/// Evaluate Poisson probability mass distribution function
fn poisson_pmf(k: usize, mu: f64) -> f64 {
  Poisson::new_unchecked(mu).pmf(&(k as u64)).into()
}

struct MappedMutations {
  pub positions: MultiMap<usize, GraphNodeKey>,
  pub mutations: MultiMap<Mutation, GraphNodeKey>,
  pub terminal_mutations: MultiMap<Mutation, GraphNodeKey>,
}

/// Gather graph node mutations and mutated positions into maps, associating mutations with nodes
fn map_mutations(graph: &AncestralGraph) -> MappedMutations {
  let mut positions = MultiMap::<usize, GraphNodeKey>::new();
  let mut mutations = MultiMap::<Mutation, GraphNodeKey>::new();
  let mut terminal_mutations = MultiMap::<Mutation, GraphNodeKey>::new();

  for node in graph.get_nodes() {
    let node = node.read();
    let key = node.key();
    let node = node.payload();
    let node = node.read();

    if node.is_root() {
      continue;
    }

    if !node.mutations.is_empty() {
      for mutation in &node.mutations {
        if mutation.is_valid() {
          positions.insert(mutation.pos, key);
          mutations.insert(mutation.clone(), key);
          if node.is_leaf() {
            terminal_mutations.insert(mutation.clone(), key);
          }
        }
      }
    }
  }

  MappedMutations {
    positions,
    mutations,
    terminal_mutations,
  }
}

/// Maps strains to homoplasic mutations they contain
fn find_homoplasic_mutations_by_strain(
  graph: &AncestralGraph,
  positions: &MultiMap<usize, GraphNodeKey>,
) -> BTreeMap<GraphNodeKey, Vec<(Mutation, usize)>> {
  let mut mutations_by_strain = BTreeMap::<GraphNodeKey, Vec<(Mutation, usize)>>::new();

  for node in graph.get_nodes() {
    let node = node.read();
    let key = node.key();
    let node = node.payload();
    let node = node.read();

    if !node.is_leaf() {
      continue;
    }

    if !node.mutations.is_empty() {
      for mutation in &node.mutations {
        let Mutation { reff, pos, qry } = &mutation;
        if mutation.is_valid() {
          if let Some(mutated_positions) = positions.get_vec(pos) {
            let num_mutations = mutated_positions.len();
            if num_mutations > 1 {
              mutations_by_strain
                .entry(key)
                .or_default()
                .push((mutation.clone(), num_mutations));
            }
          }
        }
      }
    }
  }

  mutations_by_strain
}

/// Total branch length is the expected number of substitutions
pub fn calculate_total_branch_length(graph: &AncestralGraph) -> f64 {
  graph
    .get_edges()
    .iter()
    .map(|edge| edge.read().payload().read().weight)
    .sum()
}

/// Corrected branch length is the expected number of observable substitutions
/// (probability of an odd number of substitutions at a particular site)
pub fn calculate_corrected_branch_length(graph: &AncestralGraph) -> f64 {
  graph
    .get_edges()
    .iter()
    .map(|edge| {
      let branch_length = edge.read().payload().read().weight;
      (-branch_length).exp() * branch_length.sinh()
    })
    .sum()
}

pub fn calculate_corrected_terminal_branch_length(graph: &AncestralGraph) -> f64 {
  graph
    .get_edges()
    .iter()
    .filter_map(|edge| {
      let edge = edge.read();

      let is_leaf = {
        let target = graph
          .get_node(edge.target())
          .expect("Graph node expected to exist, but not found");
        let target = target.read();
        target.is_leaf()
      };

      let branch_length = edge.payload().read().weight;

      is_leaf.then(|| (-branch_length).exp() * branch_length.sinh())
    })
    .sum()
}

pub fn calculate_histogram<K, V>(items: &MultiMap<K, V>) -> (Vec<(usize, usize)>, usize)
where
  K: Eq + Hash,
  V: Copy + Hash + Eq + Ord,
{
  let counts = items.iter_all().map(|(_, node_keys)| node_keys.len()).collect_vec();
  let multiplicities = count_occurences(counts.iter().copied());
  let total: usize = counts.iter().sum();
  (multiplicities, total)
}
