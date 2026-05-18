#[cfg(test)]
pub mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::initialize_marginal;
  use crate::gtr::gtr::GTR;
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  
  use parking_lot::RwLock;
  use std::sync::{Arc, LazyLock};
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  /// Felsenstein site likelihood for a two-taxon rooted tree: `(A:t1,B:t2)root;`.
  ///
  /// Formula (Felsenstein 1981, degenerate case with no internal nodes besides root):
  /// `L = sum_s pi[s] * P(obs_A | s, t1) * P(obs_B | s, t2)`.
  ///
  /// `P(i | j, t) = expQt[i, j]` is the transition probability from ancestral state `j`
  /// to descendant state `i` over branch length `t` in the column-stochastic convention.
  pub fn analytical_two_taxon_likelihood(gtr: &GTR, obs_a: usize, obs_b: usize, t1: f64, t2: f64) -> f64 {
    let exp_qt1 = gtr.expQt(t1);
    let exp_qt2 = gtr.expQt(t2);

    let mut likelihood = 0.0;
    for s in 0..gtr.pi.len() {
      likelihood += gtr.pi[s] * exp_qt1[[obs_a, s]] * exp_qt2[[obs_b, s]];
    }
    likelihood
  }

  /// Felsenstein site likelihood for a star tree `(A:t,B:t,C:t,D:t)root;`.
  ///
  /// Formula (generalization to multifurcation at root):
  /// `L = sum_s pi[s] * prod_i P(obs_i | s, t)`.
  pub fn analytical_star_tree_likelihood(gtr: &GTR, observations: &[usize], t: f64) -> f64 {
    let exp_qt = gtr.expQt(t);

    let mut likelihood = 0.0;
    for s in 0..gtr.pi.len() {
      let mut product = gtr.pi[s];
      for &obs in observations {
        product *= exp_qt[[obs, s]];
      }
      likelihood += product;
    }
    likelihood
  }

  /// Felsenstein site likelihood for a three-taxon tree `((A:t_a,B:t_b)AB:t_ab,C:t_c)root;`.
  ///
  /// Minimal topology with a non-trivial internal node, requiring the full recursive
  /// Felsenstein pruning pass.
  pub fn analytical_three_taxon_likelihood(
    gtr: &GTR,
    obs_a: usize,
    obs_b: usize,
    obs_c: usize,
    t_a: f64,
    t_b: f64,
    t_ab: f64,
    t_c: f64,
  ) -> f64 {
    let exp_qt_a = gtr.expQt(t_a);
    let exp_qt_b = gtr.expQt(t_b);
    let exp_qt_ab = gtr.expQt(t_ab);
    let exp_qt_c = gtr.expQt(t_c);

    let mut likelihood = 0.0;
    for s_root in 0..4 {
      for s_ab in 0..4 {
        let msg_ab = exp_qt_a[[obs_a, s_ab]] * exp_qt_b[[obs_b, s_ab]];
        likelihood += gtr.pi[s_root] * exp_qt_c[[obs_c, s_root]] * exp_qt_ab[[s_ab, s_root]] * msg_ab;
      }
    }
    likelihood
  }

  /// Map a nucleotide character to its index in the state vector: `A=0`, `C=1`, `G=2`, `T=3`.
  pub fn state_index(c: char) -> usize {
    match c {
      'A' => 0,
      'C' => 1,
      'G' => 2,
      'T' => 3,
      _ => unreachable!("Invalid nucleotide: {c}"),
    }
  }

  /// Run dense marginal ancestral reconstruction and return the total log-likelihood.
  ///
  /// Convenience wrapper: parse the Newick tree and FASTA alignment from strings,
  /// construct a single dense partition with the given GTR model, and run both passes via
  /// `initialize_marginal`.
  pub fn run_dense_marginal_get_log_lh(newick: &str, aln_str: &str, gtr: GTR) -> Result<f64, Report> {
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let aln = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let partitions = [Arc::new(RwLock::new(PartitionMarginalDense::new(0, gtr, alphabet, get_common_length(&aln)?)))];

    let log_lh = initialize_marginal(&graph, &partitions, &aln)?;
    Ok(log_lh)
  }
}
