#[cfg(test)]
pub(super) mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::gtr::gtr::GTR;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use eyre::Report;
  use treetime_graph::graph::Graph;

  use std::sync::LazyLock;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  pub(crate) fn analytical_two_taxon_likelihood(gtr: &GTR, obs_a: usize, obs_b: usize, t1: f64, t2: f64) -> f64 {
    let exp_qt1 = gtr.expQt(t1);
    let exp_qt2 = gtr.expQt(t2);

    let mut likelihood = 0.0;
    for s in 0..gtr.pi.len() {
      likelihood += gtr.pi[s] * exp_qt1[[obs_a, s]] * exp_qt2[[obs_b, s]];
    }
    likelihood
  }

  pub(crate) fn analytical_star_tree_likelihood(gtr: &GTR, observations: &[usize], t: f64) -> f64 {
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

  pub(crate) fn analytical_three_taxon_likelihood(
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

  pub(crate) fn state_index(c: char) -> usize {
    match c {
      'A' => 0,
      'C' => 1,
      'G' => 2,
      'T' => 3,
      _ => unreachable!("Invalid nucleotide: {c}"),
    }
  }

  pub(crate) fn run_dense_marginal_get_log_lh(newick: &str, aln_str: &str, gtr: GTR) -> Result<f64, Report> {
    let nwk_parsed = nwk_read_str(newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(aln_str, &*NUC_ALPHABET)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let partition = PartitionMarginalDense::new(0, alphabet, get_common_length(&aln)?);
    let node_states = partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, aln))?;
    let recon = DenseReconstruction::seeded(partition, gtr, node_states);
    let (recon, log_lh) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh = log_lh.value();
    Ok(log_lh)
  }
}
