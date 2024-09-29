use crate::representation::{graph_sparse::SparseGraph, partitions_likelihood::PartitionLikelihood};
use eyre::Report;
use statrs::function;

// The likelihood of an edge length is the product of the likelihoods of all positions of all partitions
// Lh = prod_i prod_j \sum_{ab} s^{ij}_a exp(Q_i t)_{ab} r^{ij}_b
// The log likelihood is the sum of many terms
// logLh = sum_i sum_j \log(\sum_{ab} s^{ij}_a exp(Q_i t)_{ab} r^{ij}_b)
// to effectively calculate this, we need to reformulate the likelihood in terms of the eigenvectors of the GTR matrix
// Dropping the {ij} superscripts for brevity, we can write the likelihood as
// s_a exp(Qt)_{ab} r_b = s_a \sum_c v_{ac} exp(\lambda_c t) vinv_{cb} r_b = g_c exp(\lambda_c t) h_c = k_c exp(\lambda_c t)
// the k_c can be reused for different iterations of the branch length optimization
// logLh = sum_i sum_j \log(\sum_c k_c exp(\lambda^i_c t))
// the derivative is simply
// dlogLh/dt = sum_i sum_j \sum_c k_c \lambda_c exp(\lambda^i_c t) / \sum_c k_c exp(\lambda^i_c t)

// struct SiteContribution {
//   multiplicity: f64,
//   coefficients: Vec<f64>,
// }

// struct PartitionContribution {
//   site_contributions: Vec<siteContribution>,
//   eigenvalues: f64,
// }

// // function that takes two message projections, and gtr model, and the length of branch and returns the
// // likelihood as well as its derivative with respect to the branch length
// function likelihood_and_derivative(
//   source_coefficents: &[f64],
//   target_coefficents: &[f64],
//   gtr: &GTR,
//   length: f64,
// ) -> (f64, f64) {
//   let mut likelihood = 0.0;
//   let mut derivative = 0.0;
//   for (i, (source, target)) in source_coefficents.iter().zip(target_coefficents.iter()).enumerate() {
//     let mut res = 0.0;
//     for (j, (v, v_inv)) in gtr.v.iter().zip(gtr.v_inv.iter()).enumerate() {
//       res += v * source[j] + v_inv * target[j];
//     }
//     likelihood += res.exp();
//     derivative += res.exp() * gtr.v.iter().zip(gtr.v_inv.iter()).map(|(v, v_inv)| v * source[i] + v_inv * target[i]).sum();
//   }
//   likelihood = likelihood.ln();
//   derivative = derivative / likelihood;
//   likelihood = likelihood - length * derivative;
//   (likelihood, derivative)
// }

pub fn run_optimize_sparse(graph: &SparseGraph, partitions: &[PartitionLikelihood]) -> Result<(), Report> {
  graph.get_edges().iter_mut().for_each(|edge| {
    let mut edge = edge.write_arc().payload().write_arc();

    // let mut contributions: Vec<PartitionContribution> = Vec::new();
    // for each partition, sum the parts of the likelihood messages that favor short and long branches
    // for (pi, partition) in edge.sparse_partitions.iter().enumerate() {
    //   let PartitionLikelihood { gtr, length, alphabet } = &partitions[pi];
    //   // calculate projections of the 'msg_to_child' and 'msg_to_parent' onto the left and right eigenvectors
    //   let mut p_contribution = PartitionContribution {
    //     site_contributions: Vec::new(),
    //     eigenvalues: gtr.eigenvalues,
    //   };

    //   for (pos, site) in partition.msg_to_parent.variable.iter() {

    //   }
    //   source_coefficents.push(gtr.v.iter().map(|x| {
    //     let mut res = 0.0;
    //     for v in partition.msg_to_child.variable.values() {
    //       res +=  v.dis.dot(x);
    //     }
    //     for state in alphabet.canonical() {
    //       res += partition.msg_to_child.fixed[&state] * partition.msg_to_child.fixed_counts[&state].dot(x);
    //     }
    //     res
    //   }).collect());

    //   target_coefficents.push(gtr.v_inv.iter().map(|x| {
    //     let mut res = 0.0;
    //     for v in partition.msg_to_parent.variable.values() {
    //       res +=  v.dis.dot(x);
    //     }
    //     for state in alphabet.canonical() {
    //       res += partition.msg_to_parent.fixed[&state] * partition.msg_to_parent.fixed_counts[&state].dot(x);
    //     }
    //     res
    //   }).collect());

    //   // balance these optimally
    // edge.branch_length = Some(1.1 * edge.branch_length.unwrap_or(0.001));
  });
  Ok(())
}
