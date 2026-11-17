#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::mixture::__tests__::fixtures::tests::helpers::{GAP_FREE, MISSING_DATA, fixture};
  use crate::partition::marginal::mixture::site_log_lh::{SparseMixture, SparseSiteLogLh};
  use crate::partition::marginal::shared::pass::indexed_backward_with_site_log_lh;
  use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate};
  use crate::pretty_assert_abs_diff_eq;
  use crate::test_utils::{arb_simulated_case, gtr_from_model, leaf_seq_inputs, make_nonuniform_gtr};
  use eyre::Report;
  use ndarray::array;
  use proptest::prelude::{ProptestConfig, proptest};
  use rstest::rstest;
  use treetime_graph::key_map::EdgeMap;
  use treetime_simulation::branch_lengths::exact;

  #[rustfmt::skip]
  #[rstest]
  #[case::gap_free(    &GAP_FREE)]
  #[case::missing_data(&MISSING_DATA)]
  #[trace]
  fn test_dense_site_log_lh_sums_to_the_pass_log_lh(#[case] alignment: &[(&str, &str)]) -> Result<(), Report> {
    let fixture = fixture(alignment)?;
    let gtr = make_nonuniform_gtr()?;
    let MarginalUpdate { log_lh, .. } =
      fixture.dense.marginal_update(&gtr, &fixture.graph, &fixture.branch_lengths, &())?;

    let (_, site_log_lh) = indexed_backward_with_site_log_lh(
      &fixture.dense,
      &gtr,
      &fixture.graph,
      &fixture.branch_lengths,
      fixture.dense.length,
    )?;
    let actual = site_log_lh.sum();

    pretty_assert_abs_diff_eq!(log_lh.value(), actual, epsilon = 1e-10);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::gap_free(    &GAP_FREE)]
  #[case::missing_data(&MISSING_DATA)]
  #[trace]
  fn test_sparse_site_log_lh_sums_to_the_pass_log_lh(#[case] alignment: &[(&str, &str)]) -> Result<(), Report> {
    let fixture = fixture(alignment)?;
    let gtr = make_nonuniform_gtr()?;
    let root = fixture.graph.root_key()?;
    let MarginalUpdate {
      node_states,
      edges,
      log_lh,
    } = fixture
        .sparse
        .marginal_update(&gtr, &fixture.graph, &fixture.branch_lengths, &fixture.sparse_nodes)?;

    let sites = SparseSiteLogLh::from_pass(&node_states[&root], &edges.backward, &fixture.sparse.alphabet);
    let composition = &fixture.sparse.obs_nodes[&root].composition;
    let actual = SparseMixture::new(&[sites], &array![1.0], composition, &fixture.sparse.alphabet)?.log_lh;

    pretty_assert_abs_diff_eq!(log_lh.value(), actual.value(), epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_dense_and_sparse_site_log_lh_agree_on_gap_free_data() -> Result<(), Report> {
    let fixture = fixture(&GAP_FREE)?;
    let gtr = make_nonuniform_gtr()?;
    let root = fixture.graph.root_key()?;
    let sparse =
      fixture
        .sparse
        .marginal_update(&gtr, &fixture.graph, &fixture.branch_lengths, &fixture.sparse_nodes)?;
    let (_, expected) = indexed_backward_with_site_log_lh(
      &fixture.dense,
      &gtr,
      &fixture.graph,
      &fixture.branch_lengths,
      fixture.dense.length,
    )?;

    let sites = SparseSiteLogLh::from_pass(
      &sparse.node_states[&root],
      &sparse.edges.backward,
      &fixture.sparse.alphabet,
    );
    let mixture = SparseMixture::new(
      &[sites],
      &array![1.0],
      &fixture.sparse.obs_nodes[&root].composition,
      &fixture.sparse.alphabet,
    )?;
    let root_sequence = &sparse.node_states[&root].sequence;
    let actual = (0..fixture.dense.length)
      .map(|pos| {
        mixture
          .site(pos, root_sequence[pos])
          .map_or(f64::NAN, |site| site.site_log_lh[0])
      })
      .collect::<ndarray::Array1<f64>>();

    pretty_assert_abs_diff_eq!(expected, actual, epsilon = 1e-10);
    Ok(())
  }

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    #[test]
    fn test_prop_sparse_site_log_lh_sums_to_the_pass_log_lh(spec in arb_simulated_case()) {
      let case = spec.simulate().expect("a case of the strategy simulates");
      let graph = &case.tree.graph;
      let lengths: EdgeMap<f64> = exact(&case.tree.edge_durations(), &case.model).expect("exact branch lengths");
      let records = case.sequences.leaf_records(graph, &case.tree.names);
      let alphabet = Alphabet::default();
      let fitch = create_fitch_partition(graph, alphabet, leaf_seq_inputs(graph, &case.tree.names, records))
        .expect("Fitch partition builds");
      let (partition, nodes) = fitch.into_marginal_sparse(graph).expect("sparse partition builds");
      let partition = partition.with_site_scales(true);
      let gtr = gtr_from_model(&case.model, 1.0).expect("model converts");
      let MarginalUpdate { node_states, edges, log_lh } =
        partition.marginal_update(&gtr, graph, &lengths, &nodes).expect("marginal update runs");
      let root = graph.root_key().expect("tree has a root");

      let sites = SparseSiteLogLh::from_pass(&node_states[&root], &edges.backward, &partition.alphabet);
      let actual = SparseMixture::new(&[sites], &array![1.0], &partition.obs_nodes[&root].composition, &partition.alphabet)
        .expect("mixture builds")
        .log_lh;

      pretty_assert_abs_diff_eq!(log_lh.value(), actual.value(), epsilon = 1e-8);
    }
  }
}
