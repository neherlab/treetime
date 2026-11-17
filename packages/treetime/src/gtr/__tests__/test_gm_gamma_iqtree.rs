#[cfg(test)]
mod tests {
  use crate::gtr::gamma::GammaRates;
  use crate::gtr::partition_model::PartitionModel;
  use crate::pretty_assert_abs_diff_eq;
  use eyre::Report;
  use helpers::{IqtreeInputs, IqtreeOutputs, Representation, TreetimeCase};
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( Representation::Dense)]
  #[case::sparse(Representation::Sparse)]
  #[trace]
  fn test_site_log_lh_without_rate_variation_matches_iqtree(#[case] representation: Representation) -> Result<(), Report> {
    let inputs = IqtreeInputs::load()?;
    let iqtree = IqtreeOutputs::load()?;
    let case = TreetimeCase::new(&inputs)?;
    let model = PartitionModel::from(inputs.gtr()?);

    let actual = case.site_log_lh(&model, representation)?;

    pretty_assert_abs_diff_eq!(iqtree.no_gamma.site_log_lh, actual, epsilon = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense_shape_0_2( Representation::Dense,  0)]
  #[case::dense_shape_1(   Representation::Dense,  1)]
  #[case::dense_shape_5(   Representation::Dense,  2)]
  #[case::sparse_shape_0_2(Representation::Sparse, 0)]
  #[case::sparse_shape_1(  Representation::Sparse, 1)]
  #[case::sparse_shape_5(  Representation::Sparse, 2)]
  #[trace]
  fn test_gamma_mixture_site_log_lh_matches_iqtree(
    #[case] representation: Representation,
    #[case] shape_index: usize,
  ) -> Result<(), Report> {
    let inputs = IqtreeInputs::load()?;
    let iqtree = IqtreeOutputs::load()?;
    let iqtree_gamma = &iqtree.gamma[shape_index];
    let case = TreetimeCase::new(&inputs)?;
    let model = PartitionModel::new(inputs.gtr()?, GammaRates::new(inputs.categories, iqtree_gamma.shape)?);

    let actual = case.site_log_lh(&model, representation)?;

    pretty_assert_abs_diff_eq!(iqtree_gamma.site_log_lh, actual, epsilon = 1e-6);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::shape_0_2(0)]
  #[case::shape_1(  1)]
  #[case::shape_5(  2)]
  #[trace]
  fn test_gamma_category_rates_match_iqtree(#[case] shape_index: usize) -> Result<(), Report> {
    let inputs = IqtreeInputs::load()?;
    let iqtree = IqtreeOutputs::load()?;
    let iqtree_gamma = &iqtree.gamma[shape_index];

    let actual = GammaRates::new(inputs.categories, iqtree_gamma.shape)?;

    pretty_assert_abs_diff_eq!(iqtree_gamma.category_rates, actual.rates(), epsilon = 1e-8);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::gtr::gtr::{GTR, Mu};
    use crate::gtr::partition_model::PartitionModel;
    use crate::partition::fitch::passes::create_fitch_partition;
    use crate::partition::marginal::dense::partition::PartitionMarginalDense;
    use crate::partition::marginal::mixture::dense::dense_mixture_update;
    use crate::partition::marginal::mixture::site_log_lh::mix_site;
    use crate::partition::marginal::mixture::sparse::sparse_mixture_update;
    use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
    use crate::partition::storage::sparse::SparseNodeState;
    use crate::test_utils::leaf_seq_inputs;
    use deser::Deserialize;
    use eyre::Report;
    use ndarray::{Array1, Array2};
    use std::path::Path;
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::{EdgeMap, NodeMap};
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, Seq};
    use treetime_utils::adapters::{Array2Rows, ArrayVec};
    use treetime_utils::io::json::json_read_file;
    use treetime_utils::make_internal_report;

    const FIXTURES_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/src/gtr/__tests__/__fixtures__");

    #[derive(Debug, Deserialize)]
    pub(crate) struct IqtreeInputs {
      tree: String,
      #[deser(rename = "W", deserialize_as = Array2Rows)]
      w: Array2<f64>,
      #[deser(deserialize_as = ArrayVec)]
      pi: Array1<f64>,
      mu: f64,
      pub(crate) categories: usize,
      alignment: Vec<IqtreeRecord>,
    }

    impl IqtreeInputs {
      pub(crate) fn load() -> Result<Self, Report> {
        json_read_file(Path::new(FIXTURES_DIR).join("gm_gamma_iqtree_inputs.json"))
      }

      pub(crate) fn gtr(&self) -> Result<GTR, Report> {
        GTR::builder()
          .n_states(Alphabet::default().n_canonical())
          .mu(Mu::Normalized(self.mu))
          .W(self.w.clone())
          .pi(self.pi.clone())
          .build()
      }
    }

    #[derive(Debug, Deserialize)]
    struct IqtreeRecord {
      name: String,
      seq: String,
    }

    #[derive(Debug, Deserialize)]
    pub(crate) struct IqtreeOutputs {
      pub(crate) no_gamma: IqtreeSiteLogLh,
      pub(crate) gamma: Vec<IqtreeGamma>,
    }

    impl IqtreeOutputs {
      pub(crate) fn load() -> Result<Self, Report> {
        json_read_file(Path::new(FIXTURES_DIR).join("gm_gamma_iqtree_outputs.json"))
      }
    }

    #[derive(Debug, Deserialize)]
    pub(crate) struct IqtreeSiteLogLh {
      #[deser(deserialize_as = ArrayVec)]
      pub(crate) site_log_lh: Array1<f64>,
    }

    #[derive(Debug, Deserialize)]
    pub(crate) struct IqtreeGamma {
      pub(crate) shape: f64,
      #[deser(deserialize_as = ArrayVec)]
      pub(crate) category_rates: Array1<f64>,
      #[deser(deserialize_as = ArrayVec)]
      pub(crate) site_log_lh: Array1<f64>,
    }

    #[derive(Clone, Copy, Debug)]
    pub(crate) enum Representation {
      Dense,
      Sparse,
    }

    pub(crate) struct TreetimeCase {
      graph: Graph,
      branch_lengths: EdgeMap<f64>,
      dense: PartitionMarginalDense,
      sparse: PartitionMarginalSparse,
      sparse_nodes: NodeMap<SparseNodeState>,
    }

    impl TreetimeCase {
      pub(crate) fn new(inputs: &IqtreeInputs) -> Result<Self, Report> {
        let parsed = nwk_read(inputs.tree.as_bytes())?;
        let graph = parsed.graph;
        let branch_lengths = branch_lengths_or_zero(&parsed.branch_lengths);
        let records = inputs
          .alignment
          .iter()
          .map(|record| {
            Ok(AlignmentRecord {
              name: record.name.clone(),
              seq: Seq::try_from_str(&record.seq)?,
            })
          })
          .collect::<Result<Vec<_>, Report>>()?;
        let nodes = leaf_seq_inputs(&graph, &parsed.names, records);
        let dense = PartitionMarginalDense::new(Alphabet::default(), &graph, &nodes)?;
        let (sparse, sparse_nodes) =
          create_fitch_partition(&graph, Alphabet::default(), nodes)?.into_marginal_sparse(&graph)?;
        Ok(Self {
          graph,
          branch_lengths,
          dense,
          sparse: sparse.with_site_scales(true),
          sparse_nodes,
        })
      }

      pub(crate) fn site_log_lh(
        &self,
        model: &PartitionModel,
        representation: Representation,
      ) -> Result<Array1<f64>, Report> {
        let log_weights = model.gamma().weights().mapv(f64::ln);
        let mix = |values: &mut dyn Iterator<Item = f64>| mix_site(&log_weights.view(), values).0;
        match representation {
          Representation::Dense => {
            let (_, state) = dense_mixture_update(&self.dense, model, &self.graph, &self.branch_lengths)?;
            Ok(
              state
                .mixture
                .site_log_lh
                .rows()
                .into_iter()
                .map(|row| mix(&mut row.iter().copied()))
                .collect(),
            )
          },
          Representation::Sparse => {
            let (update, state) = sparse_mixture_update(
              &self.sparse,
              model,
              &self.graph,
              &self.branch_lengths,
              &self.sparse_nodes,
            )?;
            let root_sequence = &update.node_states[&self.graph.root_key()?].sequence;
            (0..self.sparse.length)
              .map(|pos| {
                state
                  .mixture
                  .site(pos, root_sequence[pos])
                  .map(|site| mix(&mut site.site_log_lh.iter().copied()))
                  .ok_or_else(|| make_internal_report!("Sparse mixture has no value at site {pos}"))
              })
              .collect()
          },
        }
      }
    }
  }
}
