#[cfg(test)]
pub(crate) mod tests {
  pub(crate) mod helpers {
    use itertools::Itertools;
    use treetime_graph::keymap;

    use statrs::stats_tests::chisquare::chisquare;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::node_names::NodeNames;
    use treetime_utils::cast::ToF64Exact;

    pub(crate) struct ThreeTipTree {
      pub(crate) graph: Graph,
      pub(crate) names: NodeNames,
      pub(crate) root: GraphNodeKey,
      pub(crate) inner: GraphNodeKey,
      pub(crate) tip_a: GraphNodeKey,
      pub(crate) tip_b: GraphNodeKey,
      pub(crate) tip_c: GraphNodeKey,
      pub(crate) root_inner: GraphEdgeKey,
      pub(crate) root_c: GraphEdgeKey,
      pub(crate) inner_a: GraphEdgeKey,
      pub(crate) inner_b: GraphEdgeKey,
    }

    pub(crate) fn three_tip_tree() -> ThreeTipTree {
      let mut graph = Graph::new();
      let root = graph.add_node();
      let inner = graph.add_node();
      let tip_a = graph.add_node();
      let tip_b = graph.add_node();
      let tip_c = graph.add_node();
      let root_inner = graph.add_edge(root, inner).unwrap();
      let root_c = graph.add_edge(root, tip_c).unwrap();
      let inner_a = graph.add_edge(inner, tip_a).unwrap();
      let inner_b = graph.add_edge(inner, tip_b).unwrap();
      graph.build().unwrap();
      let input = keymap! {
        root => Some("root".to_owned()),
        inner => Some("inner".to_owned()),
        tip_a => Some("A".to_owned()),
        tip_b => Some("B".to_owned()),
        tip_c => Some("C".to_owned()),
      };
      let names = NodeNames::from_input(input, &graph).unwrap().names;
      ThreeTipTree {
        graph,
        names,
        root,
        inner,
        tip_a,
        tip_b,
        tip_c,
        root_inner,
        root_c,
        inner_a,
        inner_b,
      }
    }

    pub(crate) fn chisquare_p_value(observed: &[usize], probabilities: &[f64]) -> f64 {
      let probability_sum: f64 = probabilities.iter().sum();
      assert!(
        (probability_sum - 1.0).abs() < 1e-12,
        "bin probabilities sum to {probability_sum}"
      );
      let total = observed.iter().sum::<usize>().to_f64_exact();
      let expected = probabilities
        .iter()
        .map(|probability| total * probability / probability_sum)
        .collect_vec();

      let (mut merged, mut kept): (Vec<usize>, Vec<usize>) = (0..observed.len()).partition(|&bin| expected[bin] < 5.0);
      kept.sort_by(|&left, &right| expected[left].total_cmp(&expected[right]));
      while !merged.is_empty() && merged.iter().map(|&bin| expected[bin]).sum::<f64>() < 5.0 {
        merged.push(kept.remove(0));
      }

      let mut observed_bins = kept.iter().map(|&bin| observed[bin]).collect_vec();
      let mut expected_bins = kept.iter().map(|&bin| expected[bin]).collect_vec();
      if !merged.is_empty() {
        observed_bins.push(merged.iter().map(|&bin| observed[bin]).sum());
        expected_bins.push(merged.iter().map(|&bin| expected[bin]).sum());
      }
      let (_, p_value) = chisquare(&observed_bins, Some(&expected_bins), None).unwrap();
      p_value
    }
  }
}
