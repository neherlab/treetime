#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{c, edge_mutation_map};
  use app_output::DateCommentProvider;
  use app_output::EdgeMutationCommentProvider;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime::seq::mutation::Sub;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nex::{NexWriteOptions, nex_write_str_with};
  use treetime_io::nwk::{CommentProviders, NodeCommentProvider, NwkStyle, nwk_read_str};

  #[test]
  fn test_timetree_mutation_provider_produces_comments() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let edge_subs = &[(
      0,
      vec![
        Sub::new(c(b'A'), 54_usize, c(b'G'))?,
        Sub::new(c(b'T'), 92_usize, c(b'C'))?,
      ],
    )];
    let edge_mutations = edge_mutation_map(&graph, edge_subs, &[]);
    let provider = EdgeMutationCommentProvider::new(&edge_mutations, &graph);
    let leaf_key = graph.get_leaves().collect::<Vec<_>>()[0].key();
    let comments = provider.node_comments(leaf_key)?;
    assert_eq!(comments.get("mutations").map(String::as_str), Some("A55G,T93C"));
    Ok(())
  }

  #[test]
  fn test_timetree_nexus_output_includes_mutations_and_date() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let edge_subs = &[(
      0,
      vec![
        Sub::new(c(b'A'), 54_usize, c(b'G'))?,
        Sub::new(c(b'T'), 92_usize, c(b'C'))?,
      ],
    )];

    let date_times: BTreeMap<GraphNodeKey, f64> = graph.get_leaves().map(|leaf| (leaf.key(), 2003.84)).collect();

    let edge_mutations = edge_mutation_map(&graph, edge_subs, &[]);
    let provider = EdgeMutationCommentProvider::new(&edge_mutations, &graph);
    let date_provider = DateCommentProvider::new(&date_times);
    let providers = CommentProviders::new().with(&provider).with(&date_provider);
    let options = NexWriteOptions {
      style: NwkStyle::Beast,
      ..NexWriteOptions::default()
    };
    let time_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = graph.get_edges().map(|edge| (edge.key(), None)).collect();
    let nexus = nex_write_str_with(&graph, &names, &time_lengths, &options, &providers)?;
    let expected = concat!(
      indoc! {r#"
        #NEXUS
        Begin Taxa;
          Dimensions NTax=1;
          TaxLabels A;
        End;
        Begin Trees;
          Tree tree1=(A[&date=2003.84,mutations="A55G,T93C"])root;
        End;
      "#},
      "\n"
    );
    assert_eq!(nexus, expected);
    Ok(())
  }
}
