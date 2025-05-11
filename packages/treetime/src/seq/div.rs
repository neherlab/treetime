use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use maplit::btreemap;
use std::collections::BTreeMap;

#[derive(Debug, Default, Copy, Clone)]
pub struct OnlyLeaves(pub bool);

/// Calculate mapping of node name to node divergence (accumulated by summing branch lengths)
pub fn calculate_divs<N: GraphNode + Named, E: GraphEdge + Weighted, D: Send + Sync>(
  graph: &Graph<N, E, D>,
  only_leaves: OnlyLeaves,
) -> BTreeMap<String, f64> {
  let mut divs: BTreeMap<String, f64> = btreemap! {};
  let mut result: BTreeMap<String, f64> = btreemap! {};
  graph.iter_depth_first_preorder_forward(|node| {
    let name = node.payload.name().unwrap().as_ref().to_owned();
    if node.is_root {
      divs.insert(name.clone(), 0.0);
      if !only_leaves.0 {
        result.insert(name, 0.0);
      }
    } else {
      let (parent, edge) = node.get_exactly_one_parent().unwrap();
      let parent_div: f64 = {
        let parent_name = parent.read_arc().name().unwrap().as_ref().to_owned();
        divs.get(&parent_name).copied().unwrap_or_default()
      };
      let branch_length = edge.read_arc().weight().unwrap_or_default();
      let div = parent_div + branch_length;
      divs.insert(name.clone(), div);
      if node.is_leaf || !only_leaves.0 {
        result.insert(name, div);
      }
    }
  });
  result
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::commands::clock::clock_graph::ClockGraph;
  use crate::graph::graph::tests::{TestEdge, TestNode};
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_calculate_divs() -> Result<(), Report> {
    let graph: Graph<TestNode, TestEdge, ()> = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let actual = calculate_divs(&graph, OnlyLeaves(false));
    let expected = btreemap! {
      o!("A") => 0.20000000298023224,
      o!("AB") => 0.10000000149011612,
      o!("B") => 0.30000000447034836,
      o!("C") => 0.2500000037252903,
      o!("CD") => 0.05000000074505806,
      o!("D") => 0.16999999806284904,
      o!("root") => 0.0,
    };
    assert_eq!(&expected, &actual);
    Ok(())
  }

  #[test]
  fn test_calculate_divs_only_leaves() -> Result<(), Report> {
    let graph: ClockGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let actual = calculate_divs(&graph, OnlyLeaves(true));
    let expected = btreemap! {
      o!("A") => 0.20000000298023224,
      o!("B") => 0.30000000447034836,
      o!("C") => 0.2500000037252903,
      o!("D") => 0.16999999806284904,
    };
    assert_eq!(&expected, &actual);
    Ok(())
  }
}
