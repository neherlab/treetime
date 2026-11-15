#[cfg(test)]
mod tests {
  use crate::nwk::TREE_DIALECT_DEFAULT;
  use crate::tree::segment_trees_read_files;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeSet;
  use std::fs;
  use std::path::PathBuf;
  use treetime_graph::schedule::GraphSchedule;
  use treetime_utils::make_report;

  const SEGMENTS: [&str; 3] = ["seg1", "seg2", "seg3"];

  fn toy_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/toy/reassortment")
  }

  fn expected_edges() -> Result<BTreeSet<(String, String, String, String)>, Report> {
    let table = fs::read_to_string(toy_dir().join("graph.tsv"))?;
    table
      .lines()
      .skip(1)
      .map(|line| {
        let (parent, child, length, segments) = line
          .split('\t')
          .collect_tuple()
          .ok_or_else(|| make_report!("graph.tsv: malformed line '{line}'"))?;
        Ok((
          parent.to_owned(),
          child.to_owned(),
          length.to_owned(),
          segments.to_owned(),
        ))
      })
      .collect()
  }

  #[test]
  fn test_segment_trees_read_toy_reassortment_gives_the_edges_and_scopes_of_its_edge_table() -> Result<(), Report> {
    let paths = SEGMENTS.map(|segment| toy_dir().join(format!("tree_{segment}.nwk")));
    let merged = segment_trees_read_files(&paths, TREE_DIALECT_DEFAULT)?;

    let actual: BTreeSet<_> = merged
      .graph
      .get_edges()
      .map(|edge| {
        let segments = SEGMENTS
          .iter()
          .enumerate()
          .filter(|&(index, _)| merged.scopes.carries(edge.key(), index))
          .map(|(_, name)| *name)
          .join(",");
        let length = merged.branch_lengths[&edge.key()].unwrap_or(f64::NAN);
        (
          merged.names[edge.source()].to_owned(),
          merged.names[edge.target()].to_owned(),
          format!("{length:.6}"),
          segments,
        )
      })
      .collect();
    assert_eq!(expected_edges()?, actual);

    let schedule = GraphSchedule::new(&merged.graph)?;
    let roots = schedule.roots().iter().map(|&key| &merged.names[key]).collect_vec();
    assert_eq!(vec!["R", "R3"], roots);
    Ok(())
  }
}
