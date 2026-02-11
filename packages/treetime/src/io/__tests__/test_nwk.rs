#[cfg(test)]
mod tests {
  use crate::graph::graph_tests::tests::{TestEdge, TestNode};
  use crate::io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
  use eyre::Report;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_nwk_roundtrip() -> Result<(), Report> {
    let input = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root;";
    let graph = nwk_read_str::<TestNode, TestEdge, ()>(input)?;
    let output = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(input, output);
    Ok(())
  }
}
