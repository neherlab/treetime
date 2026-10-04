#[cfg(test)]
mod tests {
  use crate::graphviz::graphviz_write;
  use crate::nwk::nwk_read;
  use eyre::Report;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_graphviz_output_starts_with_the_graph_and_ends_with_one_newline() -> Result<(), Report> {
    let parse = nwk_read(b"(A:0.1,B:0.2)root;".as_slice())?;
    let names = parse.names();

    let mut buf = Vec::new();
    graphviz_write(&mut buf, &parse.graph, &names, &parse.branch_lengths)?;
    let output = String::from_utf8(buf)?;

    assert_eq!(
      (Some("digraph Phylogeny {"), true, false),
      (output.lines().next(), output.ends_with("}\n"), output.ends_with("\n\n"))
    );
    Ok(())
  }
}
