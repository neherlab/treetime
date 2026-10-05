#[cfg(test)]
mod tests {
  use crate::graphviz::graphviz_write;
  use crate::nwk::nwk_read;
  use eyre::Report;
  use indoc::indoc;
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

  #[test]
  fn test_graphviz_writes_the_root_only_in_the_roots_subgraph() -> Result<(), Report> {
    let parse = nwk_read(b"((A:0.1,B:0.2)AB:0.3,C:0.4)root;".as_slice())?;
    let names = parse.names();

    let mut buf = Vec::new();
    graphviz_write(&mut buf, &parse.graph, &names, &parse.branch_lengths)?;

    let expected = indoc! {r#"
      digraph Phylogeny {
        graph [rankdir=LR, overlap=scale, splines=ortho, nodesep=1.0, ordering=out];
        edge  [overlap=scale];
        node  [shape=box];

        subgraph roots {
          4 [label="(4) root"]

          // fake edges for alignment of nodes
          {
            rank=same
            4-> 4 [style=invis, weight=1000]
          }
        }

        subgraph internals {
          2 [label="(2) AB"]
        }

        subgraph leaves {
          0 [label="(0) A"]
          1 [label="(1) B"]
          3 [label="(3) C"]

          // fake edges for alignment of nodes
          {
            rank=same
            0-> 0 [style=invis, weight=1000]
            0-> 1 [style=invis, weight=1100]
            0-> 3 [style=invis, weight=1200]
            1-> 0 [style=invis, weight=1300]
            1-> 1 [style=invis, weight=1400]
            1-> 3 [style=invis, weight=1500]
            3-> 0 [style=invis, weight=1600]
            3-> 1 [style=invis, weight=1700]
            3-> 3 [style=invis, weight=1800]
          }
        }

        2 -> 0 [xlabel="0.1", weight="0.1"]
        2 -> 1 [xlabel="0.2", weight="0.2"]
        4 -> 2 [xlabel="0.3", weight="0.3"]
        4 -> 3 [xlabel="0.4", weight="0.4"]
      }
    "#};
    assert_eq!(expected, String::from_utf8(buf)?);
    Ok(())
  }

  #[test]
  fn test_graphviz_escapes_quotes_and_backslashes_in_labels() -> Result<(), Report> {
    let parse = nwk_read(br#"('A\B':0.1,'C""D':0.2)root;"#.as_slice())?;
    let names = parse.names();

    let mut buf = Vec::new();
    graphviz_write(&mut buf, &parse.graph, &names, &parse.branch_lengths)?;
    let output = String::from_utf8(buf)?;

    let labels: Vec<&str> = output.lines().filter(|line| line.contains("[label=")).collect();
    let expected = vec![
      r#"    2 [label="(2) root"]"#,
      r#"    0 [label="(0) A\\B"]"#,
      r#"    1 [label="(1) C\"\"D"]"#,
    ];
    assert_eq!(expected, labels);
    Ok(())
  }
}
