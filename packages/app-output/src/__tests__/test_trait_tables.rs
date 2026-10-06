#[cfg(test)]
mod tests {
  use crate::__tests__::test_tree_output::tests::helpers::{mugration_graph, mugration_setup, node_key, topology_from};
  use crate::annotated_graph::AnnotatedTreeView;
  use crate::augur_node_data_traits::build_augur_node_data_traits;
  use crate::trait_tables::{write_trait_confidence_rows, write_traits_rows};
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::gtr::gtr::GTR;
  use treetime::partition::storage::discrete::DiscreteStates;
  use treetime_io::csv::{CsvWriter, TableFormat};
  use treetime_utils::{assert_error, o};

  #[rustfmt::skip]
  #[rstest]
  #[case::plain(            ("country",       "A",     "usa"),               "node,country\nA,usa\n")]
  #[case::comma_in_value(   ("country",       "A",     "Congo, DR"),         "node,country\nA,\"Congo, DR\"\n")]
  #[case::quote_in_value(   ("country",       "A",     "Cote \"d'Ivoire\""), "node,country\nA,\"Cote \"\"d'Ivoire\"\"\"\n")]
  #[case::newline_in_value( ("country",       "A",     "a\nb"),              "node,country\nA,\"a\nb\"\n")]
  #[case::comma_in_node(    ("country",       "A/1,2", "usa"),               "node,country\n\"A/1,2\",usa\n")]
  #[case::comma_in_header(  ("region, state", "A",     "usa"),               "node,\"region, state\"\nA,usa\n")]
  #[trace]
  fn test_trait_tables_traits_csv_quotes_fields_per_rfc4180(
    #[case] (attribute, node, value): (&str, &str, &str),
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let mut setup = mugration_setup()?;
    setup.topology = topology_from("(A:0.1)root;")?;
    let leaf = node_key(&setup.topology, "A");
    let root = node_key(&setup.topology, "root");
    setup.topology.names.insert(leaf, Some(o!(node)));
    setup.values = btreemap! { leaf => Some(o!(value)), root => None };
    let graph = mugration_graph(&setup, attribute);

    let mut buf = Vec::new();
    let mut csv = CsvWriter::new(&mut buf, TableFormat::Csv);
    write_traits_rows(&AnnotatedTreeView::new(&graph)?, &mut csv)?;
    csv.into_inner()?;

    assert_eq!(expected, String::from_utf8(buf)?);
    Ok(())
  }

  #[test]
  fn test_trait_tables_traits_csv_writes_one_row_per_node_with_a_repeated_name() -> Result<(), Report> {
    let mut setup = mugration_setup()?;
    setup.topology = topology_from("(A:0.1,B:0.1)root;")?;
    let (a, b, root) = (
      node_key(&setup.topology, "A"),
      node_key(&setup.topology, "B"),
      node_key(&setup.topology, "root"),
    );
    setup.topology.names.insert(b, Some(o!("A")));
    setup.values = btreemap! { a => Some(o!("usa")), b => Some(o!("uk")), root => None };
    let graph = mugration_graph(&setup, "country");

    let mut buf = Vec::new();
    let mut csv = CsvWriter::new(&mut buf, TableFormat::Csv);
    write_traits_rows(&AnnotatedTreeView::new(&graph)?, &mut csv)?;
    csv.into_inner()?;

    assert_eq!("node,country\nA,usa\nA,uk\n", String::from_utf8(buf)?);
    Ok(())
  }

  #[test]
  fn test_trait_tables_confidence_csv_quotes_fields_per_rfc4180() -> Result<(), Report> {
    let mut setup = mugration_setup()?;
    setup.topology = topology_from("(A:0.1,B:0.1)root;")?;
    let (a, b, root) = (
      node_key(&setup.topology, "A"),
      node_key(&setup.topology, "B"),
      node_key(&setup.topology, "root"),
    );
    setup.topology.names.insert(a, Some(o!("A,1")));
    setup.states = DiscreteStates::from_values(["Congo, DR", "usa"].into_iter(), "?");
    setup.values = btreemap! { a => Some(o!("usa")), b => Some(o!("Congo, DR")), root => None };
    setup.profiles = btreemap! { a => Some(array![0.25, 0.75]), b => Some(array![1.0, 0.0]), root => None };
    let graph = mugration_graph(&setup, "country");

    let mut buf = Vec::new();
    let mut csv = CsvWriter::new(&mut buf, TableFormat::Csv);
    write_trait_confidence_rows(&AnnotatedTreeView::new(&graph)?, &mut csv)?;
    csv.into_inner()?;

    let expected = "\
      node,\"Congo, DR\",usa\n\
      \"A,1\",0.250000,0.750000\n\
      B,1.000000,0.000000\n\
    ";
    assert_eq!(expected, String::from_utf8(buf)?);
    Ok(())
  }

  #[test]
  fn test_trait_tables_node_data_rejects_non_finite_probability() -> Result<(), Report> {
    let mut setup = mugration_setup()?;
    let a = node_key(&setup.topology, "A");
    setup.profiles.insert(a, Some(array![f64::NAN, 1.0]));
    let graph = mugration_graph(&setup, "country");
    let model = GTR::builder().n_states(2).mu(1.0).pi(array![0.5, 0.5]).build()?;

    assert_error!(
      build_augur_node_data_traits(&AnnotatedTreeView::new(&graph)?, &model),
      "Node 'A' has non-finite trait state 'CH' probability=NaN"
    );
    Ok(())
  }
}
