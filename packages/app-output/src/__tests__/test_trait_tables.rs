#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{
    mugration_graph, mugration_setup, mugration_setup_from, node_key, rename_node, topology_from,
  };
  use crate::annotated_graph::AnnotatedTreeView;
  use crate::augur_node_data_traits::build_augur_node_data_traits;
  use crate::trait_tables::{write_trait_confidence_rows, write_traits_rows};
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::gtr::gtr::{GTR, Mu};
  use treetime::partition::marginal::discrete::state_posterior::StatePosterior;
  use treetime_io::csv::{CsvWriter, TableFormat};
  use treetime_utils::assert_error;

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
    let mut topology = topology_from("A;")?;
    let leaf = node_key(&topology, "A");
    rename_node(&mut topology, leaf, node)?;
    let setup = mugration_setup_from(topology, &[value, "other"], &btreemap! { leaf => value })?;
    let graph = mugration_graph(&setup, attribute);

    let mut buf = Vec::new();
    let mut csv = CsvWriter::new(&mut buf, TableFormat::Csv);
    write_traits_rows(&AnnotatedTreeView::new(&graph)?, &mut csv)?;
    csv.into_inner()?;

    assert_eq!(expected, String::from_utf8(buf)?);
    Ok(())
  }

  #[test]
  fn test_trait_tables_traits_csv_writes_unique_names_for_a_repeated_name() -> Result<(), Report> {
    let mut topology = topology_from("(A:0.1)B;")?;
    let (a, b) = (node_key(&topology, "A"), node_key(&topology, "B"));
    rename_node(&mut topology, b, "A")?;
    let setup = mugration_setup_from(topology, &["uk", "usa"], &btreemap! { a => "usa", b => "uk" })?;
    let graph = mugration_graph(&setup, "country");

    let mut buf = Vec::new();
    let mut csv = CsvWriter::new(&mut buf, TableFormat::Csv);
    write_traits_rows(&AnnotatedTreeView::new(&graph)?, &mut csv)?;
    csv.into_inner()?;

    assert_eq!("node,country\nA,usa\nA_2,uk\n", String::from_utf8(buf)?);
    Ok(())
  }

  #[test]
  fn test_trait_tables_confidence_csv_quotes_fields_per_rfc4180() -> Result<(), Report> {
    let mut topology = topology_from("(A:0.1)B;")?;
    let (a, b) = (node_key(&topology, "A"), node_key(&topology, "B"));
    rename_node(&mut topology, a, "A,1")?;
    let mut setup = mugration_setup_from(
      topology,
      &["Congo, DR", "usa"],
      &btreemap! { a => "usa", b => "Congo, DR" },
    )?;
    setup.posteriors = btreemap! {
      a => StatePosterior::new(array![0.25, 0.75])?,
      b => StatePosterior::new(array![1.0, 0.0])?,
    };
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
    setup.posteriors.insert(a, StatePosterior::new(array![f64::NAN, 1.0])?);
    let graph = mugration_graph(&setup, "country");
    let model = GTR::builder()
      .n_states(2)
      .mu(Mu::FromMatrix(1.0))
      .pi(array![0.5, 0.5])
      .build()?;

    assert_error!(
      build_augur_node_data_traits(&AnnotatedTreeView::new(&graph)?, &model),
      "Node 'A' has non-finite trait state 'CH' probability=NaN"
    );
    Ok(())
  }
}
