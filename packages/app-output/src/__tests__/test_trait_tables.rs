#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{
    ambiguity_for, mugration_graph, mugration_graph_with_ambiguity, mugration_setup, mugration_setup_from, node_key,
    rename_node, topology_from,
  };
  use crate::annotated_graph::AnnotatedTreeView;
  use crate::augur_node_data_traits::build_augur_node_data_traits;
  use crate::output_plan::OutputSelection;
  use crate::trait_tables::{write_trait_confidence_rows, write_traits_rows};
  use deser_value::Value;
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
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
      \"A,1\",0.25,0.75\n\
      B,1.0,0.0\n\
    ";
    assert_eq!(expected, String::from_utf8(buf)?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::carried(    &[OutputSelection::Traits],     "node,country,country_ambiguous\nA,Samoa,\nB,Samoa,American Samoa|Samoa\n")]
  #[case::not_carried(&[OutputSelection::NodeData], "node,country\nA,Samoa\nB,Samoa\n")]
  #[trace]
  fn test_trait_tables_traits_csv_writes_ambiguous_states_when_it_carries_the_report(
    #[case] carriers: &[OutputSelection],
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let setup = helpers::samoa_setup()?;
    let ambiguity = ambiguity_for(carriers, 0.125)?;
    let graph = mugration_graph_with_ambiguity(&setup, "country", &ambiguity);

    let mut buf = Vec::new();
    let mut csv = CsvWriter::new(&mut buf, TableFormat::Csv);
    write_traits_rows(&AnnotatedTreeView::new(&graph)?, &mut csv)?;
    csv.into_inner()?;

    assert_eq!(expected, String::from_utf8(buf)?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::carried(    &[OutputSelection::NodeData], Some("American Samoa|Samoa"))]
  #[case::not_carried(&[OutputSelection::Traits],     None)]
  #[trace]
  fn test_trait_tables_node_data_writes_ambiguous_states_when_it_carries_the_report(
    #[case] carriers: &[OutputSelection],
    #[case] expected: Option<&str>,
  ) -> Result<(), Report> {
    let setup = helpers::samoa_setup()?;
    let ambiguity = ambiguity_for(carriers, 0.125)?;
    let graph = mugration_graph_with_ambiguity(&setup, "country", &ambiguity);

    let node_data = build_augur_node_data_traits(&AnnotatedTreeView::new(&graph)?, &helpers::model(3)?)?;

    let ambiguous = |name: &str| node_data.nodes[name].fields.get("country_ambiguous").cloned();
    assert_eq!(None, ambiguous("A"));
    assert_eq!(expected.map(Value::from), ambiguous("B"));
    Ok(())
  }

  #[test]
  fn test_trait_tables_node_data_rejects_non_finite_probability() -> Result<(), Report> {
    let mut setup = mugration_setup()?;
    let a = node_key(&setup.topology, "A");
    setup.posteriors.insert(a, StatePosterior::new(array![f64::NAN, 1.0])?);
    let graph = mugration_graph(&setup, "country");
    let model = helpers::model(2)?;

    assert_error!(
      build_augur_node_data_traits(&AnnotatedTreeView::new(&graph)?, &model),
      "Node 'A' has non-finite trait state 'CH' probability=NaN"
    );
    Ok(())
  }

  mod helpers {
    use crate::__tests__::test_support::tests::{MugrationSetup, mugration_setup_from, node_key, topology_from};
    use eyre::Report;
    use maplit::btreemap;
    use ndarray::{Array1, array};
    use treetime::gtr::gtr::{GTR, Mu};
    use treetime::partition::marginal::discrete::state_posterior::StatePosterior;
    use treetime_utils::cast::ToF64Exact;

    pub(super) fn samoa_setup() -> Result<MugrationSetup, Report> {
      let topology = topology_from("(A:0.1)B;")?;
      let (a, b) = (node_key(&topology, "A"), node_key(&topology, "B"));
      let mut setup = mugration_setup_from(
        topology,
        &["American Samoa", "Samoa", "Tonga"],
        &btreemap! { a => "Samoa", b => "Samoa" },
      )?;
      setup.posteriors = btreemap! {
        a => StatePosterior::new(array![0.0, 1.0, 0.0])?,
        b => StatePosterior::new(array![0.4, 0.45, 0.15])?,
      };
      Ok(setup)
    }

    pub(super) fn model(n_states: usize) -> Result<GTR, Report> {
      GTR::builder()
        .n_states(n_states)
        .mu(Mu::FromMatrix(1.0))
        .pi(Array1::from_elem(n_states, 1.0 / n_states.to_f64_exact()))
        .build()
    }
  }
}
