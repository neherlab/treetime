#[cfg(test)]
mod tests {
  use crate::reroot::__tests__::test_midpoint::tests::helpers::midpoint;
  use eyre::Report;
  use helpers::{dataset_tree, load_gm_reroot_midpoint_outputs, root_sides};
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::pretty_assert_relative_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_20(        "zika/20")]
  #[case::ebola_20(       "ebola/20")]
  #[case::lassa_l_20(     "lassa/L/20")]
  #[case::tb_20(          "tb/20")]
  #[case::mpox_clade_2_20("mpox/clade-ii/20")]
  #[trace]
  fn test_gm_reroot_midpoint_matches_biopython_root_at_midpoint(#[case] dataset: &str) -> Result<(), Report> {
    let outputs = load_gm_reroot_midpoint_outputs()?;
    let expected = &outputs.datasets[dataset];

    let actual = midpoint(&dataset_tree(dataset)?, 1)?;

    assert_eq!(expected.sides, root_sides(&actual.newick)?);
    assert_eq!(
      expected.root_to_tip.keys().collect::<Vec<_>>(),
      actual.root_to_tip.keys().collect::<Vec<_>>()
    );
    pretty_assert_relative_eq!(
      Array1::from_iter(expected.root_to_tip.values().copied()),
      Array1::from_iter(actual.root_to_tip.values().copied()),
      max_relative = 1e-12
    );
    Ok(())
  }

  mod helpers {
    use deser::Deserialize;
    use eyre::Report;
    use itertools::Itertools;
    use std::collections::BTreeMap;
    use std::fs::read_to_string;
    use treetime_graph::clades::{TipIds, graph_clade_keys};
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::nwk_read;
    use treetime_utils::io::json::json_read_str;

    #[derive(Debug, Deserialize)]
    pub(super) struct GmOutputs {
      pub(super) datasets: BTreeMap<String, GmMidpointRoot>,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct GmMidpointRoot {
      pub(super) sides: Vec<Vec<String>>,
      pub(super) root_to_tip: BTreeMap<String, f64>,
    }

    pub(super) fn load_gm_reroot_midpoint_outputs() -> Result<GmOutputs, Report> {
      let path = format!(
        "{}/src/reroot/__tests__/__fixtures__/gm_reroot_midpoint_outputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      json_read_str(&read_to_string(path)?)
    }

    pub(super) fn dataset_tree(dataset: &str) -> Result<String, Report> {
      Ok(read_to_string(format!(
        "{}/../../data/{dataset}/tree.nwk",
        env!("CARGO_MANIFEST_DIR")
      ))?)
    }

    pub(super) fn root_sides(newick: &str) -> Result<Vec<Vec<String>>, Report> {
      let parsed = nwk_read(newick.as_bytes())?;
      let mut ids = TipIds::default();
      let clades = graph_clade_keys(&parsed.graph, &parsed.names, &mut ids)?;
      let schedule = GraphSchedule::new(&parsed.graph)?;
      let root_key = parsed.graph.get_exactly_one_root()?.key();
      Ok(
        schedule
          .children(root_key)
          .iter()
          .map(|(child_key, _)| {
            clades[child_key]
              .iter()
              .map(|&id| ids.name(id).to_owned())
              .sorted()
              .collect_vec()
          })
          .sorted()
          .collect(),
      )
    }
  }
}
