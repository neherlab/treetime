#[cfg(test)]
mod tests {
  use crate::results::auspice::display_auspice;
  use crate::results::tree::result_colorings;
  use eyre::Report;
  use helpers::{MUTED, auspice, scales};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use treetime::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::trait_with_two_states(  ("country",    "categorical", 2),  MUTED[..2].to_vec())]
  #[case::trait_with_nine_states( ("country",    "categorical", 9),  MUTED.to_vec())]
  #[case::trait_with_ten_states(  ("country",    "categorical", 10), vec![])]
  #[case::genotype(               ("gt",         "categorical", 2),  vec![])]
  #[case::clock_outliers(         ("bad_branch", "categorical", 2),  vec![])]
  #[case::continuous(             ("rate",       "continuous",  2),  vec![])]
  #[trace]
  fn test_display_auspice_sets_the_app_palette_on_small_categorical_traits(
    #[case] (key, kind, states): (&str, &str, usize),
    #[case] expected: Vec<&str>,
  ) -> Result<(), Report> {
    let document = display_auspice(auspice(key, kind, states, None))?;

    let colors = document.0.data.meta.colorings[0]
      .scale
      .iter()
      .map(|[_, color]| color.as_str())
      .collect::<Vec<_>>();
    assert_eq!(expected, colors);
    Ok(())
  }

  #[test]
  fn test_display_auspice_scale_pairs_the_sorted_states_with_the_palette() -> Result<(), Report> {
    let document = display_auspice(auspice("country", "categorical", 3, None))?;

    assert_eq!(
      vec![
        [o!("s0"), o!("#332288")],
        [o!("s1"), o!("#88ccee")],
        [o!("s2"), o!("#44aa99")],
      ],
      document.0.data.meta.colorings[0].scale
    );
    Ok(())
  }

  #[test]
  fn test_display_auspice_keeps_a_scale_the_file_sets_on_a_coloring_it_leaves_alone() -> Result<(), Report> {
    let own = json!([["s0", "#000000"], ["s1", "#ffffff"]]);

    let document = display_auspice(auspice("gt", "categorical", 2, Some(own)))?;

    assert_eq!(
      vec![[o!("s0"), o!("#000000")], [o!("s1"), o!("#ffffff")]],
      document.0.data.meta.colorings[0].scale
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::trait_with_two_states(("country",    "categorical", 2,  None))]
  #[case::trait_with_ten_states(("country",    "categorical", 10, None))]
  #[case::genotype_with_scale(  ("gt",         "categorical", 2,  Some(json!([["s0", "#000000"]]))))]
  #[case::clock_outliers(       ("bad_branch", "categorical", 2,  None))]
  #[trace]
  fn test_display_auspice_scales_match_the_result_colorings(
    #[case] (key, kind, states, own): (&str, &str, usize, Option<Value>),
  ) -> Result<(), Report> {
    let input = auspice(key, kind, states, own);
    let colorings = result_colorings(&input)?;

    let document = display_auspice(input)?;

    assert_eq!(
      colorings
        .into_iter()
        .map(|coloring| coloring.scale)
        .collect::<Vec<_>>(),
      scales(&document.0)
    );
    Ok(())
  }

  #[test]
  fn test_display_auspice_changes_nothing_but_the_color_scales() -> Result<(), Report> {
    let input = auspice("country", "categorical", 2, None);
    let mut expected = serde_json::to_value(&input)?;
    expected["meta"]["colorings"][0]["scale"] = json!([["s0", "#332288"], ["s1", "#88ccee"]]);

    let document = display_auspice(input)?;

    assert_eq!(expected, serde_json::to_value(&document)?);
    Ok(())
  }

  mod helpers {
    use crate::results::tree::StateColor;
    use serde_json::{Value, json};
    use treetime_io::auspice_types::AuspiceTree;

    pub(super) const MUTED: [&str; 9] = [
      "#332288", "#88ccee", "#44aa99", "#117733", "#999933", "#ddcc77", "#cc6677", "#882255", "#aa4499",
    ];

    pub(super) fn auspice(key: &str, kind: &str, states: usize, scale: Option<Value>) -> AuspiceTree {
      let children = (0..states)
        .map(|index| json!({ "name": format!("tip{index}"), "node_attrs": { key: { "value": format!("s{index}") } } }))
        .collect::<Vec<_>>();
      let mut coloring = json!({ "key": key, "title": key, "type": kind });
      if let Some(scale) = scale {
        coloring["scale"] = scale;
      }
      serde_json::from_value(json!({
        "version": "v2",
        "meta": { "colorings": [coloring], "display_defaults": { "color_by": key } },
        "tree": { "name": "root", "node_attrs": { "div": 0.0 }, "children": children }
      }))
      .unwrap()
    }

    pub(super) fn scales(tree: &AuspiceTree) -> Vec<Vec<StateColor>> {
      tree
        .data
        .meta
        .colorings
        .iter()
        .map(|coloring| {
          coloring
            .scale
            .iter()
            .map(|[state, color]| StateColor {
              state: state.clone(),
              color: color.clone(),
            })
            .collect()
        })
        .collect()
    }
  }
}
