#[cfg(test)]
pub(crate) mod tests {
  use crate::results::mugration::{AncestorState, MugrationResults, StateChange, mugration_results};
  use crate::results::mutations::{AncestralResults, BranchMutations, RecurrentSite, ancestral_results};
  use crate::results::tree::{DateInterval, ResultColoring, ResultNode, ResultTree, StateColor};
  use eyre::Report;
  use helpers::fixture;
  use pretty_assertions::assert_eq;
  use treetime::o;
  use treetime_utils::{assert_error, vec_of_owned};

  #[test]
  fn test_result_tree_lists_nodes_in_preorder_with_their_values() -> Result<(), Report> {
    let tree = ResultTree::from_auspice(&fixture())?;

    let expected = ResultTree {
      nodes: vec![
        ResultNode {
          name: o!("root"),
          parent: None,
          children: vec![1, 4],
          tips: 3,
          div: Some(0.0),
          date: Some(2010.0),
          date_interval: Some(DateInterval {
            lower: 2009.0,
            upper: 2011.0,
            days: 730.0,
          }),
          excluded: None,
          mutations: vec![],
        },
        ResultNode {
          name: o!("AB"),
          parent: Some(0),
          children: vec![2, 3],
          tips: 2,
          div: Some(0.01),
          date: Some(2012.0),
          date_interval: Some(DateInterval {
            lower: 2011.5,
            upper: 2012.5,
            days: 365.5,
          }),
          excluded: None,
          mutations: vec_of_owned!["A10G", "C20T"],
        },
        ResultNode {
          name: o!("A"),
          parent: Some(1),
          children: vec![],
          tips: 1,
          div: Some(0.02),
          date: Some(2015.0),
          date_interval: None,
          excluded: Some(false),
          mutations: vec_of_owned!["G10A"],
        },
        ResultNode {
          name: o!("B"),
          parent: Some(1),
          children: vec![],
          tips: 1,
          div: Some(0.015),
          date: Some(2014.0),
          date_interval: None,
          excluded: Some(false),
          mutations: vec![],
        },
        ResultNode {
          name: o!("C"),
          parent: Some(0),
          children: vec![],
          tips: 1,
          div: Some(0.03),
          date: Some(2016.0),
          date_interval: None,
          excluded: Some(true),
          mutations: vec_of_owned!["C20T", "G30-"],
        },
      ],
      colorings: vec![
        ResultColoring {
          key: o!("num_date"),
          title: o!("Date"),
          kind: o!("continuous"),
          states: vec![],
          scale: vec![],
        },
        ResultColoring {
          key: o!("country"),
          title: o!("Country"),
          kind: o!("categorical"),
          states: vec_of_owned!["brazil", "china"],
          scale: vec![
            StateColor {
              state: o!("brazil"),
              color: o!("#332288"),
            },
            StateColor {
              state: o!("china"),
              color: o!("#88ccee"),
            },
          ],
        },
        ResultColoring {
          key: o!("gt"),
          title: o!("Genotype"),
          kind: o!("categorical"),
          states: vec![],
          scale: vec![],
        },
      ],
      default_color_by: Some(o!("country")),
    };
    assert_eq!(expected, tree);
    Ok(())
  }

  #[test]
  fn test_result_tree_drops_an_empty_date_interval() -> Result<(), Report> {
    let mut auspice = fixture();
    auspice.tree.node_attrs.num_date.as_mut().unwrap().confidence = Some([2010.0, 2010.0]);

    let tree = ResultTree::from_auspice(&auspice)?;

    assert_eq!(None, tree.root().date_interval);
    Ok(())
  }

  #[test]
  fn test_ancestral_results_rank_branches_and_recurrent_sites() -> Result<(), Report> {
    let tree = ResultTree::from_auspice(&fixture())?;

    let expected = AncestralResults {
      mutations: 5,
      branches: vec![
        BranchMutations {
          name: o!("AB"),
          tips: 2,
          mutations: vec_of_owned!["A10G", "C20T"],
        },
        BranchMutations {
          name: o!("C"),
          tips: 1,
          mutations: vec_of_owned!["C20T", "G30-"],
        },
        BranchMutations {
          name: o!("A"),
          tips: 1,
          mutations: vec_of_owned!["G10A"],
        },
      ],
      recurrent_sites: vec![
        RecurrentSite {
          position: 10,
          branches: 2,
        },
        RecurrentSite {
          position: 20,
          branches: 2,
        },
      ],
    };
    assert_eq!(expected, ancestral_results(Some(&tree))?);
    Ok(())
  }

  #[test]
  fn test_ancestral_results_reject_a_mutation_with_two_positions() -> Result<(), Report> {
    let mut auspice = fixture();
    auspice
      .tree
      .branch_attrs
      .mutations
      .insert(o!("nuc"), vec_of_owned!["A1B2"]);
    let tree = ResultTree::from_auspice(&auspice)?;

    assert_error!(
      ancestral_results(Some(&tree)),
      "mutation `A1B2` names 2 sequence positions, expected one"
    );
    Ok(())
  }

  #[test]
  fn test_mugration_results_count_changes_and_uncertain_ancestors() -> Result<(), Report> {
    let auspice = fixture();
    let tree = ResultTree::from_auspice(&auspice)?;

    let root = AncestorState {
      name: o!("root"),
      tips: 3,
      state: o!("brazil"),
      probability: 0.7,
    };
    let expected = MugrationResults {
      attribute: o!("country"),
      states: 2,
      state_changes: vec![StateChange {
        from: o!("brazil"),
        to: o!("china"),
        branches: 2,
      }],
      changed_branches: 2,
      uncertain_below: 0.8,
      uncertain_ancestors: vec![root.clone()],
      root: Some(root),
    };
    assert_eq!(expected, mugration_results(Some(&auspice), Some(&tree), "country")?);
    Ok(())
  }

  pub(crate) mod helpers {
    use indoc::indoc;
    use treetime_io::auspice_types::AuspiceTree;
    use treetime_utils::io::json::json_read_str;

    pub(crate) fn fixture() -> AuspiceTree {
      json_read_str(indoc! {r#"{
        "version": "v2",
        "meta": {
          "colorings": [
            { "key": "num_date", "title": "Date", "type": "continuous" },
            { "key": "country", "title": "Country", "type": "categorical" },
            { "key": "gt", "title": "Genotype", "type": "categorical" }
          ],
          "display_defaults": { "color_by": "country" }
        },
        "tree": {
          "name": "root",
          "node_attrs": {
            "div": 0.0,
            "num_date": { "value": 2010.0, "confidence": [2009.0, 2011.0] },
            "country": { "value": "brazil", "confidence": { "brazil": 0.7, "china": 0.3 } }
          },
          "children": [
            {
              "name": "AB",
              "branch_attrs": { "mutations": { "nuc": ["A10G", "C20T"] } },
              "node_attrs": {
                "div": 0.01,
                "num_date": { "value": 2012.0, "confidence": [2011.5, 2012.5] },
                "country": { "value": "brazil", "confidence": { "brazil": 0.9, "china": 0.1 } }
              },
              "children": [
                {
                  "name": "A",
                  "branch_attrs": { "mutations": { "nuc": ["G10A"] } },
                  "node_attrs": {
                    "div": 0.02,
                    "num_date": { "value": 2015.0 },
                    "bad_branch": { "value": "No" },
                    "country": { "value": "brazil", "confidence": { "brazil": 1.0 } }
                  }
                },
                {
                  "name": "B",
                  "node_attrs": {
                    "div": 0.015,
                    "num_date": { "value": 2014.0 },
                    "bad_branch": { "value": "No" },
                    "country": { "value": "china", "confidence": { "china": 1.0 } }
                  }
                }
              ]
            },
            {
              "name": "C",
              "branch_attrs": { "mutations": { "nuc": ["C20T", "G30-"] } },
              "node_attrs": {
                "div": 0.03,
                "num_date": { "value": 2016.0 },
                "bad_branch": { "value": "Yes" },
                "country": { "value": "china" }
              }
            }
          ]
        }
      }"#})
      .unwrap()
    }
  }
}
