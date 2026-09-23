use crate::ancestral::__tests__::prop_generators::branch_length::arb_branch_length;
use proptest::prelude::*;
use std::collections::BTreeSet;
use treetime_graph::graph::Graph;
use treetime_io::nwk::nwk_read_str;

fn format_subtree(subtree: &str, bl: f64) -> String {
  if subtree.contains(',') {
    format!("({subtree}):{bl}")
  } else {
    format!("{subtree}:{bl}")
  }
}

fn is_valid_tree(newick: &str, expected_taxa: &[String]) -> bool {
  for taxon in expected_taxa {
    let count = newick.matches(taxon.as_str()).count();
    if count != 1 {
      return false;
    }
  }

  let bytes = newick.as_bytes();
  let mut i = 0;
  while i < bytes.len() {
    if bytes[i] == b'(' {
      let start = i;
      let mut depth = 1;
      let mut has_comma_at_depth_1 = false;
      i += 1;
      while i < bytes.len() && depth > 0 {
        match bytes[i] {
          b'(' => depth += 1,
          b')' => depth -= 1,
          b',' if depth == 1 => has_comma_at_depth_1 = true,
          _ => {},
        }
        i += 1;
      }
      if start > 0 && !has_comma_at_depth_1 {
        return false;
      }
    } else {
      i += 1;
    }
  }

  true
}

fn arb_tree_bisection_inner(taxa: Vec<String>) -> BoxedStrategy<String> {
  let n = taxa.len();
  if n == 1 {
    Just(taxa[0].clone()).boxed()
  } else if n == 2 {
    (arb_branch_length(), arb_branch_length())
      .prop_map(move |(bl1, bl2)| {
        let left = format_subtree(&taxa[0], bl1);
        let right = format_subtree(&taxa[1], bl2);
        format!("{left},{right}")
      })
      .boxed()
  } else {
    (1..n)
      .prop_flat_map(move |split| {
        let left_taxa: Vec<String> = taxa[..split].to_vec();
        let right_taxa: Vec<String> = taxa[split..].to_vec();
        (
          arb_tree_bisection_inner(left_taxa),
          arb_tree_bisection_inner(right_taxa),
          arb_branch_length(),
          arb_branch_length(),
        )
          .prop_map(|(left, right, bl_left, bl_right)| {
            let left_str = format_subtree(&left, bl_left);
            let right_str = format_subtree(&right, bl_right);
            format!("{left_str},{right_str}")
          })
      })
      .boxed()
  }
}

fn arb_tree_joining_inner(taxa: Vec<String>) -> BoxedStrategy<String> {
  let n = taxa.len();
  if n == 1 {
    Just(taxa[0].clone()).boxed()
  } else if n == 2 {
    (arb_branch_length(), arb_branch_length())
      .prop_map(move |(bl1, bl2)| {
        let left = format_subtree(&taxa[0], bl1);
        let right = format_subtree(&taxa[1], bl2);
        format!("{left},{right}")
      })
      .boxed()
  } else {
    let indices: Vec<usize> = (0..n).collect();
    prop::sample::subsequence(indices, 2)
      .prop_flat_map(move |pair| {
        let i = pair[0].min(pair[1]);
        let j = pair[0].max(pair[1]);
        let taxa_clone = taxa.clone();

        (arb_branch_length(), arb_branch_length()).prop_flat_map(move |(bl1, bl2)| {
          let left = format_subtree(&taxa_clone[i], bl1);
          let right = format_subtree(&taxa_clone[j], bl2);
          let new_subtree = format!("{left},{right}");

          let mut remaining: Vec<String> = Vec::with_capacity(n - 1);
          for (idx, t) in taxa_clone.iter().enumerate() {
            if idx != i && idx != j {
              remaining.push(t.clone());
            }
          }
          remaining.push(new_subtree);

          arb_tree_joining_inner(remaining)
        })
      })
      .boxed()
  }
}

fn arb_tree_caterpillar_inner(taxa: Vec<String>) -> BoxedStrategy<String> {
  let n = taxa.len();
  if n == 1 {
    Just(taxa[0].clone()).boxed()
  } else if n == 2 {
    (arb_branch_length(), arb_branch_length())
      .prop_map(move |(bl1, bl2)| {
        let left = format_subtree(&taxa[0], bl1);
        let right = format_subtree(&taxa[1], bl2);
        format!("{left},{right}")
      })
      .boxed()
  } else {
    let first = taxa[0].clone();
    let rest: Vec<String> = taxa[1..].to_vec();

    (arb_branch_length(), arb_tree_caterpillar_inner(rest))
      .prop_flat_map(move |(bl, subtree)| {
        let first_clone = first.clone();
        arb_branch_length().prop_map(move |bl2| {
          let first_str = format_subtree(&first_clone, bl);
          let rest_str = format_subtree(&subtree, bl2);
          format!("{first_str},{rest_str}")
        })
      })
      .boxed()
  }
}

pub fn arb_tree_topology(taxa: Vec<String>) -> BoxedStrategy<String> {
  let taxa_for_filter = taxa.clone();
  prop_oneof![
    2 => arb_tree_bisection_inner(taxa.clone()),
    2 => arb_tree_joining_inner(taxa.clone()),
    1 => arb_tree_caterpillar_inner(taxa),
  ]
  .prop_filter("valid tree topology", move |tree| {
    let newick = format!("({tree})root:0.001;");
    is_valid_tree(&newick, &taxa_for_filter)
  })
  .boxed()
}

fn arb_newick(n_taxa: usize) -> impl Strategy<Value = String> {
  let taxa: Vec<String> = (0..n_taxa).map(|i| format!("T{i}")).collect();
  arb_tree_topology(taxa).prop_map(|tree| format!("({tree})root:0.001;"))
}

pub fn taxa_names(n: usize) -> Vec<String> {
  (0..n).map(|i| format!("T{i}")).collect()
}

#[cfg(test)]
mod tests {
  use super::*;
  use proptest::proptest;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn test_prop_tree_arb_newick_valid_syntax(newick in arb_newick(4)) {
      prop_assert!(newick.ends_with(';'), "Newick should end with semicolon: {newick}");

      for i in 0..4 {
        prop_assert!(newick.contains(&format!("T{i}")), "Missing taxon T{i} in {newick}");
      }

      let opens = newick.chars().filter(|&c| c == '(').count();
      let closes = newick.chars().filter(|&c| c == ')').count();
      prop_assert_eq!(opens, closes, "Unbalanced parentheses in {}", newick);
    }

    #[test]
    fn test_prop_tree_arb_newick_has_branch_lengths(newick in arb_newick(3)) {
      prop_assert!(newick.contains(':'), "Newick should contain branch lengths: {newick}");
    }

    #[test]
    fn test_prop_tree_arb_newick_all_taxa_present(newick in arb_newick(5)) {
      for i in 0..5 {
        let taxon = format!("T{i}");
        let count = newick.matches(&taxon).count();
        prop_assert_eq!(count, 1, "Taxon {} should appear exactly once in {}", taxon, newick);
      }
    }

    #[test]
    fn test_prop_tree_arb_newick_no_double_branch_lengths(newick in arb_newick(4)) {
      let chars: Vec<char> = newick.chars().collect();
      for i in 0..chars.len() {
        if chars[i] == ':' {
          let mut j = i + 1;
          while j < chars.len() && chars[j] != ',' && chars[j] != ')' && chars[j] != ';' && chars[j] != '(' {
            j += 1;
          }
          if j < chars.len() && chars[j] == ':' {
            prop_assert!(false, "Double branch length at position {} in {}", i, newick);
          }
        }
      }
    }

    #[test]
    #[expect(clippy::string_slice, reason = "indices come from char_indices")]
    fn test_prop_tree_arb_newick_no_single_element_parens(newick in arb_newick(4)) {
      let mut depth = 0;
      let mut start = 0;
      for (i, c) in newick.char_indices() {
        match c {
          '(' => {
            if depth == 0 {
              start = i;
            }
            depth += 1;
          }
          ')' => {
            depth -= 1;
            if depth == 0 {
              let content = &newick[start + 1..i];
              if start > 0 {
                prop_assert!(
                  content.contains(','),
                  "Found single-element parentheses at {}..{} in {}: ({})",
                  start, i, newick, content
                );
              }
            }
          }
          _ => {}
        }
      }
    }

    #[test]
    fn test_prop_tree_arb_newick_parseable_and_leaf_names_exact(newick in arb_newick(6)) {
      let nwk_parsed = nwk_read_str(&newick).unwrap();
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let graph: Graph = graph;

      let mut actual_names = Vec::new();
      for leaf in graph.get_leaves() {
        let maybe_name = names.get(&leaf.key()).cloned().flatten();
        prop_assert!(maybe_name.is_some(), "Leaf node is missing name in Newick: {newick}");
        if let Some(name) = maybe_name {
          actual_names.push(name);
        }
      }

      let expected_names: BTreeSet<String> = (0..6).map(|index| format!("T{index}")).collect();
      let actual_name_count = actual_names.len();
      let actual_names: BTreeSet<String> = actual_names.into_iter().collect();

      prop_assert_eq!(actual_name_count, 6, "Expected 6 leaves in parsed tree: {}", newick);
      prop_assert_eq!(
        actual_names,
        expected_names,
        "Parsed leaf names should match generated taxa exactly: {}",
        newick
      );
    }
  }
}
