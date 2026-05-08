//! Proptest generators for phylogenetic tree topologies.
//!
//! Three strategies for generating random trees:
//! - Bisection: tends toward balanced trees
//! - Joining: uniform over tree histories
//! - Caterpillar: maximally unbalanced (for edge cases)

use crate::ancestral::__tests__::prop_generators::branch_length::arb_branch_length;
use crate::representation::payload::ancestral::GraphAncestral;
use proptest::prelude::*;
use std::collections::BTreeSet;
use treetime_graph::node::Named;
use treetime_io::nwk::nwk_read_str;

/// Format a subtree with branch length, wrapping in parens only if it's a compound subtree.
fn format_subtree(subtree: &str, bl: f64) -> String {
  if subtree.contains(',') {
    // Compound subtree - needs parentheses
    format!("({subtree}):{bl}")
  } else {
    // Single taxon - no parentheses
    format!("{subtree}:{bl}")
  }
}

/// Validate that a tree contains all taxa exactly once and has no single-child nodes.
///
/// Returns true if tree is valid, false otherwise.
fn is_valid_tree(newick: &str, expected_taxa: &[String]) -> bool {
  // Check all taxa present exactly once
  for taxon in expected_taxa {
    let count = newick.matches(taxon.as_str()).count();
    if count != 1 {
      return false;
    }
  }

  // Check for single-child nodes: "(X)" where X has no comma at depth 0
  // This is a simplified check that catches the common case
  let bytes = newick.as_bytes();
  let mut i = 0;
  while i < bytes.len() {
    if bytes[i] == b'(' {
      // Find matching close paren and check for comma at depth 1
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
      // If this is not the outermost root paren and has no comma, it's invalid
      if start > 0 && !has_comma_at_depth_1 {
        return false;
      }
    } else {
      i += 1;
    }
  }

  true
}

/// Generate a random tree using bisection (tends toward balanced trees).
///
/// Recursively splits taxa into two non-empty groups and joins them.
fn arb_tree_bisection_inner(taxa: Vec<String>) -> BoxedStrategy<String> {
  let n = taxa.len();
  if n == 1 {
    // Single taxon - just return the name (no branch length at this level)
    Just(taxa[0].clone()).boxed()
  } else if n == 2 {
    // Two taxa - create a simple cherry
    (arb_branch_length(), arb_branch_length())
      .prop_map(move |(bl1, bl2)| {
        let left = format_subtree(&taxa[0], bl1);
        let right = format_subtree(&taxa[1], bl2);
        format!("{left},{right}")
      })
      .boxed()
  } else {
    // Split into two non-empty groups
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

/// Generate a random tree using joining (uniform over histories).
///
/// Starts with n leaf nodes and repeatedly joins two random subtrees.
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
    // Pick two distinct indices to join
    let indices: Vec<usize> = (0..n).collect();
    prop::sample::subsequence(indices, 2)
      .prop_flat_map(move |pair| {
        let i = pair[0].min(pair[1]);
        let j = pair[0].max(pair[1]);
        let taxa_clone = taxa.clone();

        (arb_branch_length(), arb_branch_length()).prop_flat_map(move |(bl1, bl2)| {
          // Create new subtree from joining taxa[i] and taxa[j]
          let left = format_subtree(&taxa_clone[i], bl1);
          let right = format_subtree(&taxa_clone[j], bl2);
          let new_subtree = format!("{left},{right}");

          // Create remaining taxa list with the new subtree replacing the joined pair
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

/// Generate a caterpillar tree (maximally unbalanced).
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
    // Take first taxon and recurse on the rest
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

/// Generate a random tree topology as Newick string (inner part, without root).
///
/// Combines bisection, joining, and caterpillar strategies.
/// Filters out invalid trees that proptest shrinking might produce.
pub fn arb_tree_topology(taxa: Vec<String>) -> BoxedStrategy<String> {
  let taxa_for_filter = taxa.clone();
  prop_oneof![
    2 => arb_tree_bisection_inner(taxa.clone()),
    2 => arb_tree_joining_inner(taxa.clone()),
    1 => arb_tree_caterpillar_inner(taxa),
  ]
  .prop_filter("valid tree topology", move |tree| {
    // Wrap in root format for validation (matching arb_newick output)
    let newick = format!("({tree})root:0.001;");
    is_valid_tree(&newick, &taxa_for_filter)
  })
  .boxed()
}

/// Generate a complete Newick string with root and semicolon.
pub fn arb_newick(n_taxa: usize) -> impl Strategy<Value = String> {
  let taxa: Vec<String> = (0..n_taxa).map(|i| format!("T{i}")).collect();
  arb_tree_topology(taxa).prop_map(|tree| format!("({tree})root:0.001;"))
}

/// Generate taxa names for a given count.
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
      // Newick should end with semicolon
      prop_assert!(newick.ends_with(';'), "Newick should end with semicolon: {newick}");

      // Should contain all taxa names
      for i in 0..4 {
        prop_assert!(newick.contains(&format!("T{i}")), "Missing taxon T{i} in {newick}");
      }

      // Balanced parentheses
      let opens = newick.chars().filter(|&c| c == '(').count();
      let closes = newick.chars().filter(|&c| c == ')').count();
      prop_assert_eq!(opens, closes, "Unbalanced parentheses in {}", newick);
    }

    #[test]
    fn test_prop_tree_arb_newick_has_branch_lengths(newick in arb_newick(3)) {
      // Should contain colons for branch lengths
      prop_assert!(newick.contains(':'), "Newick should contain branch lengths: {newick}");
    }

    #[test]
    fn test_prop_tree_arb_newick_all_taxa_present(newick in arb_newick(5)) {
      // Verify all 5 taxa appear exactly once
      for i in 0..5 {
        let taxon = format!("T{i}");
        let count = newick.matches(&taxon).count();
        prop_assert_eq!(count, 1, "Taxon {} should appear exactly once in {}", taxon, newick);
      }
    }

    #[test]
    fn test_prop_tree_arb_newick_no_double_branch_lengths(newick in arb_newick(4)) {
      // Should not have patterns like ":0.001:0.002" (double branch lengths)
      let chars: Vec<char> = newick.chars().collect();
      for i in 0..chars.len() {
        if chars[i] == ':' {
          // Find the end of this branch length (next comma, paren, or semicolon)
          let mut j = i + 1;
          while j < chars.len() && chars[j] != ',' && chars[j] != ')' && chars[j] != ';' && chars[j] != '(' {
            j += 1;
          }
          // Check if next char after branch length is another colon
          if j < chars.len() && chars[j] == ':' {
            prop_assert!(false, "Double branch length at position {} in {}", i, newick);
          }
        }
      }
    }

    #[test]
    #[allow(clippy::string_slice)] // Newick is ASCII-only, byte indices are safe
    fn test_prop_tree_arb_newick_no_single_element_parens(newick in arb_newick(4)) {
      // Should not have patterns like "(Tx)" without comma inside
      // Find all parenthesized groups and check they have commas
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
              // Root is special - it's allowed to have just the tree content
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
      let graph: GraphAncestral = nwk_read_str(&newick).unwrap();

      let mut actual_names = Vec::new();
      for leaf in graph.get_leaves() {
        let leaf = leaf.read_arc();
        let payload = leaf.payload().read_arc();
        let maybe_name = payload.name().map(|name| name.as_ref().to_owned());
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
