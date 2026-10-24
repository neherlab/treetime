#[cfg(test)]
mod tests {
  use crate::reroot::__tests__::test_midpoint::tests::helpers::{diameter, midpoint};
  use generators::arb_newick;
  use proptest::prelude::*;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::prop_assert_ulps_eq;

  proptest! {
    #[test]
    fn test_prop_midpoint_roots_at_half_the_unchanged_diameter(newick in arb_newick()) {
      let parsed = nwk_read(newick.as_bytes()).unwrap();
      let before = diameter(&parsed.graph, &parsed.branch_lengths).unwrap();

      let once = midpoint(&newick, 1).unwrap();
      let twice = midpoint(&newick, 2).unwrap();

      let farthest = once.root_to_tip.values().copied().fold(f64::NEG_INFINITY, f64::max);
      prop_assert_ulps_eq!(before / 2.0, farthest, max_ulps = 16);
      prop_assert_ulps_eq!(before, once.diameter, max_ulps = 16);
      prop_assert_eq!(
        (once.root_name, once.root_children, once.root_to_tip),
        (twice.root_name, twice.root_children, twice.root_to_tip)
      );
    }
  }

  mod generators {
    use proptest::prelude::*;

    pub(super) fn arb_newick() -> impl Strategy<Value = String> {
      (2_usize..10)
        .prop_flat_map(|tips| {
          prop::collection::vec((any::<usize>(), any::<usize>(), 1e-3..1.0_f64, 1e-3..1.0_f64), tips - 1)
        })
        .prop_map(|merges| {
          let tips = merges.len() + 1;
          let subtrees = (0..tips).map(|index| format!("T{index}")).collect::<Vec<_>>();
          let root = merges.into_iter().fold(
            subtrees,
            |mut subtrees, (first, second, first_length, second_length)| {
              let left = subtrees.swap_remove(first % subtrees.len());
              let right = subtrees.swap_remove(second % subtrees.len());
              subtrees.push(format!("({left}:{first_length},{right}:{second_length})"));
              subtrees
            },
          );
          format!("{}R;", root.concat())
        })
    }
  }
}
