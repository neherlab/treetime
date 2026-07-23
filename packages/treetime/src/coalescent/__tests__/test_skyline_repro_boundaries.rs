#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::load_date_constraints;
  use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
  use crate::partition::timetree::GraphTimetree;
  use eyre::Report;
  use maplit::btreemap;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::nwk::nwk_read_str;

  const SMALL_TREE_NWK: &str = "((leaf1:1.0,leaf2:1.0)internal1:1.0,leaf3:1.0)root:1.0;";

  fn small_tree_dates() -> DatesMap {
    btreemap! {
      "root".to_owned() => Some(DateConstraint::exact(2000.0)),
      "internal1".to_owned() => Some(DateConstraint::exact(2005.0)),
      "leaf1".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf2".to_owned() => Some(DateConstraint::exact(2010.0)),
      "leaf3".to_owned() => Some(DateConstraint::exact(2012.0)),
    }
  }

  fn small_tree() -> Result<GraphTimetree, Report> {
    let graph = nwk_read_str(SMALL_TREE_NWK)?;
    load_date_constraints(&small_tree_dates(), &graph)?;
    Ok(graph)
  }

  // Counterexample for "every segment owns mergers" / ascending boundaries.
  //
  // The small tree has only two distinct merger (internal-node) times, 2000 and
  // 2005, but this requests 5 segments. `merger_quantile_boundaries` clamps each
  // duplicate candidate against the previous boundary instead of coalescing merger
  // events, so the returned grid is `[2000, 2000, 2000, 2002.5, 2005, 2012]`: two
  // zero-width segments carrying no mergers. Zero-width segments are unreachable via
  // `segment_index`, yet their `z_i` still participate in the smoothing penalty, so
  // the fitted trajectory depends on phantom transitions.
  //
  // Expected: strictly increasing segment boundaries, or a rejected `n_points` when
  // the tree cannot support that many merger-owning segments.
  #[test]
  fn test_skyline_segment_boundaries_strictly_increasing() -> Result<(), Report> {
    let params = SkylineParams {
      n_points: 5,
      ..SkylineParams::default()
    };

    let result = optimize_skyline(&small_tree()?, &params)?;
    let bounds = result.segment_boundaries.to_vec();
    let strictly_increasing = bounds.windows(2).all(|w| w[1] > w[0]);

    assert!(
      strictly_increasing,
      "segment boundaries must be strictly increasing (every segment spans a positive interval), got {bounds:?}"
    );
    Ok(())
  }
}
