#[cfg(test)]
mod tests {
  use crate::clock::divergence::root_to_node_divergences;
  use crate::o;
  use crate::pretty_assert_abs_diff_eq;
  use crate::test_utils::root_clock_fit;
  use crate::timetree::init::StartRateFloor;
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::cast::ToF64Exact;

  const TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  #[rstest]
  #[case::all_tips_dated(&[])]
  #[case::one_tip_undated(&["D"])]
  #[trace]
  fn test_init_start_rate_floor_is_a_tenth_of_the_mean_divergence_per_year_of_sampling(
    #[case] undated: &[&str],
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let branch_lengths = nwk_parsed.branch_lengths;
    let divergences = root_to_node_divergences(&nwk_parsed.graph, |edge_key| {
      branch_lengths[&edge_key].unwrap_or_default()
    })?;
    let divs: BTreeMap<String, f64> = nwk_parsed
      .graph
      .get_leaves()
      .map(|leaf| (names[leaf.key()].to_owned(), divergences[&leaf.key()]))
      .collect();
    let dates: BTreeMap<String, f64> = helpers::leaf_dates()
      .into_iter()
      .filter(|(name, _)| !undated.contains(&name.as_str()))
      .collect();

    let floor = StartRateFloor::of_fit(&root_clock_fit(TREE, &dates)?);

    let mean_div = dates.keys().map(|name| divs[name]).sum::<f64>() / dates.len().to_f64_exact();
    let (first, last) = dates.values().copied().minmax().into_option().expect("dated tips");
    let StartRateFloor::Signal(actual) = floor else {
      panic!("expected a signal floor, got {floor:?}");
    };
    pretty_assert_abs_diff_eq!(0.1 * mean_div / (last - first), actual, epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_init_start_rate_floor_requires_only_a_positive_rate_without_divergence() -> Result<(), Report> {
    let fit = root_clock_fit("((A:0,B:0)AB:0,(C:0,D:0)CD:0)root:0;", &helpers::leaf_dates())?;

    assert_eq!(StartRateFloor::Positive, StartRateFloor::of_fit(&fit));
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn leaf_dates() -> BTreeMap<String, f64> {
      btreemap! {
        o!("A") => 2013.0,
        o!("B") => 2022.0,
        o!("C") => 2017.0,
        o!("D") => 2005.0,
      }
    }
  }
}
