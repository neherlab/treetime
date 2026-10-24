#[cfg(test)]
mod tests {
  use crate::clock::gaussian::{GaussianClock, GaussianTree};
  use crate::pretty_assert_ulps_eq;
  use crate::timetree::inference::result::NodeTimes;
  use eyre::Report;
  use std::collections::BTreeMap;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::assert_error;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_relative_eq};

  const SEQUENCE_LENGTH: usize = 1000;
  const RATE: f64 = 1e-3;

  fn run(
    nwk: &str,
    dates: &[(&str, f64)],
    excluded: &[&str],
  ) -> Result<(GaussianClock, BTreeMap<String, (f64, f64)>), Report> {
    let parsed = nwk_read(nwk.as_bytes())?;
    let names = parsed.names.clone();
    let key_of: BTreeMap<&str, GraphNodeKey> = parsed
      .graph
      .get_nodes()
      .map(|node| (&names[node.key()], node.key()))
      .collect();
    let mut tip_dates: NodeTimes = parsed.graph.get_nodes().map(|node| (node.key(), None)).collect();
    for (name, date) in dates {
      tip_dates.insert(key_of[name], Some(*date));
    }
    let excluded = excluded.iter().map(|name| (key_of[name], true)).collect();
    let clock = GaussianTree::new(
      &parsed.graph,
      &parsed.branch_lengths,
      SEQUENCE_LENGTH,
      &tip_dates,
      &excluded,
    )?
    .clock(RATE)?;
    let times = clock
      .node_times
      .iter()
      .map(|(key, t)| (names[*key].to_owned(), (t.mean, t.std)))
      .collect();
    Ok((clock, times))
  }

  const TREE_3: &str = "((A:0.010,B:0.020)V:0.015,C:0.030)R;";
  const DATES_3: [(&str, f64); 3] = [("A", 2020.0), ("B", 2015.0), ("C", 2022.0)];

  struct Tree3 {
    nu: [f64; 4],
    t: [f64; 3],
  }

  impl Tree3 {
    fn new() -> Self {
      let nu = [0.010, 0.020, 0.015, 0.030].map(|d| d * SEQUENCE_LENGTH.to_f64_exact() + 1.0);
      Self {
        nu,
        t: DATES_3.map(|(_, t)| t),
      }
    }

    fn eta(&self, rate: f64) -> f64 {
      let d: f64 = self.nu.iter().map(|nu| nu - 1.0).sum();
      rate * SEQUENCE_LENGTH.to_f64_exact() * (d + 4.0) / d
    }

    fn joint(&self, eta: f64) -> ([[f64; 2]; 2], [f64; 2]) {
      let [na, nb, nv, nc] = self.nu;
      let [ta, tb, tc] = self.t;
      let e2 = eta * eta;
      let q = [
        [e2 * (1.0 / nc + 1.0 / nv), -e2 / nv],
        [-e2 / nv, e2 * (1.0 / na + 1.0 / nb + 1.0 / nv)],
      ];
      let h = [
        e2 * ((tc - nc / eta) / nc - 1.0 / eta),
        e2 * ((ta - na / eta) / na + (tb - nb / eta) / nb + 1.0 / eta),
      ];
      (q, h)
    }

    fn log_likelihood(&self, rate: f64, t_r: f64) -> f64 {
      let eta = self.eta(rate);
      let [na, nb, nv, nc] = self.nu;
      let [ta, tb, tc] = self.t;
      let targets = [(ta - na / eta, na), (tb - nb / eta, nb), (t_r + nv / eta, nv)];
      let t_v = targets.iter().map(|(c, nu)| c / nu).sum::<f64>() / targets.iter().map(|(_, nu)| 1.0 / nu).sum::<f64>();
      let term = |tau: f64, nu: f64| (eta * tau - nu).powi(2) / (2.0 * nu);
      3.0 * eta.ln() - term(ta - t_v, na) - term(tb - t_v, nb) - term(t_v - t_r, nv) - term(tc - t_r, nc)
    }
  }

  fn invert(m: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    [[m[1][1] / det, -m[0][1] / det], [-m[1][0] / det, m[0][0] / det]]
  }

  #[test]
  fn test_gaussian_clock_node_marginals_match_joint_gaussian() -> Result<(), Report> {
    let tree = Tree3::new();
    let (q, h) = tree.joint(tree.eta(RATE));
    let cov = invert(q);
    let mean = [cov[0][0] * h[0] + cov[0][1] * h[1], cov[1][0] * h[0] + cov[1][1] * h[1]];

    let (_, times) = run(TREE_3, &DATES_3, &[])?;
    pretty_assert_ulps_eq!(times["R"].0, mean[0], max_ulps = 64);
    pretty_assert_ulps_eq!(times["R"].1, cov[0][0].sqrt(), max_ulps = 64);
    pretty_assert_ulps_eq!(times["V"].0, mean[1], max_ulps = 64);
    pretty_assert_ulps_eq!(times["V"].1, cov[1][1].sqrt(), max_ulps = 64);
    Ok(())
  }

  #[test]
  fn test_gaussian_clock_rate_and_root_std_match_finite_difference_hessian() -> Result<(), Report> {
    let tree = Tree3::new();
    let (clock, times) = run(TREE_3, &DATES_3, &[])?;
    let t_r = times["R"].0;

    let (hm, ht) = (1e-4 * RATE, 1e-2);
    let ll = |dm: f64, dt: f64| tree.log_likelihood(RATE + dm, t_r + dt);
    let l00 = ll(0.0, 0.0);
    let i_mm = -(ll(hm, 0.0) - 2.0 * l00 + ll(-hm, 0.0)) / (hm * hm);
    let i_tt = -(ll(0.0, ht) - 2.0 * l00 + ll(0.0, -ht)) / (ht * ht);
    let i_cross = -(ll(hm, ht) - ll(hm, -ht) - ll(-hm, ht) + ll(-hm, -ht)) / (4.0 * hm * ht);
    let cov = invert([[i_mm, i_cross], [i_cross, i_tt]]);

    pretty_assert_relative_eq!(clock.rate_std.unwrap(), cov[0][0].sqrt(), max_relative = 1e-7);
    pretty_assert_relative_eq!(clock.root_time_std.unwrap(), cov[1][1].sqrt(), max_relative = 1e-7);
    Ok(())
  }

  #[test]
  fn test_gaussian_rate_maximizes_the_profile_likelihood() -> Result<(), Report> {
    let tree = Tree3::new();
    let parsed = nwk_read(TREE_3.as_bytes())?;
    let names = parsed.names.clone();
    let tip_dates: NodeTimes = parsed
      .graph
      .get_nodes()
      .map(|node| {
        let date = DATES_3
          .iter()
          .find(|(tip, _)| *tip == &names[node.key()])
          .map(|(_, t)| *t);
        (node.key(), date)
      })
      .collect();
    let rate = GaussianTree::new(
      &parsed.graph,
      &parsed.branch_lengths,
      SEQUENCE_LENGTH,
      &tip_dates,
      &BTreeMap::new(),
    )?
    .rate()?;

    let profile = |rate: f64| {
      let (q, h) = tree.joint(tree.eta(rate));
      let cov = invert(q);
      tree.log_likelihood(rate, cov[0][0] * h[0] + cov[0][1] * h[1])
    };
    let step = 1e-4 * rate;
    let slope = (profile(rate + step) - profile(rate - step)) / (2.0 * step);
    let curvature = (profile(rate + step) - 2.0 * profile(rate) + profile(rate - step)) / (step * step);
    pretty_assert_abs_diff_eq!(0.0, slope / (curvature.abs() * rate), epsilon = 1e-6);
    assert!(profile(rate) > profile(0.99 * rate) && profile(rate) > profile(1.01 * rate));
    Ok(())
  }

  #[test]
  fn test_gaussian_clock_excluded_or_undated_tip_equals_pruned_tip() -> Result<(), Report> {
    let full = "((A:0.010,B:0.020)X:0.015,(C:0.030,D:0.010)Y:0.020)R;";
    let pruned = "((A:0.010,B:0.020)X:0.015,(C:0.030)Y:0.020)R;";
    let dates = [("A", 2020.0), ("B", 2015.0), ("C", 2022.0)];
    let with_d = [dates.as_slice(), &[("D", 2018.0)]].concat();

    let (expected, expected_times) = run(pruned, &dates, &[])?;
    for (clock, times) in [run(full, &with_d, &["D"])?, run(full, &dates, &[])?] {
      pretty_assert_ulps_eq!(clock.rate_std.unwrap(), expected.rate_std.unwrap(), max_ulps = 16);
      pretty_assert_ulps_eq!(
        clock.root_time_std.unwrap(),
        expected.root_time_std.unwrap(),
        max_ulps = 16
      );
      assert_eq!(
        times.keys().filter(|name| *name != "D").collect::<Vec<_>>(),
        expected_times.keys().collect::<Vec<_>>()
      );
      for (name, (mean, std)) in &expected_times {
        pretty_assert_ulps_eq!(times[name].0, *mean, max_ulps = 16);
        pretty_assert_ulps_eq!(times[name].1, *std, max_ulps = 16);
      }
    }
    Ok(())
  }

  #[test]
  fn test_gaussian_clock_rejects_tree_without_dated_tips() {
    assert_error!(
      run(TREE_3, &[], &[]),
      "Gaussian clock needs substitutions on the dated part of the tree"
    );
  }
}
