#[cfg(test)]
mod tests {
  use crate::rtt_chart::__tests__::test_support::tests::{estimated_model, model_with_rate_std, tip};
  use crate::rtt_chart::data::RttChartData;
  use crate::rtt_chart::frame::RttFrame;
  use crate::rtt_chart::layout::{CHART_SIZE, LabelRole, LegendSymbol, Rect, RttLayout};
  use crate::rtt_chart::summary::RttSummary;
  use helpers::{Scenario, inside, intersects, layout_of};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_charts::font::register_chart_fonts;
  use treetime_utils::cast::ToF64Exact;

  #[rustfmt::skip]
  #[rstest]
  #[case::root_close_with_outliers(Scenario::RootClose)]
  #[case::root_far_before(          Scenario::RootFar)]
  #[case::negative_rate(            Scenario::NegativeRate)]
  #[case::large_counts(             Scenario::LargeCounts)]
  #[trace]
  fn test_rtt_layout_texts_do_not_overlap_each_other_or_the_plot(#[case] scenario: Scenario) {
    let layout = layout_of(scenario);
    let boxes = layout
      .labels
      .iter()
      .map(|label| (label.text.clone(), label.bounds()))
      .chain(layout.legend.iter().flat_map(|entry| {
        [
          (format!("{:?} symbol", entry.symbol), entry.symbol_box),
          (entry.label.text.clone(), entry.label.bounds()),
        ]
      }))
      .collect_vec();
    let overlapping = boxes
      .iter()
      .tuple_combinations()
      .filter(|((_, a), (_, b))| intersects(a, b))
      .map(|((a, _), (b, _))| (a.clone(), b.clone()))
      .collect_vec();
    let on_plot = layout
      .labels
      .iter()
      .filter(|label| label.role != LabelRole::RootMarker && intersects(&label.bounds(), &layout.plot))
      .map(|label| label.text.clone())
      .chain(
        layout
          .legend
          .iter()
          .filter(|entry| {
            intersects(&entry.label.bounds(), &layout.plot) || intersects(&entry.symbol_box, &layout.plot)
          })
          .map(|entry| entry.label.text.clone()),
      )
      .collect_vec();
    let canvas = Rect {
      left: 0.0,
      top: 0.0,
      right: f64::from(CHART_SIZE.0),
      bottom: f64::from(CHART_SIZE.1),
    };
    let off_canvas = boxes
      .iter()
      .filter(|(_, bounds)| !inside(bounds, &canvas))
      .map(|(text, _)| text.clone())
      .collect_vec();
    assert_eq!(
      (Vec::<(String, String)>::new(), Vec::<String>::new(), Vec::<String>::new()),
      (overlapping, on_plot, off_canvas)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::root_close_with_outliers(Scenario::RootClose)]
  #[case::root_far_before(          Scenario::RootFar)]
  #[case::negative_rate(            Scenario::NegativeRate)]
  #[case::large_counts(             Scenario::LargeCounts)]
  #[trace]
  fn test_rtt_layout_draws_every_sample_marker_whole_inside_the_plot(#[case] scenario: Scenario) {
    let (data, layout) = (scenario.data(), layout_of(scenario));
    let radius = f64::from(layout.marker_radius + 1);
    let cut = data
      .dated()
      .map(|point| (layout.dates.map(point.date), layout.divergences.map(point.div)))
      .filter(|&(x, y)| {
        x - radius < layout.plot.left
          || x + radius > layout.plot.right
          || y - radius < layout.plot.top
          || y + radius > layout.plot.bottom
      })
      .collect_vec();
    assert_eq!(Vec::<(f64, f64)>::new(), cut);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::root_close_with_outliers(Scenario::RootClose,    vec![LegendSymbol::Sample, LegendSymbol::Excluded, LegendSymbol::Fit, LegendSymbol::Extrapolation])]
  #[case::root_far_before(          Scenario::RootFar,      vec![LegendSymbol::Sample, LegendSymbol::Fit, LegendSymbol::Extrapolation])]
  #[case::negative_rate(            Scenario::NegativeRate, vec![LegendSymbol::Sample, LegendSymbol::Fit])]
  #[trace]
  fn test_rtt_layout_legend_lists_only_what_the_chart_draws(
    #[case] scenario: Scenario,
    #[case] expected: Vec<LegendSymbol>,
  ) {
    let symbols = layout_of(scenario).legend.iter().map(|entry| entry.symbol).collect_vec();
    assert_eq!(expected, symbols);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::root_close_with_outliers(Scenario::RootClose)]
  #[case::large_counts(             Scenario::LargeCounts)]
  #[trace]
  fn test_rtt_layout_keeps_the_root_label_off_the_drawn_lines(#[case] scenario: Scenario) {
    let (data, model) = (scenario.data(), scenario.model());
    let frame = RttFrame::new(&data, &model);
    let layout = layout_of(scenario);
    let label = layout
      .labels
      .iter()
      .find(|label| label.role == LabelRole::RootMarker)
      .unwrap()
      .bounds();
    let crossing = frame
      .extrapolation
      .iter()
      .chain([&frame.fit])
      .flat_map(|&[(x0, y0), (x1, y1)]| {
        (0..=100).map(move |step| {
          let t = step.to_f64_exact() / 100.0;
          (x0 + t * (x1 - x0), y0 + t * (y1 - y0))
        })
      })
      .map(|(date, div)| (layout.dates.map(date), layout.divergences.map(div)))
      .filter(|&(x, y)| x >= label.left && x <= label.right && y >= label.top && y <= label.bottom)
      .collect_vec();
    assert_eq!(Vec::<(f64, f64)>::new(), crossing);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::root_close_with_outliers(Scenario::RootClose,    Some("Root 2016.71"))]
  #[case::root_far_before(          Scenario::RootFar,      Some("\u{2190} Root date 1814.2"))]
  #[case::negative_rate(            Scenario::NegativeRate, None)]
  #[trace]
  fn test_rtt_layout_marks_the_root(#[case] scenario: Scenario, #[case] expected: Option<&str>) {
    let layout = layout_of(scenario);
    let root = layout
      .labels
      .iter()
      .find(|label| matches!(label.role, LabelRole::RootMarker | LabelRole::RootNote))
      .map(|label| label.text.as_str());
    assert_eq!(expected, root);
  }

  mod helpers {
    use super::*;
    use std::iter::repeat_with;
    use treetime::clock::clock_model::ClockModel;
    use treetime::clock::rtt::ClockRegressionResult;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Scenario {
      RootClose,
      RootFar,
      NegativeRate,
      LargeCounts,
    }

    impl Scenario {
      pub(super) fn data(self) -> RttChartData {
        RttChartData::new(&self.results()).unwrap()
      }

      fn results(self) -> Vec<ClockRegressionResult> {
        let model = self.model();
        let (first, span, count) = match self {
          Self::RootClose => (2017.77, 6.85, 460),
          Self::RootFar | Self::NegativeRate => (1944.5, 73.2, 48),
          Self::LargeCounts => (2019.9, 4.5, 12_000),
        };
        let date = |index: usize| first + span * index.to_f64_exact() / count.to_f64_exact();
        let noise = |index: usize| ((index * 7919) % 101).to_f64_exact() / 100.0 - 0.5;
        let fitted = (0..=count).map(|index| {
          let div = model.div(date(index)).abs() * (1.0 + 0.1 * noise(index));
          tip(Some(date(index)), div, false)
        });
        let (excluded, undated) = match self {
          Self::RootClose => (2, 3),
          Self::RootFar | Self::NegativeRate => (0, 0),
          Self::LargeCounts => (1234, 54_321),
        };
        let excluded = (0..excluded).map(|index| tip(Some(date(index * 7)), model.div(date(index * 7)) * 0.4, true));
        let undated = repeat_with(|| tip(None, 0.0, false)).take(undated);
        fitted.chain(excluded).chain(undated).collect()
      }

      pub(super) fn model(self) -> ClockModel {
        match self {
          Self::RootClose => estimated_model(5.072_949e-5, -0.102_306_463_7, 0.96),
          Self::RootFar => estimated_model(8.47e-4, -1.5366, 0.95),
          Self::NegativeRate => estimated_model(-2.0e-4, 0.55, -0.4),
          Self::LargeCounts => model_with_rate_std(5.07e-4, -1.0236, 3.1e-5),
        }
      }
    }

    pub(super) fn intersects(a: &Rect, b: &Rect) -> bool {
      a.left < b.right && b.left < a.right && a.top < b.bottom && b.top < a.bottom
    }

    pub(super) fn inside(bounds: &Rect, canvas: &Rect) -> bool {
      bounds.left >= canvas.left
        && bounds.right <= canvas.right
        && bounds.top >= canvas.top
        && bounds.bottom <= canvas.bottom
    }

    pub(super) fn layout_of(scenario: Scenario) -> RttLayout {
      register_chart_fonts().unwrap();
      let (data, model) = (scenario.data(), scenario.model());
      let frame = RttFrame::new(&data, &model);
      let summary = RttSummary::new(&data, &model);
      RttLayout::new(&data, &frame, &summary).unwrap()
    }
  }
}
