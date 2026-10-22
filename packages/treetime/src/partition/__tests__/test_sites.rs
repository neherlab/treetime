#[cfg(test)]
mod tests {
  use crate::partition::sites::SiteRange;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use treetime_utils::assert_error;
  use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};

  #[rustfmt::skip]
  #[rstest]
  #[case::one_site(     "7",       vec![6]            )]
  #[case::range(        "3-6",     vec![2, 3, 4, 5]   )]
  #[case::step(         "2-10\\3", vec![1, 4, 7]      )]
  #[case::step_on_last( "1-7\\3",  vec![0, 3, 6]      )]
  #[case::step_of_one(  "4-5\\1",  vec![3, 4]         )]
  #[case::one_site_step("5-5\\2",  vec![4]            )]
  #[trace]
  fn test_site_range_columns_are_zero_based(#[case] text: &str, #[case] expected: Vec<usize>) -> Result<(), Report> {
    let range: SiteRange = text.parse()?;
    assert_eq!(expected, range.columns().collect_vec());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::position_zero(   "0-5",                     "site range '0-5': sites are numbered from 1")]
  #[case::first_after_last("9-5",                     "site range '9-5': the first site 9 is after the last site 5")]
  #[case::step_zero(       "1-9\\0",                  "site range '1-9\\0': the step must be at least 1")]
  #[case::word(            "abc",                     "site range 'abc': write a (one site), a-b (the sites a to b) or a-b\\k (every k-th site from a to b), with whole numbers a, b and k")]
  #[case::open_end(        "1-",                      "site range '1-': write a (one site), a-b (the sites a to b) or a-b\\k (every k-th site from a to b), with whole numbers a, b and k")]
  #[case::step_without_end("5\\2",                    "site range '5\\2': write a (one site), a-b (the sites a to b) or a-b\\k (every k-th site from a to b), with whole numbers a, b and k")]
  #[case::slash_step(      "1-9/3",                   "site range '1-9/3': write a (one site), a-b (the sites a to b) or a-b\\k (every k-th site from a to b), with whole numbers a, b and k")]
  #[case::spaces(          " 1-9",                    "site range ' 1-9': write a (one site), a-b (the sites a to b) or a-b\\k (every k-th site from a to b), with whole numbers a, b and k")]
  #[case::empty(           "",                        "site range '': write a (one site), a-b (the sites a to b) or a-b\\k (every k-th site from a to b), with whole numbers a, b and k")]
  #[case::too_large(       "1-99999999999999999999999", "site range '1-99999999999999999999999': write a (one site), a-b (the sites a to b) or a-b\\k (every k-th site from a to b), with whole numbers a, b and k")]
  #[trace]
  fn test_site_range_errors_name_the_text(#[case] text: &str, #[case] expected: &str) {
    assert_error!(text.parse::<SiteRange>().map_err(Report::new), expected);
  }

  #[test]
  fn test_site_range_serializes_as_written() -> Result<(), Report> {
    let ranges: Vec<SiteRange> = json_read_str(r#"["7", "3-6", "2-10\\3"]"#)?;
    assert_eq!(r#"["7","3-6","2-10\\3"]"#, json_write_str(&ranges, JsonPretty(false))?);
    Ok(())
  }

  #[test]
  fn test_site_range_deserialization_reports_the_range_error() {
    assert_error!(
      json_read_str::<SiteRange>(r#""9-5""#),
      "When parsing JSON: Unexpected: invalid value: site range '9-5': the first site 9 is after the last site 5 at line \
       1 column 1"
    );
  }

  proptest! {
    #[test]
    fn test_prop_site_range_text_round_trips(text in generators::site_range_text()) {
      let range: SiteRange = text.parse().unwrap();
      prop_assert_eq!(&text, &range.to_string());
      prop_assert_eq!(range, range.to_string().parse::<SiteRange>().unwrap());
    }

    #[test]
    fn test_prop_site_range_columns_lie_between_first_and_last_with_the_step(
      (first, last, step) in generators::bounds(),
    ) {
      let range: SiteRange = format!("{first}-{last}\\{step}").parse().unwrap();
      let expected = (first..=last).step_by(step).map(|site| site - 1).collect_vec();
      prop_assert_eq!(expected, range.columns().collect_vec());
    }
  }

  mod generators {
    use proptest::prelude::*;

    pub(super) fn bounds() -> impl Strategy<Value = (usize, usize, usize)> {
      (1_usize..1_000_000, 0_usize..1_000, 1_usize..10).prop_map(|(first, extra, step)| (first, first + extra, step))
    }

    pub(super) fn site_range_text() -> impl Strategy<Value = String> {
      bounds().prop_map(|(first, last, step)| match (first == last, step) {
        (true, 1) => format!("{first}"),
        (false, 1) => format!("{first}-{last}"),
        _ => format!("{first}-{last}\\{step}"),
      })
    }
  }
}
