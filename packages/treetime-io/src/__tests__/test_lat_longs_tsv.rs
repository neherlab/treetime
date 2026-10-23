#[cfg(test)]
mod tests {
  use crate::lat_longs_tsv::{Coordinates, lat_longs_read};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::tabs(               "country\tUnited States\t38.5\t-97.25\n",         ("country", "United States"), Some((38.5, -97.25)))]
  #[case::white_space(        "country  brazil -10.5   -55\n",                  ("country", "brazil"),        Some((-10.5, -55.0)))]
  #[case::comments_and_blanks("# resolution name lat long\n\n  \ncountry peru -9 -75\n", ("country", "peru"), Some((-9.0, -75.0)))]
  #[case::lower_case_match(   "Country\tBrazil\t-10\t-55\n",                    ("COUNTRY", "BRAZIL"),        Some((-10.0, -55.0)))]
  #[case::later_line_wins(    "country\tperu\t1\t2\ncountry\tPeru\t-9\t-75\n", ("country", "peru"),          Some((-9.0, -75.0)))]
  #[case::other_resolution(   "region\tperu\t-9\t-75\n",                        ("country", "peru"),          None)]
  #[case::range_ends(         "country\tpole\t-90\t180\n",                      ("country", "pole"),          Some((-90.0, 180.0)))]
  #[trace]
  fn test_lat_longs_read_finds_the_coordinates_of_a_location(
    #[case] input: &str,
    #[case] (resolution, name): (&str, &str),
    #[case] expected: Option<(f64, f64)>,
  ) -> Result<(), Report> {
    let expected = expected.map(|(latitude, longitude)| Coordinates { latitude, longitude });

    let actual = lat_longs_read(input.as_bytes())?.get(resolution, name);

    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::three_fields(        "country\tperu\t-9\n",                 "Line 1: expected 4 fields (resolution, name, latitude and longitude), found 3: 'country\tperu\t-9'")]
  #[case::spaces_in_name(      "country United States 38 -97\n",     "Line 1: expected 4 fields (resolution, name, latitude and longitude), found 5: 'country United States 38 -97'")]
  #[case::latitude_not_number( "country\tperu\tnine\t-75\n",          "Line 1: latitude 'nine' is not a number: invalid float literal")]
  #[case::latitude_too_large(  "country\tperu\t-9\t-75\nx\ty\t95\t0\n", "Line 2: latitude 95 is outside the range -90 to 90")]
  #[case::longitude_too_small( "country\tperu\t-9\t-181\n",           "Line 1: longitude -181 is outside the range -180 to 180")]
  #[case::latitude_not_finite( "country\tperu\tnan\t0\n",             "Line 1: latitude nan is outside the range -90 to 90")]
  #[trace]
  fn test_lat_longs_read_stops_at_an_invalid_line(#[case] input: &str, #[case] expected: &str) {
    assert_error!(lat_longs_read(input.as_bytes()), expected);
  }
}
