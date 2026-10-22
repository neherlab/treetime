#[cfg(test)]
mod tests {
  use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
  use crate::commands::timetree::coalescent::{CoalescentArg, parse_coalescent};
  use crate::config::schema::command_schema;
  use crate::config::schema_check::schema_diagnostics;
  use crate::yaml::yaml_read_str;
  use clap::Parser;
  use clap::error::ErrorKind;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use treetime::timetree::coalescent_timescale::{CoalescentMode, FixedTimescale};
  use treetime_utils::assert_error;
  use treetime_utils::io::json::{json_read_str, to_json_value};

  #[rustfmt::skip]
  #[rstest]
  #[case::substitutions(   "0.01",    CoalescentMode::Fixed(FixedTimescale::SubstitutionsPerSite(0.01)))]
  #[case::exponent(        "1e-2",    CoalescentMode::Fixed(FixedTimescale::SubstitutionsPerSite(0.01)))]
  #[case::years(           "10y",     CoalescentMode::Fixed(FixedTimescale::Years(10.0)))]
  #[case::fractional_years("2.56y",   CoalescentMode::Fixed(FixedTimescale::Years(2.56)))]
  #[case::weeks(           "26w",     CoalescentMode::Fixed(FixedTimescale::Years(26.0 / (365.25 / 7.0))))]
  #[case::largest_weeks(   "1e308w",  CoalescentMode::Fixed(FixedTimescale::Years(1e308 / (365.25 / 7.0))))]
  #[case::days(            "90d",     CoalescentMode::Fixed(FixedTimescale::Years(90.0 / 365.25)))]
  #[case::opt(             "opt",     CoalescentMode::Constant)]
  #[case::skyline(         "skyline", CoalescentMode::Skyline)]
  #[trace]
  fn test_coalescent_flag_and_config_parse_into_one_mode(
    #[case] text: &str,
    #[case] expected: CoalescentMode,
  ) -> Result<(), Report> {
    let from_flag = TreetimeTimetreeArgsRaw::try_parse_from(["treetime", "--coalescent", text])?;
    let from_config: TreetimeTimetreeArgsRaw = yaml_read_str(&format!("coalescent: {text}\n"))?;

    assert_eq!(Some(expected), from_flag.coalescent.map(CoalescentArg::mode));
    assert_eq!(Some(expected), from_config.coalescent.map(CoalescentArg::mode));
    Ok(())
  }

  #[test]
  fn test_coalescent_unset_disables_the_prior() -> Result<(), Report> {
    let args = TreetimeTimetreeArgsRaw::try_parse_from(["treetime"])?;

    assert_eq!(None, args.coalescent);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(            "0")]
  #[case::space_before_unit("10 y")]
  #[case::unknown_unit(    "10m")]
  #[case::negative(        "-1y")]
  #[case::infinite(        "inf")]
  #[case::unit_only(       "y")]
  #[case::empty(           "")]
  #[case::capitalized(     "Opt")]
  #[case::days_underflow(  "1e-323d")]
  #[case::weeks_underflow( "1e-323w")]
  #[trace]
  fn test_coalescent_rejects_values_outside_the_accepted_forms(#[case] text: &str) {
    let expected = format!(
      "'{text}' is not a coalescent prior. Give a fixed time scale as a positive number in substitutions per site (0.01, \
       1e-2) or as a positive number directly followed by the time unit y, w or d (10y, 26w, 90d); opt to fit a \
       constant time scale; or skyline to fit a piecewise-constant one"
    );

    assert_error!(parse_coalescent(text).map_err(Report::new), expected);
    let flag_error = TreetimeTimetreeArgsRaw::try_parse_from(["treetime".to_owned(), format!("--coalescent={text}")])
      .err()
      .map(|error| error.kind());
    assert_eq!(Some(ErrorKind::ValueValidation), flag_error);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::substitutions(        "0.01",    (true,  true))]
  #[case::exponent(             "1e-2",    (true,  true))]
  #[case::leading_dot(          ".5",      (true,  true))]
  #[case::trailing_dot(         "5.",      (true,  true))]
  #[case::plus_sign(            "+1y",     (true,  true))]
  #[case::zero_before_digit(    "0.05d",   (true,  true))]
  #[case::leading_zeros(        "007w",    (true,  true))]
  #[case::zero(                 "0",       (false, false))]
  #[case::zero_years(           "0y",      (false, false))]
  #[case::zero_with_exponent(   "0.0e5d",  (false, false))]
  #[case::zeros_with_dot(       "00.000",  (false, false))]
  #[case::zero_dot(             "0.",      (false, false))]
  #[case::dot_zero_weeks(       ".0w",     (false, false))]
  #[case::negative(             "-1y",     (false, false))]
  #[case::infinite(             "inf",     (false, false))]
  #[case::space_before_unit(    "10 y",    (false, false))]
  #[case::unknown_unit(         "10m",     (false, false))]
  #[case::overflow(             "1e999y",  (false, true))]
  #[case::underflow(            "1e-400",  (false, true))]
  #[trace]
  fn test_coalescent_schema_pattern_agrees_with_the_parser_except_out_of_range_numbers(
    #[case] text: &str,
    #[case] expected: (bool, bool),
  ) {
    let schema = command_schema::<TreetimeTimetreeArgsRaw>();

    let parser_accepts = parse_coalescent(text).is_ok();
    let schema_accepts = schema_diagnostics(&json!({ "coalescent": text }), &schema, false).is_empty();

    assert_eq!(expected, (parser_accepts, schema_accepts));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(    json!(0.0),  "0")]
  #[case::negative(json!(-0.5), "-0.5")]
  #[case::boolean( json!(true), "true")]
  #[trace]
  fn test_coalescent_config_rejects_values_outside_the_accepted_forms(#[case] value: Value, #[case] text: &str) {
    let config = json!({ "coalescent": value }).to_string();
    let expected = format!(
      "When parsing JSON: Unexpected: invalid value: '{text}' is not a coalescent prior. Give a fixed time scale as a \
       positive number in substitutions per site (0.01, 1e-2) or as a positive number directly followed by the time \
       unit y, w or d (10y, 26w, 90d); opt to fit a constant time scale; or skyline to fit a piecewise-constant one at \
       line 1 column 15 (path: coalescent)"
    );

    assert_error!(json_read_str::<TreetimeTimetreeArgsRaw>(config), expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::number_stays_a_number( json!(0.01),      json!(0.01))]
  #[case::number_text(           json!("0.01"),    json!(0.01))]
  #[case::time_keeps_its_unit(   json!("26w"),     json!("26w"))]
  #[case::option(                json!("skyline"), json!("skyline"))]
  #[trace]
  fn test_coalescent_config_writes_the_form_the_user_gave(
    #[case] value: Value,
    #[case] expected: Value,
  ) -> Result<(), Report> {
    let args: TreetimeTimetreeArgsRaw = json_read_str(json!({ "coalescent": value }).to_string())?;

    let written = to_json_value(&args)?;

    assert_eq!(Some(&expected), written.get("coalescent"));
    Ok(())
  }
}
