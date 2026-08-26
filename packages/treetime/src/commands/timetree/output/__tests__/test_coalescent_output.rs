#[cfg(test)]
mod tests {
  use crate::commands::timetree::output::coalescent::{
    CoalescentBand, CoalescentInputs, CoalescentOutput, CoalescentOutputMode, CoalescentSegmentRow, CoalescentSolve,
    Estimate, SegmentInterval, coalescent_delimited_str, coalescent_json_str,
  };
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;
  use treetime_utils::io::json::json_read_str;

  // Oracle for N_e: N_e = T_c * gen_per_year, and the band scales by the same factor
  // (packages/treetime/src/coalescent/population_size.rs).
  const GEN_PER_YEAR: f64 = 50.0;

  #[test]
  fn test_coalescent_output_rows_map_rich_to_flat() -> Result<(), Report> {
    // Skyline: two segments with a band. The rich -> flat projection shifts the 0-based index to
    // 1-based and flattens the nested segment/T_c/N_e objects to dotted scalar columns.
    let output = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Skyline,
        gen_per_year: GEN_PER_YEAR,
        confidence_n_std: Some(2.0),
      },
      &CoalescentSolve {
        segment_boundaries: &[0.0, 5.0, 10.0],
        tc_values: &[2.0, 3.0],
        band: Some(CoalescentBand {
          lower: &[1.0, 1.5],
          upper: &[4.0, 6.0],
        }),
      },
    )?;

    // Rich indices stay 0-based.
    assert_eq!(
      vec![0, 1],
      output.outputs.segments.iter().map(|s| s.index).collect::<Vec<_>>()
    );

    let expected = vec![
      CoalescentSegmentRow {
        index: 1,
        segment_start: 0.0,
        segment_end: 5.0,
        tc_value: 2.0,
        tc_lower: Some(1.0),
        tc_upper: Some(4.0),
        ne_value: 100.0,
        ne_lower: Some(50.0),
        ne_upper: Some(200.0),
      },
      CoalescentSegmentRow {
        index: 2,
        segment_start: 5.0,
        segment_end: 10.0,
        tc_value: 3.0,
        tc_lower: Some(1.5),
        tc_upper: Some(6.0),
        ne_value: 150.0,
        ne_lower: Some(75.0),
        ne_upper: Some(300.0),
      },
    ];
    assert_eq!(expected, output.rows());
    Ok(())
  }

  #[test]
  fn test_coalescent_output_ne_scales_tc_by_gen_per_year() -> Result<(), Report> {
    let output = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Constant,
        gen_per_year: GEN_PER_YEAR,
        confidence_n_std: Some(2.0),
      },
      &CoalescentSolve {
        segment_boundaries: &[0.0, 10.0],
        tc_values: &[2.0],
        band: Some(CoalescentBand {
          lower: &[1.0],
          upper: &[4.0],
        }),
      },
    )?;

    let segment = &output.outputs.segments[0];
    // N_e = T_c * gen_per_year for the point and both band bounds.
    assert_eq!(Estimate::with_band(2.0, 1.0, 4.0), segment.tc);
    assert_eq!(Estimate::with_band(100.0, 50.0, 200.0), segment.ne);
    Ok(())
  }

  #[test]
  fn test_coalescent_output_fixed_has_no_band() -> Result<(), Report> {
    // A fixed, user-supplied T_c carries no band, so T_c and the derived N_e report a value only.
    let output = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Fixed,
        gen_per_year: GEN_PER_YEAR,
        confidence_n_std: None,
      },
      &CoalescentSolve {
        segment_boundaries: &[0.0, 10.0],
        tc_values: &[2.0],
        band: None,
      },
    )?;

    let segment = &output.outputs.segments[0];
    assert_eq!(SegmentInterval { start: 0.0, end: 10.0 }, segment.segment);
    assert_eq!(Estimate::point(2.0), segment.tc);
    assert_eq!(Estimate::point(100.0), segment.ne);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::fixed_one_segment(   CoalescentOutputMode::Fixed,    &[0.0, 10.0],            &[2.0],           None,                                          1)]
  #[case::constant_one_segment(CoalescentOutputMode::Constant, &[0.0, 10.0],            &[2.0],           Some((&[1.0][..], &[4.0][..])),                1)]
  #[case::skyline_three(       CoalescentOutputMode::Skyline,  &[0.0, 3.0, 6.0, 9.0],   &[2.0, 3.0, 4.0], Some((&[1.0, 1.0, 1.0][..], &[4.0, 5.0, 6.0][..])), 3)]
  #[trace]
  fn test_coalescent_output_segment_count(
    #[case] mode: CoalescentOutputMode,
    #[case] boundaries: &[f64],
    #[case] tc_values: &[f64],
    #[case] band: Option<(&[f64], &[f64])>,
    #[case] expected_rows: usize,
  ) -> Result<(), Report> {
    let output = CoalescentOutput::new(
      CoalescentInputs {
        mode,
        gen_per_year: GEN_PER_YEAR,
        confidence_n_std: band.map(|_| 2.0),
      },
      &CoalescentSolve {
        segment_boundaries: boundaries,
        tc_values,
        band: band.map(|(lower, upper)| CoalescentBand { lower, upper }),
      },
    )?;
    assert_eq!(expected_rows, output.rows().len());
    Ok(())
  }

  #[test]
  fn test_coalescent_output_json_roundtrips_and_renames_keys() -> Result<(), Report> {
    let output = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Constant,
        gen_per_year: GEN_PER_YEAR,
        confidence_n_std: Some(2.0),
      },
      &CoalescentSolve {
        segment_boundaries: &[0.0, 10.0],
        tc_values: &[2.0],
        band: Some(CoalescentBand {
          lower: &[1.0],
          upper: &[4.0],
        }),
      },
    )?;

    let json = coalescent_json_str(&output)?;
    // Serde renames apply: the rich document uses T_c/N_e, never the Rust field names.
    assert!(json.contains("\"inputs\""));
    assert!(json.contains("\"outputs\""));
    assert!(json.contains("\"segments\""));
    assert!(json.contains("\"T_c\""));
    assert!(json.contains("\"N_e\""));
    assert!(!json.contains("\"tc\""));
    assert!(!json.contains("\"ne\""));

    let parsed: CoalescentOutput = json_read_str(&json)?;
    assert_eq!(output, parsed);
    Ok(())
  }

  #[test]
  fn test_coalescent_output_json_fixed_omits_band_keys() -> Result<(), Report> {
    let output = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Fixed,
        gen_per_year: GEN_PER_YEAR,
        confidence_n_std: None,
      },
      &CoalescentSolve {
        segment_boundaries: &[0.0, 10.0],
        tc_values: &[2.0],
        band: None,
      },
    )?;

    let json = coalescent_json_str(&output)?;
    // No band, so the optional bound keys and confidence are omitted entirely.
    assert!(!json.contains("\"lower\""));
    assert!(!json.contains("\"upper\""));
    assert!(!json.contains("\"confidence_n_std\""));
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::csv(b',',  "index,segment.start,segment.end,T_c.value,T_c.lower,T_c.upper,N_e.value,N_e.lower,N_e.upper")]
  #[case::tsv(b'\t', "index\tsegment.start\tsegment.end\tT_c.value\tT_c.lower\tT_c.upper\tN_e.value\tN_e.lower\tN_e.upper")]
  #[trace]
  fn test_coalescent_output_delimited_headers_are_dotted(
    #[case] delimiter: u8,
    #[case] expected_header: &str,
  ) -> Result<(), Report> {
    let output = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Constant,
        gen_per_year: GEN_PER_YEAR,
        confidence_n_std: Some(2.0),
      },
      &CoalescentSolve {
        segment_boundaries: &[0.0, 10.0],
        tc_values: &[2.0],
        band: Some(CoalescentBand {
          lower: &[1.0],
          upper: &[4.0],
        }),
      },
    )?;

    let serialized = coalescent_delimited_str(&output, delimiter)?;
    let header = serialized.lines().next().unwrap();
    assert_eq!(expected_header, header);
    Ok(())
  }

  #[test]
  fn test_coalescent_output_delimited_rows_are_1_based_with_empty_band() -> Result<(), Report> {
    // Constant: full band. Fixed: empty band columns (None serializes to an empty field).
    let constant = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Constant,
        gen_per_year: GEN_PER_YEAR,
        confidence_n_std: Some(2.0),
      },
      &CoalescentSolve {
        segment_boundaries: &[0.0, 10.0],
        tc_values: &[2.0],
        band: Some(CoalescentBand {
          lower: &[1.0],
          upper: &[4.0],
        }),
      },
    )?;
    let csv = coalescent_delimited_str(&constant, b',')?;
    assert_eq!("1,0.0,10.0,2.0,1.0,4.0,100.0,50.0,200.0", csv.lines().nth(1).unwrap());

    let fixed = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Fixed,
        gen_per_year: GEN_PER_YEAR,
        confidence_n_std: None,
      },
      &CoalescentSolve {
        segment_boundaries: &[0.0, 10.0],
        tc_values: &[2.0],
        band: None,
      },
    )?;
    let csv = coalescent_delimited_str(&fixed, b',')?;
    assert_eq!("1,0.0,10.0,2.0,,,100.0,,", csv.lines().nth(1).unwrap());
    Ok(())
  }

  #[test]
  fn test_coalescent_output_rejects_boundary_length_mismatch() {
    let result = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Skyline,
        gen_per_year: GEN_PER_YEAR,
        confidence_n_std: Some(2.0),
      },
      &CoalescentSolve {
        segment_boundaries: &[0.0, 5.0], // needs 3 boundaries for 2 segments
        tc_values: &[2.0, 3.0],
        band: None,
      },
    );
    assert_error!(
      result,
      "Coalescent output expects 3 segment boundaries for 2 segment(s), got 2"
    );
  }

  #[test]
  fn test_coalescent_output_rejects_band_length_mismatch() {
    let result = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Skyline,
        gen_per_year: GEN_PER_YEAR,
        confidence_n_std: Some(2.0),
      },
      &CoalescentSolve {
        segment_boundaries: &[0.0, 5.0, 10.0],
        tc_values: &[2.0, 3.0],
        band: Some(CoalescentBand {
          lower: &[1.0], // one bound for two segments
          upper: &[4.0, 6.0],
        }),
      },
    );
    assert_error!(
      result,
      "Coalescent band bounds must have one entry per segment (2), got lower=1, upper=2"
    );
  }
}
