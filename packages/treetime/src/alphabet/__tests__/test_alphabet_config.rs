#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{FILL_CHAR, NON_CHAR, VARIABLE_CHAR};
  use crate::alphabet::alphabet_config::AlphabetConfig;
  use crate::vec_u8;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use indexmap::indexmap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
  use treetime_primitives::AsciiChar;

  fn make_valid_config() -> AlphabetConfig {
    AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'G', 'T'],
      ambiguous: indexmap! {
        b'R' => vec_u8!['A', 'G'],
        b'Y' => vec_u8!['C', 'T'],
      },
      unknown: b'N',
      gap: b'-',
    }
  }

  #[test]
  fn test_alphabet_config_validate_valid() {
    let config = make_valid_config();
    let result = config.validate();
    result.unwrap();
  }

  #[test]
  fn test_alphabet_config_validate_duplicate_canonical() {
    let config = AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'A', 'T'],
      ambiguous: indexmap! {},
      unknown: b'N',
      gap: b'-',
    };
    let result = config.validate();
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("duplicate"));
  }

  #[test]
  fn test_alphabet_config_validate_canonical_ambiguous_overlap() {
    let config = AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'G', 'T'],
      ambiguous: indexmap! {
        b'A' => vec_u8!['C', 'G'],
      },
      unknown: b'N',
      gap: b'-',
    };
    let result = config.validate();
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("disjoint"));
  }

  #[test]
  fn test_alphabet_config_validate_canonical_contains_gap() {
    let config = AlphabetConfig {
      canonical: vec_u8!['A', 'C', '-', 'T'],
      ambiguous: indexmap! {},
      unknown: b'N',
      gap: b'-',
    };
    let result = config.validate();
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("gap"));
  }

  #[test]
  fn test_alphabet_config_validate_canonical_contains_unknown() {
    let config = AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'N', 'T'],
      ambiguous: indexmap! {},
      unknown: b'N',
      gap: b'-',
    };
    let result = config.validate();
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("unknown"));
  }

  #[test]
  fn test_alphabet_config_validate_ambiguous_contains_gap() {
    let config = AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'G', 'T'],
      ambiguous: indexmap! {
        b'-' => vec_u8!['A', 'G'],
      },
      unknown: b'N',
      gap: b'-',
    };
    let result = config.validate();
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("gap"));
  }

  #[test]
  fn test_alphabet_config_validate_ambiguous_maps_to_noncanonical() {
    let config = AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'G', 'T'],
      ambiguous: indexmap! {
        b'R' => vec_u8!['A', 'X'],
      },
      unknown: b'N',
      gap: b'-',
    };
    let result = config.validate();
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("canonical"));
  }

  #[rstest]
  #[case::non_char(u8::from(NON_CHAR))]
  #[case::variable_char(u8::from(VARIABLE_CHAR))]
  #[case::fill_char(u8::from(FILL_CHAR))]
  #[trace]
  fn test_alphabet_config_validate_reserved_in_canonical(#[case] reserved: u8) {
    let config = AlphabetConfig {
      canonical: vec![b'A', b'C', reserved, b'T'],
      ambiguous: indexmap! {},
      unknown: b'N',
      gap: b'-',
    };
    let result = config.validate();
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("reserved"));
  }

  #[rstest]
  #[case::non_char(u8::from(NON_CHAR))]
  #[case::variable_char(u8::from(VARIABLE_CHAR))]
  #[case::fill_char(u8::from(FILL_CHAR))]
  #[trace]
  fn test_alphabet_config_validate_reserved_in_ambiguous_key(#[case] reserved: u8) {
    let config = AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'G', 'T'],
      ambiguous: indexmap! {
        reserved => vec_u8!['A', 'G'],
      },
      unknown: b'N',
      gap: b'-',
    };
    let result = config.validate();
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("reserved"));
  }

  #[rstest]
  #[case::non_char(u8::from(NON_CHAR))]
  #[case::variable_char(u8::from(VARIABLE_CHAR))]
  #[case::fill_char(u8::from(FILL_CHAR))]
  #[trace]
  fn test_alphabet_config_validate_reserved_in_ambiguous_value(#[case] reserved: u8) {
    let config = AlphabetConfig {
      canonical: vec![b'A', b'C', b'G', b'T', reserved],
      ambiguous: indexmap! {
        b'R' => vec![b'A', reserved],
      },
      unknown: b'N',
      gap: b'-',
    };
    let result = config.validate();
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("reserved"));
  }

  #[test]
  fn test_alphabet_config_create_profile_map() {
    let config = make_valid_config();
    let profile_map = config.create_profile_map().unwrap();

    let profile_a = &profile_map[&AsciiChar::from_byte_unchecked(b'A')];
    let expected_a = array![1.0, 0.0, 0.0, 0.0];
    for (actual, expected) in profile_a.iter().zip(expected_a.iter()) {
      pretty_assert_ulps_eq!(*expected, *actual, max_ulps = 4);
    }

    let profile_r = &profile_map[&AsciiChar::from_byte_unchecked(b'R')];
    let expected_r = array![1.0, 0.0, 1.0, 0.0];
    for (actual, expected) in profile_r.iter().zip(expected_r.iter()) {
      pretty_assert_ulps_eq!(*expected, *actual, max_ulps = 4);
    }

    let profile_n = &profile_map[&AsciiChar::from_byte_unchecked(b'N')];
    let expected_n = array![1.0, 1.0, 1.0, 1.0];
    for (actual, expected) in profile_n.iter().zip(expected_n.iter()) {
      pretty_assert_ulps_eq!(*expected, *actual, max_ulps = 4);
    }
  }

  #[test]
  fn test_alphabet_config_create_profile_map_gap_matches_unknown() {
    let config = make_valid_config();
    let profile_map = config.create_profile_map().unwrap();

    let profile_gap = &profile_map[&AsciiChar::from_byte_unchecked(b'-')];
    let profile_n = &profile_map[&AsciiChar::from_byte_unchecked(b'N')];
    for (gap, unk) in profile_gap.iter().zip(profile_n.iter()) {
      pretty_assert_ulps_eq!(*unk, *gap, max_ulps = 4);
    }
  }

  #[test]
  fn test_alphabet_config_create_profile_map_identity_matrix_for_canonical() {
    let config = AlphabetConfig {
      canonical: vec_u8!['X', 'Y', 'Z'],
      ambiguous: indexmap! {},
      unknown: b'?',
      gap: b'-',
    };
    let profile_map = config.create_profile_map().unwrap();

    let profile_x = &profile_map[&AsciiChar::from_byte_unchecked(b'X')];
    let expected_x = array![1.0, 0.0, 0.0];
    for (actual, expected) in profile_x.iter().zip(expected_x.iter()) {
      pretty_assert_ulps_eq!(*expected, *actual, max_ulps = 4);
    }

    let profile_y = &profile_map[&AsciiChar::from_byte_unchecked(b'Y')];
    let expected_y = array![0.0, 1.0, 0.0];
    for (actual, expected) in profile_y.iter().zip(expected_y.iter()) {
      pretty_assert_ulps_eq!(*expected, *actual, max_ulps = 4);
    }

    let profile_z = &profile_map[&AsciiChar::from_byte_unchecked(b'Z')];
    let expected_z = array![0.0, 0.0, 1.0];
    for (actual, expected) in profile_z.iter().zip(expected_z.iter()) {
      pretty_assert_ulps_eq!(*expected, *actual, max_ulps = 4);
    }
  }

  #[test]
  fn test_alphabet_config_serde_roundtrip() -> Result<(), Report> {
    let config = make_valid_config();
    let json = json_write_str(&config, JsonPretty(false))?;
    let deserialized: AlphabetConfig = json_read_str(&json)?;
    assert_eq!(config, deserialized);
    Ok(())
  }

  #[rstest]
  #[case::non_char(u8::from(NON_CHAR))]
  #[case::variable_char(u8::from(VARIABLE_CHAR))]
  #[case::fill_char(u8::from(FILL_CHAR))]
  #[trace]
  fn test_alphabet_config_validate_reserved_as_unknown(#[case] reserved: u8) {
    let config = AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'G', 'T'],
      ambiguous: indexmap! {},
      unknown: reserved,
      gap: b'-',
    };
    let result = config.validate();
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("reserved"));
  }

  #[rstest]
  #[case::non_char(u8::from(NON_CHAR))]
  #[case::variable_char(u8::from(VARIABLE_CHAR))]
  #[case::fill_char(u8::from(FILL_CHAR))]
  #[trace]
  fn test_alphabet_config_validate_reserved_as_gap(#[case] reserved: u8) {
    let config = AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'G', 'T'],
      ambiguous: indexmap! {},
      unknown: b'N',
      gap: reserved,
    };
    let result = config.validate();
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("reserved"));
  }

  #[test]
  fn test_alphabet_config_validate_ambiguous_empty_values() {
    let config = AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'G', 'T'],
      ambiguous: indexmap! {
        b'R' => vec![],
      },
      unknown: b'N',
      gap: b'-',
    };
    let result = config.validate();
    result.unwrap();
  }

  #[test]
  fn test_alphabet_config_validate_empty_canonical() {
    let config = AlphabetConfig {
      canonical: vec![],
      ambiguous: indexmap! {},
      unknown: b'N',
      gap: b'-',
    };
    let result = config.validate();
    result.unwrap();
  }

  #[test]
  fn test_alphabet_config_create_profile_map_all_ambiguous() {
    let config = AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'G', 'T'],
      ambiguous: indexmap! {
        b'R' => vec_u8!['A', 'G'],
        b'Y' => vec_u8!['C', 'T'],
        b'S' => vec_u8!['C', 'G'],
        b'W' => vec_u8!['A', 'T'],
      },
      unknown: b'N',
      gap: b'-',
    };
    let profile_map = config.create_profile_map().unwrap();

    let profile_s = &profile_map[&AsciiChar::from_byte_unchecked(b'S')];
    let expected_s = array![0.0, 1.0, 1.0, 0.0];
    for (actual, expected) in profile_s.iter().zip(expected_s.iter()) {
      pretty_assert_ulps_eq!(*expected, *actual, max_ulps = 4);
    }

    let profile_w = &profile_map[&AsciiChar::from_byte_unchecked(b'W')];
    let expected_w = array![1.0, 0.0, 0.0, 1.0];
    for (actual, expected) in profile_w.iter().zip(expected_w.iter()) {
      pretty_assert_ulps_eq!(*expected, *actual, max_ulps = 4);
    }
  }

  #[test]
  fn test_alphabet_config_equality() {
    let config1 = make_valid_config();
    let config2 = make_valid_config();
    assert_eq!(config1, config2);

    let config3 = AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'G'],
      ambiguous: indexmap! {},
      unknown: b'N',
      gap: b'-',
    };
    assert_ne!(config1, config3);
  }

  #[test]
  fn test_alphabet_config_clone() {
    let config = make_valid_config();
    let cloned = config.clone();
    assert_eq!(config, cloned);
  }

  #[test]
  fn test_alphabet_config_validate_gap_equals_unknown() {
    let config = AlphabetConfig {
      canonical: vec_u8!['A', 'C', 'G', 'T'],
      ambiguous: indexmap! {},
      unknown: b'-',
      gap: b'-',
    };
    let result = config.validate();
    result.unwrap();
  }
}
