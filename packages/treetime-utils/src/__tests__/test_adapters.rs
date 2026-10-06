#[cfg(test)]
mod tests {
  use crate::adapters::{Array2Rows, ArrayVec, TrueOrNull};
  use crate::assert_error;
  use crate::io::json::{JsonPretty, json_read_str, json_write_str};
  use deser::{Deserialize, Serialize};
  use eyre::Report;
  use ndarray::{Array1, Array2, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[test]
  fn test_adapters_array_vec_round_trips_as_a_list() -> Result<(), Report> {
    let value = Vector {
      values: array![0.5, -1.0, 2.25],
    };
    let expected = r#"{"values":[0.5,-1.0,2.25]}"#;

    let written = json_write_str(&value, JsonPretty(false))?;

    assert_eq!(expected, written);
    assert_eq!(value, json_read_str::<Vector>(&written)?);
    Ok(())
  }

  #[test]
  fn test_adapters_array2_rows_writes_the_rows_of_a_non_contiguous_array() -> Result<(), Report> {
    let value = Matrix {
      values: array![[1, 2], [3, 4]].reversed_axes(),
    };
    let expected = r#"{"values":[[1,3],[2,4]]}"#;

    let written = json_write_str(&value, JsonPretty(false))?;

    assert_eq!(expected, written);
    assert_eq!(value, json_read_str::<Matrix>(&written)?);
    Ok(())
  }

  #[test]
  fn test_adapters_array2_rows_rejects_rows_of_different_lengths() {
    assert_error!(
      json_read_str::<Matrix>(r#"{"values":[[1,2],[3]]}"#),
      "When parsing JSON: Unexpected: invalid value: ShapeError/OutOfBounds: out of bounds indexing at line 1 column 21 (path: values)"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::set(  Flags { outlier: true  }, r#"{"outlier":true}"#)]
  #[case::unset(Flags { outlier: false }, r#"{"outlier":null}"#)]
  #[trace]
  fn test_adapters_true_or_null_round_trips(#[case] value: Flags, #[case] expected: &str) -> Result<(), Report> {
    let written = json_write_str(&value, JsonPretty(false))?;

    assert_eq!(expected, written);
    assert_eq!(value, json_read_str::<Flags>(&written)?);
    Ok(())
  }

  #[test]
  fn test_adapters_true_or_null_reads_a_missing_value_as_false() -> Result<(), Report> {
    assert_eq!(Flags { outlier: false }, json_read_str::<Flags>("{}")?);
    Ok(())
  }

  #[derive(Debug, PartialEq, Serialize, Deserialize)]
  struct Vector {
    #[deser(as = ArrayVec)]
    values: Array1<f64>,
  }

  #[derive(Debug, PartialEq, Serialize, Deserialize)]
  struct Matrix {
    #[deser(as = Array2Rows)]
    values: Array2<i32>,
  }

  #[derive(Debug, PartialEq, Serialize, Deserialize)]
  struct Flags {
    #[deser(default, as = TrueOrNull)]
    outlier: bool,
  }
}
