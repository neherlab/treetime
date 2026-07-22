#[cfg(test)]
mod tests {
  use crate::array::serde::array2_as_vec;
  use crate::io::json::{JsonPretty, json_write_str};
  use ndarray::{Array2, array};
  use pretty_assertions::assert_eq;
  use serde::Serialize;

  #[test]
  fn serde_serializes_non_contiguous_array2() {
    let value = Matrix {
      values: array![[1, 2], [3, 4]].reversed_axes(),
    };
    let expected = r#"{"values":[[1,3],[2,4]]}"#;

    let actual = json_write_str(&value, JsonPretty(false)).unwrap();

    assert_eq!(expected, actual);
  }

  #[derive(Serialize)]
  struct Matrix {
    #[serde(serialize_with = "array2_as_vec")]
    values: Array2<i32>,
  }
}
