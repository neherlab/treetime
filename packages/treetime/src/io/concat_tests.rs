#![allow(clippy::iter_on_single_items, clippy::redundant_type_annotations)]
use crate::io::concat::*;
use rstest::rstest;
use std::io::Read;

#[rstest]
fn test_concatenate_with_delimiter_both_with_trailing_newline() {
  let r1: &[u8] = b"First\nreader\n";
  let r2: &[u8] = b"Second\nreader\n";

  let mut concat = Concat::with_delimiter([r1, r2].into_iter(), Some(b"\n".to_vec()));
  let mut result = String::new();
  concat.read_to_string(&mut result).unwrap();

  assert_eq!(result, "First\nreader\n\nSecond\nreader\n\n");
}

#[rstest]
fn test_concatenate_with_delimiter_one_with_trailing_newline() {
  let r1: &[u8] = b"First\nreader\n";

  let mut concat = Concat::with_delimiter([r1].into_iter(), Some(b"\n".to_vec()));
  let mut result = String::new();
  concat.read_to_string(&mut result).unwrap();

  assert_eq!(result, "First\nreader\n\n");
}

#[rstest]
fn test_concatenate_with_delimiter_one_without_trailing_newline() {
  let r1: &[u8] = b"First\nreader\nwithout\ntrailing\nnewline";

  let mut concat = Concat::with_delimiter([r1].into_iter(), Some(b"\n".to_vec()));
  let mut result = String::new();
  concat.read_to_string(&mut result).unwrap();

  assert_eq!(result, "First\nreader\nwithout\ntrailing\nnewline\n");
}

#[rstest]
fn test_concatenate_with_delimiter_one_empty() {
  let r1: &[u8] = b"";

  let mut concat = Concat::with_delimiter([r1].into_iter(), Some(b"\n".to_vec()));
  let mut result = String::new();
  concat.read_to_string(&mut result).unwrap();

  assert_eq!(result, "\n");
}

#[rstest]
fn test_concatenate_with_delimiter_one_empty_with_trailing_newline() {
  let r1: &[u8] = b"\n";

  let mut concat = Concat::with_delimiter([r1].into_iter(), Some(b"\n".to_vec()));
  let mut result = String::new();
  concat.read_to_string(&mut result).unwrap();

  assert_eq!(result, "\n\n");
}

#[rstest]
fn test_concatenate_with_delimiter_no_newline() {
  let r1: &[u8] = b"No\ntrailing\nnewline";
  let r2: &[u8] = b"And\nneither\nhere";

  let mut concat = Concat::with_delimiter([r1, r2].into_iter(), Some(b"\n".to_vec()));
  let mut result = String::new();
  concat.read_to_string(&mut result).unwrap();

  assert_eq!(result, "No\ntrailing\nnewline\nAnd\nneither\nhere\n");
}

#[rstest]
fn test_concatenate_with_delimiter_first_empty_no_newline() {
  let r1: &[u8] = b"";
  let r2: &[u8] = b"Second";

  let mut concat = Concat::with_delimiter([r1, r2].into_iter(), Some(b"\n".to_vec()));
  let mut result = String::new();
  concat.read_to_string(&mut result).unwrap();

  assert_eq!(result, "\nSecond\n");
}

#[rstest]
fn test_concatenate_with_delimiter_second_empty_no_newline() {
  let r1: &[u8] = b"First";
  let r2: &[u8] = b"";

  let mut concat = Concat::with_delimiter([r1, r2].into_iter(), Some(b"\n".to_vec()));
  let mut result = String::new();
  concat.read_to_string(&mut result).unwrap();

  assert_eq!(result, "First\n\n");
}

#[rstest]
fn test_concatenate_with_delimiter_both_empty() {
  let r1: &[u8] = b"";
  let r2: &[u8] = b"";

  let mut concat = Concat::with_delimiter([r1, r2].into_iter(), Some(b"\n".to_vec()));
  let mut result = String::new();
  concat.read_to_string(&mut result).unwrap();

  assert_eq!(result, "\n\n");
}
