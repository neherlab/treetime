#[cfg(test)]
mod tests {
  use crate::annotation::{write_beast_attrs, write_nhx_attrs, write_raw_comments};
  use crate::types::NewickValue;
  use helpers::{parse_node_attrs, write_beast, write_nhx};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

  #[rustfmt::skip]
  #[rstest]
  #[case::canonical_number(     "[&a=0.95]",         NewickValue::Number(0.95))]
  #[case::canonical_integer(    "[&a=12]",           NewickValue::Number(12.0))]
  #[case::leading_zero(         "[&a=0123]",         NewickValue::NumberText("0123".to_owned()))]
  #[case::trailing_zero(        "[&a=2020.50]",      NewickValue::NumberText("2020.50".to_owned()))]
  #[case::exponent(             "[&a=1e5]",          NewickValue::NumberText("1e5".to_owned()))]
  #[case::infinity(             "[&a=-infinity]",    NewickValue::String("-infinity".to_owned()))]
  #[case::single_quoted_comma(  "[&a='x,y']",        NewickValue::String("x,y".to_owned()))]
  #[case::double_quoted_comma(  r#"[&a="x,y"]"#,     NewickValue::String("x,y".to_owned()))]
  #[case::apostrophe_in_word(   "[&a=Cote d'Ivoire]", NewickValue::String("Cote d'Ivoire".to_owned()))]
  #[case::empty_array(          "[&a={}]",           NewickValue::Array(Vec::new()))]
  #[case::quoted_brace(         r#"[&a={"x}",y}]"#,  NewickValue::Array(vec![NewickValue::String("x}".to_owned()), NewickValue::String("y".to_owned())]))]
  #[case::boolean(              "[&a=True]",         NewickValue::Boolean(true))]
  #[trace]
  fn test_annotation_beast_value(#[case] comment: &str, #[case] expected: NewickValue) {
    let attrs = parse_node_attrs(comment);

    assert_eq!(BTreeMap::from([("a".to_owned(), expected)]), attrs);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::comma_in_value(    r#"[&a="x,y",b=1]"#,      vec![("a", NewickValue::String("x,y".to_owned())), ("b", NewickValue::Number(1.0))])]
  #[case::quoted_key_equals( r#"[&"a=b"=1,c=2]"#,      vec![("a=b", NewickValue::Number(1.0)), ("c", NewickValue::Number(2.0))])]
  #[case::nhx_lowercase(     "[&&nhx:S=human]",        vec![("S", NewickValue::String("human".to_owned()))])]
  #[case::nhx_bare_tag(      "[&&NHX:D:S=human]",      vec![("D", NewickValue::Boolean(true)), ("S", NewickValue::String("human".to_owned()))])]
  #[trace]
  fn test_annotation_pairs(#[case] comment: &str, #[case] expected: Vec<(&str, NewickValue)>) {
    let expected: BTreeMap<String, NewickValue> = expected.into_iter().map(|(key, value)| (key.to_owned(), value)).collect();

    assert_eq!(expected, parse_node_attrs(comment));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::key_with_equals(   "a=b",   NewickValue::Number(1.0))]
  #[case::key_with_brace(    "a{b",   NewickValue::Number(1.0))]
  #[case::key_with_quote(    "a\"b",  NewickValue::Number(1.0))]
  #[case::key_with_amp(      "&NHX",  NewickValue::Number(1.0))]
  #[case::value_brackets(    "k",     NewickValue::String("x]y[z".to_owned()))]
  #[case::value_quotes(      "k",     NewickValue::String("say \"hi\", 'you'".to_owned()))]
  #[case::number_huge(       "k",     NewickValue::Number(1e300))]
  #[case::number_text(       "k",     NewickValue::NumberText("1.50".to_owned()))]
  #[case::nested_array(      "k",     NewickValue::Array(vec![NewickValue::Array(vec![NewickValue::Boolean(false)])]))]
  #[trace]
  fn test_annotation_beast_roundtrip(#[case] key: &str, #[case] value: NewickValue) {
    let expected = BTreeMap::from([(key.to_owned(), value), ("z".to_owned(), NewickValue::Number(2.0))]);

    let written = write_beast(&expected).unwrap();

    assert_eq!(expected, parse_node_attrs(&written), "written: {written}");
  }

  #[test]
  fn test_annotation_beast_huge_number_uses_exponent() {
    let attrs = BTreeMap::from([("k".to_owned(), NewickValue::Number(1e300))]);

    assert_eq!("[&k=1.0e300]", write_beast(&attrs).unwrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nan(            NewickValue::Number(f64::NAN),                    "Newick cannot represent the number NaN")]
  #[case::infinite(       NewickValue::Number(f64::INFINITY),               "Newick cannot represent the number inf")]
  #[case::bad_text(       NewickValue::NumberText("1,5".to_owned()),        r#"The annotation value "1,5" is marked as a number, but it is not a finite number"#)]
  #[case::infinite_text(  NewickValue::NumberText("inf".to_owned()),        r#"The annotation value "inf" is marked as a number, but it is not a finite number"#)]
  #[trace]
  fn test_annotation_beast_rejects_unrepresentable_value(#[case] value: NewickValue, #[case] expected: &str) {
    let attrs = BTreeMap::from([("k".to_owned(), value)]);

    assert_eq!(expected, format!("{:#}", write_beast(&attrs).unwrap_err()));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::colon_value(     ("k",   "a:b"),  "NHX cannot represent a value containing a reserved character (':', '=', '[', ']' or '\"'): a:b")]
  #[case::bracket_value(   ("k",   "x[y"),  "NHX cannot represent a value containing a reserved character (':', '=', '[', ']' or '\"'): x[y")]
  #[case::quote_value(     ("k",   "5\""),  "NHX cannot represent a value containing a reserved character (':', '=', '[', ']' or '\"'): 5\"")]
  #[case::colon_key(       ("a:b", "x"),    "NHX cannot represent a key containing a reserved character (':', '=', '[', ']' or '\"'): a:b")]
  #[case::empty_key(       ("",    "x"),    "An NHX annotation needs a non-empty key")]
  #[trace]
  fn test_annotation_nhx_rejects_reserved_characters(#[case] (key, value): (&str, &str), #[case] expected: &str) {
    let attrs = BTreeMap::from([(key.to_owned(), NewickValue::String(value.to_owned()))]);

    assert_eq!(expected, format!("{:#}", write_nhx(&attrs).unwrap_err()));
  }

  #[test]
  fn test_annotation_beast_rejects_empty_key() {
    let attrs = BTreeMap::from([(String::new(), NewickValue::Boolean(true))]);

    assert_eq!(
      "A BEAST annotation needs a non-empty key",
      format!("{:#}", write_beast(&attrs).unwrap_err())
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::valid(          "[note [nested]]",  None)]
  #[case::no_brackets(    "note",             Some(r#"The raw comment "note" is not a Newick comment: it must start with '[', end with the matching ']', and close every quote"#))]
  #[case::unbalanced(     "[a]b]",            Some(r#"The raw comment "[a]b]" is not a Newick comment: it must start with '[', end with the matching ']', and close every quote"#))]
  #[case::open_quote(     "[say \"]",         Some(r#"The raw comment "[say \"]" is not a Newick comment: it must start with '[', end with the matching ']', and close every quote"#))]
  #[trace]
  fn test_annotation_raw_comment_validation(#[case] comment: &str, #[case] expected: Option<&str>) {
    let mut buffer = Vec::new();

    let actual = write_raw_comments(&mut buffer, &[comment.to_owned()]).err().map(|error| format!("{error:#}"));

    assert_eq!(expected.map(str::to_owned), actual);
  }

  #[test]
  fn test_annotation_writers_skip_empty_attrs() {
    let mut buffer = Vec::new();

    write_beast_attrs(&mut buffer, []).unwrap();
    write_nhx_attrs(&mut buffer, []).unwrap();

    assert!(buffer.is_empty());
  }

  mod helpers {
    use crate::annotation::{classify_comment, write_beast_attrs, write_nhx_attrs};
    use crate::types::NewickValue;
    use eyre::Report;
    use std::collections::BTreeMap;

    pub(super) fn parse_node_attrs(comment: &str) -> BTreeMap<String, NewickValue> {
      let mut attrs = BTreeMap::new();
      let mut raw = Vec::new();
      classify_comment(comment, &mut attrs, &mut raw);
      assert_eq!(Vec::<String>::new(), raw);
      attrs
    }

    pub(super) fn write_beast(attrs: &BTreeMap<String, NewickValue>) -> Result<String, Report> {
      let mut buffer = Vec::new();
      write_beast_attrs(&mut buffer, attrs.iter().map(|(key, value)| (key.as_str(), value)))?;
      Ok(String::from_utf8(buffer)?)
    }

    pub(super) fn write_nhx(attrs: &BTreeMap<String, NewickValue>) -> Result<String, Report> {
      let mut buffer = Vec::new();
      write_nhx_attrs(&mut buffer, attrs.iter().map(|(key, value)| (key.as_str(), value)))?;
      Ok(String::from_utf8(buffer)?)
    }
  }
}
