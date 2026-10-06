use crate::number::{format_shortest, parse_number_text};
use crate::parse::is_comment_token;
use crate::types::NewickValue;
use eyre::{Report, eyre};
use std::collections::BTreeMap;
use std::io;

const NHX_PREFIX: &str = "&&NHX";
const NHX_RESERVED: [char; 5] = [':', '=', '[', ']', '"'];

pub(crate) fn classify_comment(comment: &str, attrs: &mut BTreeMap<String, NewickValue>, raw: &mut Vec<String>) {
  let inner = comment
    .strip_prefix('[')
    .and_then(|text| text.strip_suffix(']'))
    .unwrap_or(comment);
  let nhx_body = inner
    .get(..NHX_PREFIX.len())
    .filter(|prefix| prefix.eq_ignore_ascii_case(NHX_PREFIX))
    .and_then(|_| inner.get(NHX_PREFIX.len()..));
  if let Some(body) = nhx_body {
    parse_nhx_attrs(body, attrs);
  } else if let Some(body) = inner.strip_prefix('&') {
    parse_beast_attrs(body, attrs);
  } else {
    raw.push(comment.to_owned());
  }
}

pub fn write_beast_attrs<'a>(
  writer: &mut impl io::Write,
  attrs: impl IntoIterator<Item = (&'a str, &'a NewickValue)>,
) -> Result<(), Report> {
  let mut attrs = attrs.into_iter().peekable();
  if attrs.peek().is_none() {
    return Ok(());
  }
  let mut text = String::from("[&");
  for (i, (key, value)) in attrs.enumerate() {
    if i > 0 {
      text.push(',');
    }
    if key.is_empty() {
      return Err(eyre!("A BEAST annotation needs a non-empty key"));
    }
    if needs_beast_quoting(key) {
      push_double_quoted(&mut text, key);
    } else {
      text.push_str(key);
    }
    text.push('=');
    push_beast_value(&mut text, value)?;
  }
  text.push(']');
  writer.write_all(text.as_bytes())?;
  Ok(())
}

pub fn write_nhx_attrs<'a>(
  writer: &mut impl io::Write,
  attrs: impl IntoIterator<Item = (&'a str, &'a NewickValue)>,
) -> Result<(), Report> {
  let mut attrs = attrs.into_iter().peekable();
  if attrs.peek().is_none() {
    return Ok(());
  }
  let mut text = String::from("[");
  text.push_str(NHX_PREFIX);
  for (key, value) in attrs {
    if key.is_empty() {
      return Err(eyre!("An NHX annotation needs a non-empty key"));
    }
    check_nhx_text("key", key)?;
    text.push(':');
    text.push_str(key);
    text.push('=');
    push_nhx_value(&mut text, value)?;
  }
  text.push(']');
  writer.write_all(text.as_bytes())?;
  Ok(())
}

pub(crate) fn write_raw_comments(writer: &mut impl io::Write, comments: &[String]) -> Result<(), Report> {
  for comment in comments {
    if !is_comment_token(comment) {
      return Err(eyre!(
        "The raw comment {comment:?} is not a Newick comment: it must start with '[', end with the matching ']', and close every quote"
      ));
    }
    writer.write_all(comment.as_bytes())?;
  }
  Ok(())
}

fn push_beast_value(text: &mut String, value: &NewickValue) -> Result<(), Report> {
  match value {
    NewickValue::Boolean(true) => text.push_str("TRUE"),
    NewickValue::Boolean(false) => text.push_str("FALSE"),
    NewickValue::Number(number) => text.push_str(&format_shortest(*number)?),
    NewickValue::NumberText(number) => text.push_str(checked_number_text(number)?),
    NewickValue::String(string) => push_double_quoted(text, string),
    NewickValue::Array(values) => {
      text.push('{');
      for (i, value) in values.iter().enumerate() {
        if i > 0 {
          text.push(',');
        }
        push_beast_value(text, value)?;
      }
      text.push('}');
    },
  }
  Ok(())
}

fn push_nhx_value(text: &mut String, value: &NewickValue) -> Result<(), Report> {
  match value {
    NewickValue::Boolean(boolean) => text.push_str(if *boolean { "true" } else { "false" }),
    NewickValue::Number(number) => text.push_str(&format_shortest(*number)?),
    NewickValue::NumberText(number) => text.push_str(checked_number_text(number)?),
    NewickValue::String(string) => {
      check_nhx_text("value", string)?;
      text.push_str(string);
    },
    NewickValue::Array(values) => {
      for (i, value) in values.iter().enumerate() {
        if i > 0 {
          text.push('>');
        }
        push_nhx_value(text, value)?;
      }
    },
  }
  Ok(())
}

fn check_nhx_text(what: &str, text: &str) -> Result<(), Report> {
  if text.contains(NHX_RESERVED) {
    return Err(eyre!(
      "NHX cannot represent a {what} containing a reserved character (':', '=', '[', ']' or '\"'): {text}"
    ));
  }
  Ok(())
}

fn checked_number_text(text: &str) -> Result<&str, Report> {
  match parse_number_text(text) {
    Some(_) => Ok(text),
    None => Err(eyre!(
      "The annotation value {text:?} is marked as a number, but it is not a finite number"
    )),
  }
}

fn needs_beast_quoting(key: &str) -> bool {
  key.contains(|c: char| matches!(c, ',' | '=' | '[' | ']' | '{' | '}' | '\'' | '"' | '&') || c.is_whitespace())
    || key.eq_ignore_ascii_case("true")
    || key.eq_ignore_ascii_case("false")
}

fn push_double_quoted(text: &mut String, value: &str) {
  text.push('"');
  text.push_str(&value.replace('"', "\"\""));
  text.push('"');
}

fn parse_nhx_attrs(body: &str, attrs: &mut BTreeMap<String, NewickValue>) {
  for tag in body.split(':').map(str::trim).filter(|tag| !tag.is_empty()) {
    match tag.split_once('=') {
      Some((key, value)) => attrs.insert(key.trim().to_owned(), NewickValue::String(value.trim().to_owned())),
      None => attrs.insert(tag.to_owned(), NewickValue::Boolean(true)),
    };
  }
}

fn parse_beast_attrs(body: &str, attrs: &mut BTreeMap<String, NewickValue>) {
  for pair in split_top_level(body, ',') {
    let pair = pair.trim();
    if pair.is_empty() {
      continue;
    }
    match split_once_unquoted(pair, '=') {
      Some((key, value)) => attrs.insert(unquote(key.trim()), parse_beast_value(value.trim())),
      None => attrs.insert(unquote(pair), NewickValue::Boolean(true)),
    };
  }
}

fn parse_beast_value(text: &str) -> NewickValue {
  if let Some(inner) = text.strip_prefix('{').and_then(|text| text.strip_suffix('}')) {
    if inner.trim().is_empty() {
      return NewickValue::Array(Vec::new());
    }
    return NewickValue::Array(
      split_top_level(inner, ',')
        .into_iter()
        .map(|element| parse_beast_value(element.trim()))
        .collect(),
    );
  }
  if text.eq_ignore_ascii_case("true") {
    return NewickValue::Boolean(true);
  }
  if text.eq_ignore_ascii_case("false") {
    return NewickValue::Boolean(false);
  }
  match parse_number_text(text) {
    Some(number) if format_shortest(number).is_ok_and(|canonical| canonical == text) => NewickValue::Number(number),
    Some(_) => NewickValue::NumberText(text.to_owned()),
    None => NewickValue::String(unquote(text)),
  }
}

#[expect(
  clippy::string_slice,
  reason = "the indices come from char_indices of the same string"
)]
fn split_top_level(text: &str, separator: char) -> Vec<&str> {
  let mut parts = Vec::new();
  let mut depth = 0_usize;
  let mut start = 0;
  for (i, c) in unquoted_chars(text) {
    match c {
      '{' => depth += 1,
      '}' => depth = depth.saturating_sub(1),
      _ if c == separator && depth == 0 => {
        parts.push(&text[start..i]);
        start = i + c.len_utf8();
      },
      _ => {},
    }
  }
  parts.push(&text[start..]);
  parts
}

#[expect(
  clippy::string_slice,
  reason = "the index comes from char_indices of the same string"
)]
fn split_once_unquoted(text: &str, separator: char) -> Option<(&str, &str)> {
  let (i, c) = unquoted_chars(text).into_iter().find(|&(_, c)| c == separator)?;
  Some((&text[..i], &text[i + c.len_utf8()..]))
}

fn unquoted_chars(text: &str) -> Vec<(usize, char)> {
  let mut unquoted = Vec::with_capacity(text.len());
  let mut chars = text.char_indices().peekable();
  let mut open_quote = None;
  let mut at_token_start = true;
  while let Some((i, c)) = chars.next() {
    match open_quote {
      Some(quote) if c == quote => {
        if chars.next_if(|&(_, next)| next == quote).is_none() {
          open_quote = None;
        }
      },
      Some(_) => {},
      None if at_token_start && matches!(c, '"' | '\'') => {
        open_quote = Some(c);
        at_token_start = false;
      },
      None => {
        unquoted.push((i, c));
        if matches!(c, ',' | '=' | '{') {
          at_token_start = true;
        } else if !c.is_whitespace() {
          at_token_start = false;
        }
      },
    }
  }
  unquoted
}

fn unquote(text: &str) -> String {
  for quote in ['"', '\''] {
    if let Some(inner) = text.strip_prefix(quote).and_then(|text| text.strip_suffix(quote)) {
      return inner.replace(&format!("{quote}{quote}"), &quote.to_string());
    }
  }
  text.to_owned()
}
