use crate::dialect::NewickDialect;
use crate::grammar::{Rule, matches};
use crate::model::comment::{MrBayesComment, NewickComment};
use crate::model::value::NewickValue;
use crate::nhx::{DUPLICATION_VALUES, NhxType};
use crate::number::format_shortest;
use crate::write::conversions::{Conversion, DataKind, conversion};
use eyre::{Report, eyre};
use std::fmt::Write;
use std::{iter, slice};

pub(crate) fn encode_comment(comment: &NewickComment, dialect: NewickDialect) -> Result<Option<String>, Report> {
  let data = match comment {
    NewickComment::Beast(_) => DataKind::BeastComments,
    NewickComment::Nhx(_) => DataKind::NhxComments,
    NewickComment::MrBayesMcmc(_) => DataKind::MrBayesComments,
    NewickComment::Plain(text) if text.starts_with('&') => DataKind::AmpersandPlainComments,
    NewickComment::Plain(text) => return encode_plain(text).map(Some),
  };
  match conversion(dialect, data) {
    Conversion::Drop => Ok(None),
    Conversion::Fail => Err(eyre!(
      "The comment {comment:?} cannot be written in the {dialect} dialect, which would read it as an annotation"
    )),
    Conversion::Keep => match comment {
      NewickComment::Beast(pairs) => encode_beast(pairs).map(Some),
      NewickComment::Nhx(pairs) => encode_nhx(pairs).map(Some),
      NewickComment::MrBayesMcmc(mrbayes) => encode_mrbayes(mrbayes).map(Some),
      NewickComment::Plain(text) => encode_plain(text).map(Some),
    },
  }
}

pub(crate) fn encode_beast(pairs: &[(String, NewickValue)]) -> Result<String, Report> {
  let mut text = String::from("[&");
  for (i, (key, value)) in pairs.iter().enumerate() {
    if i > 0 {
      text.push(',');
    }
    if matches(Rule::beast_bare_key_exact, key) {
      text.push_str(key);
    } else {
      push_double_quoted(&mut text, key);
    }
    text.push('=');
    push_beast_value(&mut text, value).map_err(|error| eyre!("In the BEAST annotation {key:?}: {error}"))?;
  }
  text.push(']');
  Ok(text)
}

pub(crate) fn encode_nhx(pairs: &[(String, NewickValue)]) -> Result<String, Report> {
  let mut text = String::from("[&&NHX");
  for (key, value) in pairs {
    if !matches(Rule::nhx_key_exact, key) {
      return Err(eyre!(
        "NHX cannot hold the tag {key:?}: a tag is not empty and contains none of ':', '=', '[', ']', '>'"
      ));
    }
    text.push(':');
    text.push_str(key);
    if let Some(value) = nhx_value(key, value)? {
      text.push('=');
      text.push_str(&value);
    }
  }
  text.push(']');
  Ok(text)
}

fn encode_mrbayes(comment: &MrBayesComment) -> Result<String, Report> {
  let mut text = format!("[&{} ", comment.kind.letter());
  for (i, token) in iter::once(&comment.name).chain(&comment.values).enumerate() {
    if !matches(Rule::mrbayes_token_exact, token) {
      return Err(eyre!(
        "The MrBayes comment {:?} cannot hold {token:?}: a token is not empty and contains no whitespace, '[' or ']'",
        comment.name
      ));
    }
    if i > 0 {
      text.push(' ');
    }
    text.push_str(token);
  }
  text.push(']');
  Ok(text)
}

fn encode_plain(text: &str) -> Result<String, Report> {
  let comment = format!("[{text}]");
  if !matches(Rule::plain_comment_exact, &comment) {
    return Err(eyre!(
      "The comment text {text:?} cannot be written as a comment, because its brackets do not balance"
    ));
  }
  Ok(comment)
}

fn push_beast_value(text: &mut String, value: &NewickValue) -> Result<(), Report> {
  let mut pending: Vec<(slice::Iter<'_, NewickValue>, bool)> = Vec::new();
  let mut current = Some(value);
  loop {
    if let Some(value) = current.take() {
      match value {
        NewickValue::Array(values) => {
          text.push('{');
          pending.push((values.as_slice().iter(), true));
        },
        NewickValue::Boolean(true) => text.push_str("TRUE"),
        NewickValue::Boolean(false) => text.push_str("FALSE"),
        NewickValue::Number(number) => text.push_str(&format_shortest(*number)?),
        NewickValue::NumberText(number) => text.push_str(checked_number_text(number)?),
        NewickValue::String(string) => push_double_quoted(text, string),
        NewickValue::Color([red, green, blue]) => write!(text, "#{red:02x}{green:02x}{blue:02x}")?,
      }
    }
    let Some((elements, is_first)) = pending.last_mut() else {
      return Ok(());
    };
    if let Some(next) = elements.next() {
      if !*is_first {
        text.push(',');
      }
      *is_first = false;
      current = Some(next);
    } else {
      text.push('}');
      pending.pop();
    }
  }
}

fn push_double_quoted(text: &mut String, value: &str) {
  text.push('"');
  text.push_str(&value.replace('"', "\"\""));
  text.push('"');
}

fn checked_number_text(text: &str) -> Result<&str, Report> {
  if matches(Rule::number_exact, text) {
    Ok(text)
  } else {
    Err(eyre!(
      "The value {text:?} is marked as a number, but it is not a number"
    ))
  }
}

fn nhx_value(key: &str, value: &NewickValue) -> Result<Option<String>, Report> {
  let expected = NhxType::of(key);
  let text = match (value, expected) {
    (NewickValue::Boolean(true), _) => return Ok(None),
    (NewickValue::Number(number), NhxType::Decimal | NhxType::Integer | NhxType::Text) => format_shortest(*number)?,
    (NewickValue::NumberText(number), NhxType::Decimal | NhxType::Integer | NhxType::Text) => {
      checked_number_text(number)?.to_owned()
    },
    (NewickValue::String(string), NhxType::Text) => nhx_part(key, string)?.to_owned(),
    (NewickValue::String(string), NhxType::Duplication) if DUPLICATION_VALUES.contains(&string.as_str()) => {
      string.clone()
    },
    (NewickValue::Color([red, green, blue]), NhxType::Color) => format!("{red}.{green}.{blue}"),
    (NewickValue::Array(values), NhxType::Text) if values.as_slice().len() > 1 => nhx_parts(key, values.as_slice())?,
    _ => {
      return Err(eyre!(
        "The NHX tag {key} needs {}, and NHX cannot hold the value {value:?} there",
        expected.description()
      ));
    },
  };
  if expected == NhxType::Integer && !matches(Rule::integer_exact, &text) {
    return Err(eyre!("The NHX tag {key} needs an integer, but the value is {text}"));
  }
  Ok(Some(text))
}

fn nhx_parts(key: &str, values: &[NewickValue]) -> Result<String, Report> {
  let parts = values
    .iter()
    .map(|value| match value {
      NewickValue::String(string) => Ok(nhx_part(key, string)?.to_owned()),
      NewickValue::Number(number) => format_shortest(*number),
      NewickValue::NumberText(number) => Ok(checked_number_text(number)?.to_owned()),
      NewickValue::Boolean(_) | NewickValue::Color(_) | NewickValue::Array(_) => {
        Err(eyre!("The NHX tag {key} cannot hold {value:?} as a part of its value"))
      },
    })
    .collect::<Result<Vec<_>, _>>()?;
  Ok(parts.join(">"))
}

fn nhx_part<'s>(key: &str, text: &'s str) -> Result<&'s str, Report> {
  if matches(Rule::nhx_part_exact, text) {
    Ok(text)
  } else {
    Err(eyre!(
      "NHX cannot hold the value {text:?} of the tag {key}: a value contains none of ':', '=', '[', ']', '>'"
    ))
  }
}
