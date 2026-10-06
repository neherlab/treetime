use crate::grammar::{Rule, matches, parse};
use crate::model::comment::{MrBayesComment, MrBayesKind, NewickComment};
use crate::model::value::{NewickArray, NewickValue};
use crate::nhx::{DUPLICATION_VALUES, NhxType};
use crate::number::format_shortest;
use crate::read::context::MapContext;
use crate::read::error::{NewickError, NewickErrorKind};
use pest::iterators::Pair;

pub(crate) fn read_comment(
  pair: Pair<'_, Rule>,
  context: &mut MapContext<'_, '_>,
) -> Result<NewickComment, NewickError> {
  match pair.as_rule() {
    Rule::plain_comment => Ok(NewickComment::Plain(comment_body(pair.as_str()).to_owned())),
    Rule::malformed_annotation => {
      context.tolerate(
        NewickErrorKind::Annotation,
        &pair,
        format!(
          "The annotation {} does not follow the annotation syntax of the dialect",
          pair.as_str()
        ),
      )?;
      Ok(NewickComment::Plain(comment_body(pair.as_str()).to_owned()))
    },
    Rule::beast_comment => Ok(NewickComment::Beast(
      pair
        .into_inner()
        .map(|beast_pair| read_beast_pair(beast_pair, context))
        .collect::<Result<_, _>>()?,
    )),
    Rule::nhx_comment => Ok(NewickComment::Nhx(
      pair
        .into_inner()
        .map(|tag| read_nhx_tag(tag, context))
        .collect::<Result<_, _>>()?,
    )),
    Rule::mrbayes_comment => Ok(NewickComment::MrBayesMcmc(read_mrbayes(pair))),
    _ => Err(context.error(
      NewickErrorKind::Structure,
      &pair,
      format!("Unexpected token {:?} where a comment was expected", pair.as_str()),
    )),
  }
}

pub(crate) fn is_comment_rule(rule: Rule) -> bool {
  matches!(
    rule,
    Rule::plain_comment | Rule::malformed_annotation | Rule::beast_comment | Rule::nhx_comment | Rule::mrbayes_comment
  )
}

pub(crate) fn number_value(text: &str) -> NewickValue {
  match text.parse::<f64>() {
    Ok(number) if format_shortest(number).is_ok_and(|shortest| shortest == text) => NewickValue::Number(number),
    Ok(_) | Err(_) => NewickValue::NumberText(text.to_owned()),
  }
}

pub(crate) fn unquote(text: &str, quote: char) -> String {
  let doubled: String = [quote, quote].iter().collect();
  text
    .strip_prefix(quote)
    .and_then(|inner| inner.strip_suffix(quote))
    .unwrap_or(text)
    .replace(&doubled, &quote.to_string())
}

fn comment_body(text: &str) -> &str {
  text
    .strip_prefix('[')
    .and_then(|inner| inner.strip_suffix(']'))
    .unwrap_or(text)
}

fn read_beast_pair(pair: Pair<'_, Rule>, context: &MapContext<'_, '_>) -> Result<(String, NewickValue), NewickError> {
  let position = pair.clone();
  let mut tokens = pair.into_inner();
  let key = match tokens.next() {
    Some(key) => beast_text(&key),
    None => return Err(context.error(NewickErrorKind::Structure, &position, "An annotation pair has no key")),
  };
  let mut open: Vec<NewickArray> = Vec::new();
  let mut value = None;
  for token in tokens {
    let item = match token.as_rule() {
      Rule::array_open => {
        open.push(NewickArray::default());
        continue;
      },
      Rule::array_close => match open.pop() {
        Some(array) => NewickValue::Array(array),
        None => return Err(context.error(NewickErrorKind::Structure, &token, "An array closes without opening")),
      },
      Rule::color => match parse_hex_color(&token) {
        Some(color) => NewickValue::Color(color),
        None => return Err(context.error(NewickErrorKind::Structure, &token, "A color has no hex channels")),
      },
      Rule::boolean => NewickValue::Boolean(token.as_str().eq_ignore_ascii_case("true")),
      Rule::number => number_value(token.as_str()),
      _ => NewickValue::String(beast_text(&token)),
    };
    match open.last_mut() {
      Some(array) => array.push(item),
      None => value = Some(item),
    }
  }
  Ok((key, value.unwrap_or(NewickValue::Boolean(true))))
}

fn beast_text(token: &Pair<'_, Rule>) -> String {
  match token.as_rule() {
    Rule::double_quoted => unquote(token.as_str(), '"'),
    Rule::single_quoted => unquote(token.as_str(), '\''),
    _ => token.as_str().to_owned(),
  }
}

fn parse_hex_color(token: &Pair<'_, Rule>) -> Option<[u8; 3]> {
  let text = token.as_str();
  let mut color = [0_u8; 3];
  for (channel, start) in color.iter_mut().zip([1, 3, 5]) {
    *channel = match text.get(start..start + 2).map(|hex| u8::from_str_radix(hex, 16)) {
      Some(Ok(value)) => value,
      Some(Err(_)) | None => return None,
    };
  }
  Some(color)
}

fn read_nhx_tag(tag: Pair<'_, Rule>, context: &mut MapContext<'_, '_>) -> Result<(String, NewickValue), NewickError> {
  let position = tag.clone();
  let mut tokens = tag.into_inner();
  let key = tokens.next().map(|key| key.as_str().to_owned()).unwrap_or_default();
  let Some(value) = tokens.next() else {
    return Ok((key, NewickValue::Boolean(true)));
  };
  let parts: Vec<&str> = value.into_inner().map(|part| part.as_str()).collect();
  let expected = NhxType::of(&key);
  let typed = match parts.as_slice() {
    [single] => nhx_scalar(single, expected),
    _ if expected == NhxType::Text => Some(NewickValue::Array(
      parts
        .iter()
        .map(|part| NewickValue::String((*part).to_owned()))
        .collect::<Vec<_>>()
        .into(),
    )),
    _ => None,
  };
  if let Some(value) = typed {
    Ok((key, value))
  } else {
    let text = parts.join(">");
    context.tolerate(
      NewickErrorKind::Annotation,
      &position,
      format!(
        "The NHX tag {key} needs {}, but its value is {text:?}",
        expected.description()
      ),
    )?;
    Ok((key, NewickValue::String(text)))
  }
}

fn nhx_scalar(text: &str, expected: NhxType) -> Option<NewickValue> {
  match expected {
    NhxType::Text => Some(NewickValue::String(text.to_owned())),
    NhxType::Decimal => matches(Rule::number_exact, text).then(|| number_value(text)),
    NhxType::Integer => matches(Rule::integer_exact, text).then(|| number_value(text)),
    NhxType::Duplication => DUPLICATION_VALUES
      .contains(&text)
      .then(|| NewickValue::String(text.to_owned())),
    NhxType::Color => nhx_color(text),
  }
}

fn nhx_color(text: &str) -> Option<NewickValue> {
  let channels = match parse(Rule::nhx_color_exact, text) {
    Ok(pairs) => pairs.flatten().filter(|pair| pair.as_rule() == Rule::color_channel),
    Err(_) => return None,
  };
  let mut color = [0_u8; 3];
  for (slot, channel) in color.iter_mut().zip(channels) {
    *slot = match channel.as_str().parse::<u8>() {
      Ok(value) => value,
      Err(_) => return None,
    };
  }
  Some(NewickValue::Color(color))
}

fn read_mrbayes(pair: Pair<'_, Rule>) -> MrBayesComment {
  let mut tokens = pair.into_inner();
  let kind = match tokens.next().map(|kind| kind.as_str()) {
    Some("E") => MrBayesKind::E,
    Some("B") => MrBayesKind::B,
    _ => MrBayesKind::N,
  };
  let name = tokens.next().map(|name| name.as_str().to_owned()).unwrap_or_default();
  let values = tokens.map(|value| value.as_str().to_owned()).collect();
  MrBayesComment { kind, name, values }
}
