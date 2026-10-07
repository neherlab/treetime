use crate::grammar::Rule;
use crate::model::data::NewickHybrid;
use crate::read::comments::unquote;
use crate::read::context::MapContext;
use crate::read::error::{NewickError, NewickErrorKind};
use pest::iterators::Pair;

pub(crate) struct LabelToken {
  pub(crate) text: String,
  pub(crate) name: Option<String>,
  pub(crate) support: Option<Vec<f64>>,
  pub(crate) hybrid: Option<NewickHybrid>,
  pub(crate) is_acceptor: bool,
}

pub(crate) fn read_label(label: Pair<'_, Rule>, context: &MapContext<'_, '_>) -> Result<LabelToken, NewickError> {
  let mut token = LabelToken {
    text: label.as_str().to_owned(),
    name: None,
    support: None,
    hybrid: None,
    is_acceptor: false,
  };
  for part in label.into_inner() {
    match part.as_rule() {
      Rule::support_label => {
        let values = part
          .clone()
          .into_inner()
          .map(|number| {
            let value = number.as_str().parse::<f64>().map_err(|error| {
              context.error(
                NewickErrorKind::Syntax,
                &number,
                format!("The support value {:?} is not a number: {error}", number.as_str()),
              )
            })?;
            if value.is_finite() {
              Ok(value)
            } else {
              Err(context.error(
                NewickErrorKind::Syntax,
                &number,
                format!(
                  "The support value {:?} is too large for a 64-bit floating-point number",
                  number.as_str()
                ),
              ))
            }
          })
          .collect::<Result<Vec<_>, _>>()?;
        token.support = Some(values);
        token.name = Some(part.as_str().to_owned());
      },
      Rule::quoted_label => token.name = Some(unquote(part.as_str(), '\'')),
      Rule::hybrid_tag => {
        let (hybrid, is_acceptor) = read_hybrid_tag(&part, context)?;
        token.hybrid = Some(hybrid);
        token.is_acceptor = is_acceptor;
      },
      _ => {
        let text = part.as_str();
        token.name = Some(if context.options.underscores_as_spaces {
          text.replace('_', " ")
        } else {
          text.to_owned()
        });
      },
    }
  }
  Ok(token)
}

fn read_hybrid_tag(tag: &Pair<'_, Rule>, context: &MapContext<'_, '_>) -> Result<(NewickHybrid, bool), NewickError> {
  let mut hybrid = NewickHybrid::new(None, 0);
  let mut is_acceptor = false;
  for part in tag.clone().into_inner() {
    match part.as_rule() {
      Rule::hybrid_marker => is_acceptor = part.as_str() == "##",
      Rule::hybrid_kind => hybrid.kind = Some(part.as_str().to_owned()),
      _ => {
        hybrid.index = part.as_str().parse::<u32>().map_err(|error| {
          context.error(
            NewickErrorKind::Structure,
            tag,
            format!(
              "The hybrid node index in {:?} is not a valid index: {error}",
              tag.as_str()
            ),
          )
        })?;
      },
    }
  }
  Ok((hybrid, is_acceptor))
}
