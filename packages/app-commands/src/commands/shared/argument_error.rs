use eyre::Report;

pub fn argument_error(settings: &[&str], message: impl Into<String>) -> Report {
  Report::new(ArgumentError {
    message: message.into(),
    settings: settings.iter().map(|setting| (*setting).to_owned()).collect(),
  })
}

#[derive(Debug, derive_more::Display, derive_more::Error)]
#[display("{message}")]
pub struct ArgumentError {
  pub message: String,
  pub settings: Vec<String>,
}
