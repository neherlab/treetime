use deser::Serialize;
use schemars::Schema;
use schemars::transform::Transform;
use serde_json::Value;
use treetime_utils::io::json::to_json_value;

#[cfg(feature = "clap")]
use clap::Arg;

#[expect(
  clippy::expect_used,
  reason = "a schema transform cannot return errors, and a documented default is a constant that always serializes"
)]
pub fn documented_default<T: Serialize>(value: &T) -> DocumentedDefault {
  DocumentedDefault(to_json_value(value).expect("the documented default serializes to JSON"))
}

pub struct DocumentedDefault(Value);

impl DocumentedDefault {
  #[cfg(feature = "clap")]
  pub fn help(&self, arg: Arg) -> Arg {
    let spec = format!("[default: {}]", self.text());
    let help = arg.get_help().map(|help| format!("{help} {spec}"));
    let long_help = arg
      .get_long_help()
      .or_else(|| arg.get_help())
      .map(|help| format!("{help}\n\n{spec}"));
    let arg = match help {
      Some(help) => arg.help(help),
      None => arg.help(spec.clone()),
    };
    match long_help {
      Some(long_help) => arg.long_help(long_help),
      None => arg.long_help(spec),
    }
  }

  #[cfg(feature = "clap")]
  fn text(&self) -> String {
    match &self.0 {
      Value::String(text) => text.clone(),
      value => value.to_string(),
    }
  }
}

impl Transform for DocumentedDefault {
  fn transform(&mut self, schema: &mut Schema) {
    schema.insert("default".to_owned(), self.0.clone());
  }
}
