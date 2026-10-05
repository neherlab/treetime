use aide::generate::{in_context, on_error, reset_context};
use eyre::Report;
use itertools::Itertools;
use schemars::generate::SchemaSettings;
use std::cell::RefCell;
use std::rc::Rc;
use treetime_utils::make_error;

const COMPONENTS_PREFIX: &str = "#/components/schemas/";

pub(crate) fn with_project_schemas<R>(settings: &SchemaSettings, build: impl FnOnce() -> R) -> Result<R, Report> {
  reset_context();
  let errors = Rc::new(RefCell::new(Vec::<String>::new()));
  let collected = Rc::clone(&errors);
  on_error(move |error| collected.borrow_mut().push(error.to_string()));
  in_context(|ctx| {
    ctx.schema = settings
      .clone()
      .with(|settings| settings.definitions_path = COMPONENTS_PREFIX.into())
      .into_generator();
  });
  let built = build();
  reset_context();
  let errors = errors.borrow();
  if errors.is_empty() {
    Ok(built)
  } else {
    make_error!(
      "the OpenAPI document of the router is inconsistent: {}",
      errors.iter().join("; ")
    )
  }
}
