use crate::error::OperationError;

#[derive(Debug)]
pub struct HomoplasyParams;

#[derive(Debug)]
pub struct HomoplasyInput;

#[derive(Debug)]
pub struct HomoplasyOutput;

pub fn run(_params: &HomoplasyParams, _input: HomoplasyInput) -> Result<HomoplasyOutput, OperationError> {
  Err(OperationError::NotImplemented("homoplasy"))
}
