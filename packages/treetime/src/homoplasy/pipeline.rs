use crate::error::OperationError;

/// Scientific and algorithm policy for the homoplasy operation.
///
/// No core operation computes homoplasy yet, so the params carry no fields.
#[derive(Debug)]
pub struct HomoplasyParams;

/// Parsed domain input for the homoplasy operation.
///
/// No core operation computes homoplasy yet, so the input carries no fields.
#[derive(Debug)]
pub struct HomoplasyInput;

/// Aggregate result of the homoplasy operation.
///
/// No core operation computes homoplasy yet, so the output carries no fields.
#[derive(Debug)]
pub struct HomoplasyOutput;

/// Run homoplasy inference.
///
/// No core operation computes homoplasy yet, so `run` returns the typed
/// [`OperationError::NotImplemented`] failure instead of a result.
pub fn run(_params: &HomoplasyParams, _input: HomoplasyInput) -> Result<HomoplasyOutput, OperationError> {
  Err(OperationError::NotImplemented("homoplasy"))
}
