// Public signatures with string error types are flagged. Private functions,
// trait-impl methods (signature dictated elsewhere), and main are not.

use std::borrow::Cow;
use std::str::FromStr;

pub fn public_string_err_is_flagged(x: u32) -> Result<u32, String> {
    if x > 0 { Ok(x) } else { Err("zero".to_owned()) }
}

pub fn public_str_err_is_flagged(x: u32) -> Result<u32, &'static str> {
    if x > 0 { Ok(x) } else { Err("zero") }
}

pub fn public_cow_err_is_flagged(x: u32) -> Result<u32, Cow<'static, str>> {
    if x > 0 { Ok(x) } else { Err(Cow::Borrowed("zero")) }
}

pub trait Parser {
    fn required_method_is_flagged(&self, input: &str) -> Result<u32, String>;
}

pub enum RealError {
    Zero,
}

pub fn real_error_is_fine(x: u32) -> Result<u32, RealError> {
    if x > 0 { Ok(x) } else { Err(RealError::Zero) }
}

fn private_string_err_is_fine(x: u32) -> Result<u32, String> {
    if x > 0 { Ok(x) } else { Err("zero".to_owned()) }
}

pub struct Numeric(u32);

impl FromStr for Numeric {
    type Err = String;

    // Foreign-trait impl: the signature is FromStr's, not ours.
    fn from_str(s: &str) -> Result<Self, String> {
        s.parse().map(Numeric).map_err(|e| e.to_string())
    }
}

fn main() -> Result<(), String> {
    let _ = public_string_err_is_flagged(1);
    let _ = public_str_err_is_flagged(1);
    let _ = public_cow_err_is_flagged(1);
    let _ = real_error_is_fine(1);
    let _ = private_string_err_is_fine(1);
    let _ = Numeric::from_str("1");
    Ok(())
}
