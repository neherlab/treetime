// The destruction site: a typed error collapsed to a string. Fns are private
// so stringly_error (public signatures only) stays out of this file's output.

use std::fmt;

#[derive(Debug)]
enum ParseError {
    Bad,
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "bad input")
    }
}

#[derive(Debug)]
enum Wrapped {
    Parse(ParseError),
}

fn to_string_is_flagged(r: Result<u32, ParseError>) -> Result<u32, String> {
    r.map_err(|e| e.to_string())
}

fn format_is_flagged(r: Result<u32, ParseError>) -> Result<u32, String> {
    r.map_err(|e| format!("parse failed: {e}"))
}

fn already_string_is_fine(r: Result<u32, String>) -> Result<u32, String> {
    r.map_err(|e| e.to_string())
}

fn wrapping_is_fine(r: Result<u32, ParseError>) -> Result<u32, Wrapped> {
    r.map_err(Wrapped::Parse)
}

fn main() {
    let _ = to_string_is_flagged(Err(ParseError::Bad));
    let _ = format_is_flagged(Err(ParseError::Bad));
    let _ = already_string_is_fine(Err("x".to_owned()));
    let _ = wrapping_is_fine(Err(ParseError::Bad));
}
