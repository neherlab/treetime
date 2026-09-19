#![allow(dead_code, unused_variables, clippy::all, reason = "ui fixture exercises literals that are never used")]

// FAILING: curly quotes, em dash, and emoji as literal characters in the source.
const CURLY: &str = "the ‘quoted’ word"; // should warn: curly quote
const DASH: &str = "range 1—9"; // should warn: em/en dash
const EMOJI: &str = "done 😀"; // should warn: emoji

// PASSING: plain ASCII punctuation and legitimate escapes.
const STRAIGHT: &str = "the 'quoted' word - range 1-9 done";
const NEWLINE_ESCAPE: &str = "line\tbreak\n";

fn main() {}
