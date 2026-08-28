use itertools::Itertools;
use strsim::levenshtein;

/// Closest candidate to `input` for a "did you mean?" hint, or `None` when nothing is close enough.
///
/// Uses Levenshtein edit distance with an input-length-proportional threshold so short tokens must
/// match tightly and long tokens tolerate more typos. This is the same shape of heuristic clap uses
/// for unknown-flag suggestions; it keeps unrelated candidates out of the hint.
pub fn did_you_mean(input: &str, candidates: &[&str]) -> Option<String> {
  let threshold = (input.len() / 3).max(1) + 1;
  candidates
    .iter()
    .map(|candidate| (levenshtein(input, candidate), *candidate))
    .filter(|(distance, _)| *distance <= threshold)
    .min_by_key(|(distance, candidate)| (*distance, candidate.len()))
    .map(|(_, candidate)| candidate.to_owned())
}

/// Render a sorted, comma-separated list of valid values for an error message.
pub fn valid_values(candidates: &[&str]) -> String {
  candidates.iter().sorted().map(|candidate| format!("`{candidate}`")).join(", ")
}

/// Compose the trailing hint of an error message: an optional "did you mean" plus the valid set.
pub fn suggestion_suffix(input: &str, candidates: &[&str]) -> String {
  match did_you_mean(input, candidates) {
    Some(best) => format!("did you mean `{best}`? Valid values: {}", valid_values(candidates)),
    None => format!("valid values: {}", valid_values(candidates)),
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_suggest_did_you_mean_close_typo_matches() {
    assert_eq!(Some("timetree".to_owned()), did_you_mean("timtree", &["timetree", "ancestral", "clock"]));
  }

  #[test]
  fn test_suggest_did_you_mean_unrelated_returns_none() {
    assert_eq!(None, did_you_mean("xyzzy", &["timetree", "ancestral", "clock"]));
  }

  #[test]
  fn test_suggest_did_you_mean_single_char_typo_short_token() {
    assert_eq!(Some("clock".to_owned()), did_you_mean("clok", &["timetree", "clock", "prune"]));
  }

  #[test]
  fn test_suggest_valid_values_sorted_and_quoted() {
    assert_eq!("`ancestral`, `clock`, `timetree`", valid_values(&["timetree", "ancestral", "clock"]));
  }
}
