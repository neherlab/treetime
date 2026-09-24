use itertools::Itertools;
use strsim::levenshtein;

pub(crate) fn suggestion_suffix(input: &str, candidates: &[&str]) -> String {
  match did_you_mean(input, candidates) {
    Some(best) => format!("did you mean `{best}`? Valid values: {}", valid_values(candidates)),
    None => format!("valid values: {}", valid_values(candidates)),
  }
}

#[allow(clippy::integer_division, reason = "integer division is the intended floor division")]
fn did_you_mean(input: &str, candidates: &[&str]) -> Option<String> {
  let threshold = (input.len() / 3).max(1) + 1;
  candidates
    .iter()
    .map(|candidate| (levenshtein(input, candidate), *candidate))
    .filter(|(distance, _)| *distance <= threshold)
    .min_by_key(|(distance, candidate)| (*distance, candidate.len()))
    .map(|(_, candidate)| candidate.to_owned())
}

pub(crate) fn valid_values(candidates: &[&str]) -> String {
  candidates
    .iter()
    .sorted()
    .map(|candidate| format!("`{candidate}`"))
    .join(", ")
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_suggest_did_you_mean_close_typo_matches() {
    assert_eq!(
      Some("timetree".to_owned()),
      did_you_mean("timtree", &["timetree", "ancestral", "clock"])
    );
  }

  #[test]
  fn test_suggest_did_you_mean_unrelated_returns_none() {
    assert_eq!(None, did_you_mean("xyzzy", &["timetree", "ancestral", "clock"]));
  }

  #[test]
  fn test_suggest_did_you_mean_single_char_typo_short_token() {
    assert_eq!(
      Some("clock".to_owned()),
      did_you_mean("clok", &["timetree", "clock", "prune"])
    );
  }

  #[test]
  fn test_suggest_valid_values_sorted_and_quoted() {
    assert_eq!(
      "`ancestral`, `clock`, `timetree`",
      valid_values(&["timetree", "ancestral", "clock"])
    );
  }
}
