use crate::alphabet::alphabet::{FILL_CHAR, NON_CHAR, ProfileMap, VARIABLE_CHAR};
use crate::make_error;
use color_eyre::{Section, SectionExt};
use eyre::{Report, WrapErr};
use indexmap::IndexMap;
use itertools::{Itertools, chain};
use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};
use std::iter::once;
use treetime_io::json::{JsonPretty, json_write_str};
use treetime_primitives::{AsciiChar, StateSet};
use treetime_utils::string::quote;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlphabetConfig {
  pub canonical: Vec<u8>,
  pub ambiguous: IndexMap<u8, Vec<u8>>,
  pub unknown: u8,
  pub gap: u8,
  pub treat_gap_as_unknown: bool,
}

impl AlphabetConfig {
  pub fn create_profile_map(&self) -> Result<ProfileMap, Report> {
    let AlphabetConfig {
      canonical,
      ambiguous,
      unknown,
      gap,
      treat_gap_as_unknown,
    } = self;

    let gap = AsciiChar::from(*gap);
    let unknown = AsciiChar::from(*unknown);

    self
      .validate()
      .wrap_err("When validating alphabet config")
      .with_section(|| {
        json_write_str(&self, JsonPretty(true))
          .unwrap()
          .header("Alphabet config")
      })?;

    let eye = Array2::<f64>::eye(canonical.len());

    // Add canonical to profile map
    let mut profile_map: ProfileMap = canonical
      .iter()
      .zip(eye.rows())
      .map(|(s, x)| (AsciiChar(*s), x.to_owned()))
      .collect();

    // Add unknown to profile map
    profile_map.insert(unknown, Array1::<f64>::ones(canonical.len()));

    // Add ambiguous to profile map
    ambiguous.iter().for_each(|(&key, values)| {
      let profile = canonical
        .iter()
        .enumerate()
        .map(|(i, c)| if values.contains(c) { 1.0 } else { 0.0 })
        .collect::<Array1<f64>>();
      profile_map.insert(AsciiChar(key), profile);
    });

    if *treat_gap_as_unknown {
      // Add gap to profile map
      profile_map.insert(gap, profile_map[&unknown].clone());
    }

    Ok(profile_map)
  }

  pub fn validate(&self) -> Result<(), Report> {
    let AlphabetConfig {
      canonical,
      ambiguous,
      unknown,
      gap,
      ..
    } = self;

    {
      let all = chain![
        self.canonical.iter().copied(),
        self.ambiguous.keys().copied(),
        self.ambiguous.values().flatten().copied(),
        once(self.unknown),
        once(self.gap)
      ]
      .collect_vec();

      for reserved in [NON_CHAR, VARIABLE_CHAR, FILL_CHAR] {
        if all.iter().any(|&c| c == u8::from(reserved)) {
          return make_error!("Alphabet contains reserved character: {reserved}");
        }
      }
    }

    {
      let canonical_dupes = canonical.iter().duplicates().collect_vec();
      if !canonical_dupes.is_empty() {
        let msg = canonical_dupes.into_iter().join(", ");
        return make_error!("Canonical set contains duplicates: {msg}");
      }
    }

    {
      let ambiguous_dupes = ambiguous.keys().duplicates().collect_vec();
      if !ambiguous_dupes.is_empty() {
        let msg = ambiguous_dupes.into_iter().join(", ");
        return make_error!("Ambiguous set contains duplicates: {msg}");
      }
    }

    let canonical: StateSet = canonical.iter().copied().collect();
    let ambiguous_keys: StateSet = ambiguous.keys().copied().collect();
    let ambiguous_set_map: IndexMap<u8, StateSet> = ambiguous
      .iter()
      .map(|(key, vals)| (*key, vals.iter().copied().collect()))
      .collect();
    {
      let canonical_inter_ambig: StateSet = canonical.intersection(&ambiguous_keys);
      if !canonical_inter_ambig.is_empty() {
        let msg = canonical_inter_ambig.iter().join(", ");
        return make_error!("Canonical and ambiguous sets must be disjoint, but these characters are shared: {msg}");
      }
    }

    if canonical.contains(*gap) {
      let msg = canonical.iter().map(quote).join(", ");
      return make_error!("Canonical set contains 'gap' character: {msg}");
    }

    if canonical.contains(*unknown) {
      let msg = canonical.iter().map(quote).join(", ");
      return make_error!("Canonical set contains 'unknown' character: {msg}");
    }

    if ambiguous_keys.contains(*gap) {
      let msg = ambiguous.keys().map(quote).join(", ");
      return make_error!("Ambiguous set contains 'gap' character: {msg}");
    }

    if ambiguous_keys.contains(*gap) {
      let msg = ambiguous.keys().map(quote).join(", ");
      return make_error!("Ambiguous set contains 'unknown' character: {msg}");
    }

    {
      let ambig_gaps = ambiguous_set_map
        .iter()
        .map(|(key, vals)| (key, vals.difference(&canonical)))
        .filter(|(key, extra)| !extra.is_empty())
        .collect_vec();

      if !ambig_gaps.is_empty() {
        let msg = ambig_gaps
          .iter()
          .map(|(key, vals)| format!("\"{key}\" => {}", vals.iter().map(quote).join(", ")))
          .join("; ");
        return make_error!("Ambiguous set can only contain mapping(s) to canonical characters, but found: {msg}");
      }
    }

    Ok(())
  }
}
