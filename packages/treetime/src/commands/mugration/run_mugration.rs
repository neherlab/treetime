use crate::commands::mugration::mugration_args::TreetimeMugrationArgs;
use crate::io::discrete_states_csv::read_discrete_attrs;
use crate::{make_error, make_internal_report};
use eyre::Report;
use indexmap::{IndexMap, IndexSet};
use itertools::Itertools;
use log::warn;
use ndarray::Array1;
use num_traits::ToPrimitive;
use statrs::statistics::Statistics;
use std::fmt::Display;

pub fn run_mugration(mugration_args: &TreetimeMugrationArgs) -> Result<(), Report> {
  let TreetimeMugrationArgs {
    tree,
    attribute,
    states,
    weights,
    name_column,
    confidence,
    pc,
    missing_data,
    missing_weights_threshold,
    sampling_bias_correction,
    outdir,
    seed,
  } = mugration_args;

  let (attr_values, attr_name) =
    read_discrete_attrs::<String>(states, name_column, &Some(attribute.clone()), |s| Ok(s.to_owned()))?;

  let unique_values: IndexSet<String> = attr_values.values().sorted().cloned().collect();

  let (weights, mut alphabet) = if let Some(weights_filepath) = weights {
    let weights = read_discrete_attrs::<f64>(weights_filepath, &Some(attr_name), &Some("weight".to_owned()), |s| {
      Ok(s.parse::<f64>()?)
    })?
    .0;

    let unique_values_in_weights: IndexSet<String> = weights.keys().sorted().cloned().collect();

    let unique_values: IndexSet<String> = unique_values.union(&unique_values_in_weights).cloned().collect();

    let missing_values_in_weights: IndexSet<String> = unique_values
      .difference(&unique_values_in_weights)
      .filter(|&value| value != missing_data)
      .cloned()
      .collect();

    let missing_weights_ratio = missing_values_in_weights.len() as f64 / unique_values.len() as f64;

    if !missing_values_in_weights.is_empty() {
      warn!(
        "Mugration: discrete attributes missing from weights file: {} (ratio: {missing_weights_ratio:.3})",
        missing_values_in_weights.iter().join(", ")
      );
    }

    if missing_weights_ratio > *missing_weights_threshold {
      return make_error!("Mugration: too many discrete attributes missing from the weights file. The ratio of missing values {missing_weights_ratio} is greater than the threshold {missing_weights_threshold}. Weights were read from file {weights_filepath:?}");
    }

    let alphabet = calculate_alphabet(&unique_values, missing_data)?;

    let mean_weight = weights.values().mean();

    let weights: Array1<f64> = alphabet
      .letter_to_attr
      .iter()
      .map(|(letter, value)| {
        let attr = &alphabet.letter_to_attr[letter];
        let weight = weights.get(attr).unwrap_or(&mean_weight);
        *weight
      })
      .collect();

    let sum = weights.sum();
    let weights = weights / sum;

    (Some(weights), alphabet)
  } else {
    let alphabet = calculate_alphabet(&unique_values, missing_data)?;
    (None, alphabet)
  };

  let num_attrs = alphabet.letters.len();
  if num_attrs < 2 {
    return make_error!(
      "Mugration: only {num_attrs} discrete attributes provided for mugration. This does not make sense."
    );
  }

  let missing_char = chr(65 + num_attrs)?;
  alphabet.attr_to_letter.insert(missing_data.clone(), missing_char);
  alphabet.letter_to_attr.insert(missing_char, missing_data.clone());

  // let gtr = GTR::new(GTRParams {
  //   alphabet: Alphabet::with_letters(&alphabet.letters, missing_char)?,
  //   mu: 1.0,
  //   W: Some(Array2::<f64>::ones((num_attrs, num_attrs))),
  //   pi: weights.unwrap_or_else(|| Array1::<f64>::ones(num_attrs)),
  // })?;

  Ok(())
}

pub fn chr<I: ToPrimitive + Display + Copy>(i: I) -> Result<char, Report> {
  let i_u32: u32 = i
    .to_u32()
    .ok_or_else(|| make_internal_report!("chr(): cannot convert integer {i} to a character"))?;
  let c: char = i_u32.try_into()?;
  Ok(c)
}

pub struct DiscreteAttrAlphabet {
  pub attr_to_letter: IndexMap<String, char>,
  pub letter_to_attr: IndexMap<char, String>,
  pub letters: Vec<char>,
}

pub fn calculate_alphabet(
  unique_values: &IndexSet<String>,
  missing_data: &str,
) -> Result<DiscreteAttrAlphabet, Report> {
  let attr_to_letter: IndexMap<String, char> = unique_values
    .iter()
    .filter(|&value| value != missing_data)
    .enumerate()
    .map(|(i, value)| {
      let letter = chr(i + 65)?;
      Ok((value.clone(), letter))
    })
    .collect::<Result<IndexMap<String, char>, Report>>()?;

  let letter_to_attr: IndexMap<char, String> = attr_to_letter
    .iter()
    .map(|(key, value)| (*value, key.clone()))
    .collect();

  let letters = attr_to_letter.values().copied().collect_vec();

  Ok(DiscreteAttrAlphabet {
    attr_to_letter,
    letter_to_attr,
    letters,
  })
}
