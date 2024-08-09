use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::seq_utils::normalize_profile::normalize_profile;
use crate::utils::ndarray::{argmax_axis, cumsum_axis, random};
use eyre::Report;
use ndarray::{stack, Array1, Array2, ArrayBase, Axis, Data, Ix1, Ix2};
use rand::Rng;

#[derive(Debug, Clone)]
pub struct Prof2SeqParams {
  pub should_sample_from_profile: bool,
  pub should_normalize_profile: bool,
}

#[derive(Debug, Clone)]
pub struct Prof2SeqResult {
  pub seq: Array1<char>,
  pub prof_values: Array1<f64>,
  pub seq_ii: Array1<usize>,
}

/// Convert profile to sequence and normalize profile across sites.
///
/// Parameters
/// ----------
///
///  profile : numpy 2D array
///     Profile. Shape of the profile should be (L x a), where L - sequence
///     length, a - alphabet size.
///  gtr : gtr.GTR
///     Instance of the GTR class to supply the sequence alphabet
///  collapse_prof : bool
///     Whether to convert the profile to the delta-function
///
/// Returns
/// -------
///  seq : numpy.array
///     Sequence as numpy array of length L
///  prof_values :  numpy.array
///     Values of the profile for the chosen sequence characters (length L)
///  idx : numpy.array
///     Indices chosen from profile as array of length L
pub fn prof2seq<S, R: Rng>(
  profile: &ArrayBase<S, Ix2>,
  gtr: &GTR,
  rng: &mut R,
  params: &Prof2SeqParams,
) -> Result<Prof2SeqResult, Report>
where
  S: Data<Elem = f64>,
{
  // Normalize profile such that probabilities at each site sum to one
  let profile = if params.should_normalize_profile {
    let (prof_norm, _) = normalize_profile(profile, false)?;
    prof_norm
  } else {
    profile.to_owned()
  };

  let seq_ii: Array1<usize> = if params.should_sample_from_profile {
    let randnum: Array1<f64> = random(profile.shape()[0], rng);
    sample_from_prof(&profile, &randnum)
  } else {
    argmax_axis(&profile, Axis(1))
  };

  let seq = gtr.alphabet.indices_to_sequence(seq_ii.iter().copied()).collect(); // max LH over the alphabet

  let prof_values = get_prof_values(&profile, &seq_ii);

  Ok(Prof2SeqResult {
    seq,
    prof_values,
    seq_ii,
  })
}

/// Sample sequence according to the probabilities in the profile
/// (sampling from cumulative distribution over the different states)
pub fn sample_from_prof(profile: &Array2<f64>, randnum: &Array1<f64>) -> Array1<usize> {
  assert_eq!(profile.shape()[0], randnum.shape()[0]);

  let cumdis: Array2<f64> = cumsum_axis(profile, Axis(1)).t().to_owned();

  cumdis
    .axis_iter(Axis(1))
    .enumerate()
    .map(|(i, row)| {
      for (j, val) in row.iter().enumerate() {
        if val > &randnum[i] {
          return j;
        }
      }
      row.len()
    })
    .collect()
}

pub fn get_prof_values(profile: &Array2<f64>, seq_ii: &Array1<usize>) -> Array1<f64> {
  profile
    .axis_iter(Axis(0))
    .enumerate()
    .map(|(i, row)| {
      let index = seq_ii[i];
      row[index]
    })
    .collect()
}

/// Convert the given character sequence into the profile according to the
/// alphabet specified.
///
/// Parameters
/// ----------
///
///  seq : numpy.array
///     Sequence to be converted to the profile
///
///  profile_map : dic
///     Mapping valid characters to profiles
///
/// Returns
/// -------
///
///  idx : numpy.array
///     Profile for the character. Zero array if the character not found
pub fn seq2prof<S>(seq: &ArrayBase<S, Ix1>, alphabet: &Alphabet) -> Result<Array2<f64>, Report>
where
  S: Data<Elem = char>,
{
  let prof = stack(
    Axis(0),
    seq.map(|&c| alphabet.get_profile(c).view()).as_slice().unwrap(),
  )?;
  Ok(prof)
}

#[cfg(test)]
mod tests {
  #![allow(clippy::excessive_precision, clippy::lossy_float_literal)]

  use super::*;
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::get_gtr::{jc69, JC69Params};
  use approx::assert_ulps_eq;
  use eyre::Report;
  use lazy_static::lazy_static;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rand::SeedableRng;
  use rand_isaac::Isaac64Rng;
  use rstest::rstest;

  lazy_static! {
    static ref ALPHABET: Alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
  }

  #[rstest]
  fn samples_from_profile() -> Result<(), Report> {
    let rng = &mut Isaac64Rng::seed_from_u64(42);

    let profile = array![
      [0.19356424, 0.25224431, 0.21259213, 0.19217803, 0.14942128],
      [0.19440831, 0.13170981, 0.26841564, 0.29005381, 0.11541244],
      [0.27439982, 0.18330691, 0.19687558, 0.32079767, 0.02462001],
      [0.03366488, 0.00781195, 0.32170632, 0.30066296, 0.33615390],
      [0.31185458, 0.25466645, 0.14705881, 0.24872985, 0.03769030],
      [0.24016971, 0.05380214, 0.35454510, 0.19585567, 0.15562739],
      [0.12705805, 0.37184099, 0.21907519, 0.27300161, 0.00902417]
    ];

    let randnum = array![0.6176355, 0.61209572, 0.616934, 0.94374808, 0.6818203, 0.3595079, 0.43703195];

    let sample: Array1<usize> = sample_from_prof(&profile, &randnum);

    assert_eq!(sample, array![2, 3, 2, 4, 2, 2, 1]);

    Ok(())
  }

  #[rstest]
  fn gets_prof_values() -> Result<(), Report> {
    let rng = &mut Isaac64Rng::seed_from_u64(42);

    let profile = array![
      [0.19356424, 0.25224431, 0.21259213, 0.19217803, 0.14942128],
      [0.19440831, 0.13170981, 0.26841564, 0.29005381, 0.11541244],
      [0.27439982, 0.18330691, 0.19687558, 0.32079767, 0.02462001],
      [0.03366488, 0.00781195, 0.32170632, 0.30066296, 0.33615390],
      [0.31185458, 0.25466645, 0.14705881, 0.24872985, 0.03769030],
      [0.24016971, 0.05380214, 0.35454510, 0.19585567, 0.15562739],
      [0.12705805, 0.37184099, 0.21907519, 0.27300161, 0.00902417],
    ];

    let seq_ii = array![2, 3, 2, 4, 2, 2, 1];

    let prof_values = get_prof_values(&profile, &seq_ii);

    assert_eq!(
      prof_values,
      array![0.21259213, 0.29005381, 0.19687558, 0.3361539, 0.14705881, 0.3545451, 0.37184099]
    );

    Ok(())
  }

  #[rstest]
  fn calculates_prof2seq_with_sample_without_normalize() -> Result<(), Report> {
    let rng = &mut Isaac64Rng::seed_from_u64(42);

    let gtr = jc69(JC69Params::default())?;

    let norm_prof = array![
      [0.19356424, 0.25224431, 0.21259213, 0.19217803, 0.14942128],
      [0.19440831, 0.13170981, 0.26841564, 0.29005381, 0.11541244],
      [0.27439982, 0.18330691, 0.19687558, 0.32079767, 0.02462001],
      [0.03366488, 0.00781195, 0.32170632, 0.30066296, 0.33615390],
      [0.31185458, 0.25466645, 0.14705881, 0.24872985, 0.03769030],
      [0.24016971, 0.05380214, 0.35454510, 0.19585567, 0.15562739],
      [0.12705805, 0.37184099, 0.21907519, 0.27300161, 0.00902417],
    ];

    let Prof2SeqResult {
      seq,
      prof_values,
      seq_ii,
    } = prof2seq(
      &norm_prof,
      &gtr,
      rng,
      &Prof2SeqParams {
        should_sample_from_profile: true,
        should_normalize_profile: false,
      },
    )?;

    assert_eq!(seq_ii, array![2, 3, 3, 2, 1, 3, 3]);

    assert_eq!(seq, array!['G', 'T', 'T', 'G', 'C', 'T', 'T']);

    assert_ulps_eq!(
      prof_values,
      array![0.21259213, 0.29005381, 0.32079767, 0.32170632, 0.25466645, 0.19585567, 0.27300161]
    );

    Ok(())
  }

  #[rstest]
  fn calculates_seq2prof() -> Result<(), Report> {
    let seq = array!['G', 'T', 'G', '-', 'G', 'G', 'C'];
    let prof = seq2prof(&seq, &ALPHABET)?;
    assert_eq!(
      prof,
      array![
        [0., 0., 1., 0., 0.],
        [0., 0., 0., 1., 0.],
        [0., 0., 1., 0., 0.],
        [0., 0., 0., 0., 1.],
        [0., 0., 1., 0., 0.],
        [0., 0., 1., 0., 0.],
        [0., 1., 0., 0., 0.]
      ]
    );
    Ok(())
  }
}
