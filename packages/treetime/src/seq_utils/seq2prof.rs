use crate::gtr::gtr::GTR;
use crate::seq_utils::normalize_profile::normalize_profile;
use crate::utils::ndarray::cumsum_axis;
use eyre::Report;
use itertools_num::ItertoolsNum;
use ndarray::{s, Array1, Array2, Axis};
use ndarray_rand::rand_distr::Uniform;
use ndarray_rand::RandomExt;
use ndarray_stats::QuantileExt;
use serde_json::Value::Array;

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
pub fn prof2seq(profile: &Array2<f32>, gtr: &GTR, sample_from_prof: bool, normalize: bool) -> Result<(), Report> {
  // Normalize profile such that probabilities at each site sum to one
  let profile = &if normalize {
    normalize_profile(profile, false)?
  } else {
    profile.to_owned()
  };

  // Sample sequence according to the probabilities in the profile
  // (sampling from cumulative distribution over the different states)
  let idx = if sample_from_prof {
    let cumdis: Array1<f32> = cumsum_axis(profile, Axis(1)).t();
    let randnum = Array1::<f32>::random(cumdis.shape()[1], Uniform::new(0.0, 1.0));

    // let indices = cumdis >= randnum;
    // np.argmax(cumdis >= randnum, axis = 0);
    profile.argmax()

    // profile.fold_axis(Axis(1), Array1::<usize>::zeros(profile.shape()[0]), |res, axis| {
    //   res.argmax()
    //
    //   axis.argmax()
    // })
  } else {
    // profile.argmax(Axis(1))
    profile.argmax()
  }?;

  let seq = gtr.alphabet[idx]; // max LH over the alphabet
  let prof_values = profile[s!(.., idx)];
  return (seq, prof_values, idx);
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
pub fn seq2prof(seq: String, profile_map: Array2<f32>) {
  // np.array([profile_map[k] for k in seq])
}
