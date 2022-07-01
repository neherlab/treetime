// use crate::gtr::gtr::GTR;
// use crate::seq_utils::normalize_profile::normalize_profile;
// use crate::utils::ndarray::{argmax_axis, cumsum_axis};
// use eyre::Report;
// use ndarray::{s, Array1, Array2, Axis, Zip};
// use ndarray_rand::rand_distr::Uniform;
// use ndarray_rand::RandomExt;
//
// /// Convert profile to sequence and normalize profile across sites.
// ///
// /// Parameters
// /// ----------
// ///
// ///  profile : numpy 2D array
// ///     Profile. Shape of the profile should be (L x a), where L - sequence
// ///     length, a - alphabet size.
// ///  gtr : gtr.GTR
// ///     Instance of the GTR class to supply the sequence alphabet
// ///  collapse_prof : bool
// ///     Whether to convert the profile to the delta-function
// ///
// /// Returns
// /// -------
// ///  seq : numpy.array
// ///     Sequence as numpy array of length L
// ///  prof_values :  numpy.array
// ///     Values of the profile for the chosen sequence characters (length L)
// ///  idx : numpy.array
// ///     Indices chosen from profile as array of length L
// pub fn prof2seq(profile: &Array2<f32>, gtr: &GTR, sample_from_prof: bool, normalize: bool) -> Result<(), Report> {
//   // Normalize profile such that probabilities at each site sum to one
//   let profile = &if normalize {
//     normalize_profile(profile, false)?
//   } else {
//     profile.to_owned()
//   };
//
//   // Sample sequence according to the probabilities in the profile
//   // (sampling from cumulative distribution over the different states)
//   let idx = if sample_from_prof {
//     //cumdis = tmp_profile.cumsum(axis=1).T
//     let cumdis: Array1<f32> = cumsum_axis(profile, Axis(1))?.t().to_owned();
//
//     // randnum = np.random.random(size=cumdis.shape[1])
//     let randnum: Array1<f32> = Array1::<f32>::random(cumdis.shape()[1], Uniform::new(0.0, 1.0));
//
//     // np.argmax(cumdis >= randnum, axis = 0);
//     let mask: Array1<u8> = Zip::from(&cumdis)
//       .and(&randnum)
//       .map_collect(|&c, &r| if c > r { 1 } else { 0 });
//
//     argmax_axis(&mask, Axis(0))
//   } else {
//     // tmp_profile.argmax(axis=1)
//     argmax_axis(profile, Axis(1))
//   };
//
//   let seq = &gtr.alphabet[idx]; // max LH over the alphabet
//   let prof_values = &profile[s!(.., idx)];
//
//   Ok((seq, prof_values, idx))
// }
//
// /// Convert the given character sequence into the profile according to the
// /// alphabet specified.
// ///
// /// Parameters
// /// ----------
// ///
// ///  seq : numpy.array
// ///     Sequence to be converted to the profile
// ///
// ///  profile_map : dic
// ///     Mapping valid characters to profiles
// ///
// /// Returns
// /// -------
// ///
// ///  idx : numpy.array
// ///     Profile for the character. Zero array if the character not found
// pub fn seq2prof(seq: String, profile_map: Array2<f32>) {
//   // np.array([profile_map[k] for k in seq])
// }
