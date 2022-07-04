use crate::make_report;
use crate::utils::ndarray::max_axis;
use eyre::Report;
use ndarray::prelude::*;
use ndarray_einsum_beta::einsum;
use num_traits::real::Real;
use num_traits::FloatConst;

/// Returns a normalized version of a profile matrix
///
/// Parameters
/// ----------
/// in_profile : np.array
///     shape Lxq, will be normalized to one across each row
/// log : bool, optional
///     treat the input as log probabilities
/// return_offset : bool, optional
///     return the log of the scale factor for each row
///
/// Returns
/// -------
/// tuple
///     normalized profile (fresh np object) and offset (if return_offset==True)
pub fn normalize_profile(in_prof: &Array2<f32>, log: bool) -> Result<(Array2<f32>, Array1<f32>), Report> {
  let (tmp_prof, tmp_prefactor) = if log {
    let tmp_prefactor = max_axis(in_prof, Axis(1))?;

    let tmp_prof = {
      let tmp_prof = &in_prof.t() - &tmp_prefactor;
      let tmp_prof = tmp_prof.mapv(f32::exp);
      let tmp_prof = tmp_prof.t();
      tmp_prof.to_owned()
    };

    (tmp_prof, tmp_prefactor)
  } else {
    let tmp_prefactor: Array1<f32> = Array1::zeros(in_prof.shape()[0]);
    (in_prof.to_owned(), tmp_prefactor)
  };

  let norm_vector = &tmp_prof.sum_axis(Axis(1));
  let one_div_norm_vector = 1.0 / norm_vector;

  let norm_prof = einsum("ai,a->ai", &[&tmp_prof, &one_div_norm_vector])
    .map_err(|err| make_report!("einsum: {err}"))?
    .into_dimensionality::<Ix2>()?;

  let offset: Array1<f32> = norm_vector.mapv(|x| f32::log(x, f32::E())) + tmp_prefactor;

  Ok((norm_prof, offset))
}

#[allow(clippy::excessive_precision)]
#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use rstest::rstest;

  #[rstest]
  fn normalizes_profile_non_log() -> Result<(), Report> {
    let prof = array![
      [0.54881350, 0.71518937, 0.60276338, 0.54488318, 0.42365480],
      [0.64589411, 0.43758721, 0.89177300, 0.96366276, 0.38344152],
      [0.79172504, 0.52889492, 0.56804456, 0.92559664, 0.07103606],
      [0.08712930, 0.02021840, 0.83261985, 0.77815675, 0.87001215],
      [0.97861834, 0.79915856, 0.46147936, 0.78052918, 0.11827443],
      [0.63992102, 0.14335329, 0.94466892, 0.52184832, 0.41466194],
      [0.26455561, 0.77423369, 0.45615033, 0.56843395, 0.01878980],
    ];

    let (norm_prof, offset) = normalize_profile(&prof, false)?;

    assert_ulps_eq!(
      norm_prof,
      array![
        [0.19356424, 0.25224431, 0.21259213, 0.19217803, 0.14942128],
        [0.19440831, 0.13170981, 0.26841564, 0.29005381, 0.11541244],
        [0.27439982, 0.18330691, 0.19687558, 0.32079767, 0.02462001],
        [0.03366488, 0.00781195, 0.32170632, 0.30066296, 0.33615390],
        [0.31185458, 0.25466645, 0.14705881, 0.24872985, 0.03769030],
        [0.24016971, 0.05380214, 0.35454510, 0.19585567, 0.15562739],
        [0.12705805, 0.37184099, 0.21907519, 0.27300161, 0.00902417],
      ]
    );

    assert_ulps_eq!(
      offset,
      array![1.04214924, 1.20067495, 1.05962792, 0.9509381, 1.14360473, 0.97999897, 0.73340744]
    );
    Ok(())
  }

  #[rstest]
  fn normalizes_profile_log() -> Result<(), Report> {
    let prof = array![
      [0.54881350, 0.71518937, 0.60276338, 0.54488318, 0.42365480],
      [0.64589411, 0.43758721, 0.89177300, 0.96366276, 0.38344152],
      [0.79172504, 0.52889492, 0.56804456, 0.92559664, 0.07103606],
      [0.08712930, 0.02021840, 0.83261985, 0.77815675, 0.87001215],
      [0.97861834, 0.79915856, 0.46147936, 0.78052918, 0.11827443],
      [0.63992102, 0.14335329, 0.94466892, 0.52184832, 0.41466194],
      [0.26455561, 0.77423369, 0.45615033, 0.56843395, 0.01878980],
    ];

    let (norm_prof, offset) = normalize_profile(&prof, true)?;

    assert_ulps_eq!(
      norm_prof,
      array![
        [0.19550790, 0.23089814, 0.20634523, 0.19474100, 0.17250773],
        [0.19106125, 0.15513367, 0.24431856, 0.26252930, 0.14695721],
        [0.23820265, 0.18314747, 0.19045983, 0.27232423, 0.11586582],
        [0.12156150, 0.11369386, 0.25618782, 0.24260819, 0.26594863],
        [0.27208969, 0.22739122, 0.16222638, 0.22319427, 0.11509845],
        [0.21496871, 0.13083340, 0.29155841, 0.19102795, 0.17161153],
        [0.16630235, 0.27685269, 0.20142200, 0.22535700, 0.13006596],
      ]
    );

    assert_ulps_eq!(
      offset,
      array![2.180968, 2.30105533, 2.22635853, 2.19446427, 2.28024188, 2.17718383, 2.0585034]
    );

    Ok(())
  }
}
