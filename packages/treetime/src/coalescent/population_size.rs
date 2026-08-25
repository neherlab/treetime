/// Effective population size `N_e` implied by a coalescent time scale.
///
/// The Kingman coalescent parameterises lineage merging by a time scale $T_c$: loosely, the
/// expected pairwise coalescence time in the tree's time unit. TreeTime works in calendar years,
/// so $T_c$ is already in years. Rescaling by the number of generations per year re-expresses the
/// same quantity as an effective population size in generation units, the standard vertical axis
/// of a skyline plot:
///
/// $$N_e = T_c \cdot g$$
///
/// where $g$ is `gen_per_year`. The relation is linear, so a confidence band on $T_c$ maps onto a
/// band on $N_e$ by scaling each bound with the same factor; callers that hold a band apply this
/// function to the lower and upper bounds as well as to the point estimate.
///
/// $N_e$ is a reporting quantity: it does not enter the inference.
pub fn effective_population_size(tc: f64, gen_per_year: f64) -> f64 {
  tc * gen_per_year
}
