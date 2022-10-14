use crate::commands::mugration::mugration_args::TreetimeMugrationArgs;
use crate::io::discrete_states_csv::read_discrete_states;
use eyre::Report;

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
    sampling_bias_correction,
    outdir,
    seed,
  } = mugration_args;

  let states = read_discrete_states(states, name_column, &Some(attribute.clone()))?;

  dbg!(&states);

  Ok(())
}
