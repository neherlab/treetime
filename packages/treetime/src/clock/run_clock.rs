use crate::clock::timetree_args::TreetimeClockArgs;
use eyre::Report;
use crate::io::dates::read_dates;

pub fn run_clock(clock_args: &TreetimeClockArgs) -> Result<(), Report> {
  let TreetimeClockArgs {
    input_fastas,
    aln,
    tree,
    vcf_reference,
    dates,
    name_column,
    date_column,
    sequence_length,
    clock_filter,
    reroot,
    keep_root,
    tip_slack,
    covariation,
    allow_negative_rate,
    plot_rtt,
    outdir,
    seed,
  } = clock_args;

  let dates = read_dates(dates, name_column, date_column)?;

  dbg!(dates);

  Ok(())
}
