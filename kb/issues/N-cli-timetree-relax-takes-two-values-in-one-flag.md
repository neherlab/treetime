# `timetree --relax` takes two values in one flag

`timetree --relax` reads two space-separated values, the slack and the coupling of the relaxed clock (`--relax 1.0 0.5`), through `num_args = 2` with the value names `SLACK` and `COUPLING` ([packages/app-commands/src/commands/timetree/args.rs#L712-L719](../../packages/app-commands/src/commands/timetree/args.rs#L712-L719)). The config key `relax` takes a list of two numbers.

Most command-line tools take one value per flag, or several values separated by commas. A flag with two space-separated values is easy to misread in a command line, and a workflow rule (Snakemake, Nextflow) that builds arguments from a map of flag to one value needs a special case for it.

## Approaches

- Two flags, one per parameter, for example `--relax-slack` and `--relax-coupling`, with a config key for each
- One flag with comma-separated values, `--relax 1.0,0.5`, with the same list in the config
- Keep the two space-separated values, which TreeTime 0.x also uses (`--relax 1.0 0.5`)

> [!IMPORTANT]
> **Decision required.** The flag shape of the relaxed clock is open. The options are the three approaches above. The first two change the CLI spelling of TreeTime 0.x, so the CLI must then accept the 0.x form or stop with a message that names the replacement.
