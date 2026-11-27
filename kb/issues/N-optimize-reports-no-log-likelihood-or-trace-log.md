# Optimize reports no log-likelihood and has no trace log

`treetime optimize` never writes the log-likelihood it maximizes. The model document (`--output-model`) holds W, pi and mu per partition but no likelihood, there is no `--output-tracelog`, and the only likelihood in the console output is the convergence message of the loop ("the log-likelihood is within 0.1 of the iteration before the previous one", or "the log-likelihood decreased"), which names no value. `--verbosity debug` adds no likelihood line.

Two comparisons that the partition work needs are therefore impossible from the outputs:

- **Scheme and sharing choices**: `--partition-by codon-positions`, `genes` and `rate-categories` differ in the number of free parameters, and `--shared-params` changes it by 8 per partition for a nucleotide alphabet (5 exchange rates and 3 frequencies). Without a likelihood no AIC, AICc or BIC can be computed, so [kb/reports/auto-partitioning.md](../reports/auto-partitioning.md) S6 ("no default criterion has been approved") cannot be settled on data. The only comparison left is the total tree length, which is not a fit statistic
- **Convergence diagnosis**: [kb/issues/M-optimize-loop-stops-on-decreased-likelihood-with-long-branches.md](M-optimize-loop-stops-on-decreased-likelihood-with-long-branches.md) reports that the loop often stops on a decreased likelihood, and that rate buckets and gamma make it worse. A user cannot see the size of the decrease, the iteration it happened in, or how far the returned best iteration sits from the last one

`timetree` already has the output and the infrastructure. `--output-tracelog` writes one row per iteration with the columns `n_diff`, `n_resolved`, `max_time_change`, `rms_time_change`, `log_lh_seq`, `log_lh_time`, `log_lh_coal` (`struct TraceCsvSink`, [packages/app-commands/src/commands/timetree/trace.rs#L14-L17](../../packages/app-commands/src/commands/timetree/trace.rs#L14-L17), selected by `OutputSelection::Tracelog` in [packages/app-commands/src/commands/timetree/args.rs#L237](../../packages/app-commands/src/commands/timetree/args.rs#L237)). `log_lh_seq` is the quantity `optimize` maximizes. The optimize loop computes it every iteration for its stop rule ([packages/treetime/src/optimize/run_loop.rs](../../packages/treetime/src/optimize/run_loop.rs)) and discards it.

With `--gamma-categories` the value to report is the mixture log-likelihood, which already replaces the root log-likelihood everywhere a reconstruction reports one ([kb/decisions/discrete-gamma-rate-variation.md](../decisions/discrete-gamma-rate-variation.md)).

> [!IMPORTANT]
> **Decision required.** Where the value belongs, and what counts as one number for a multi-partition run.
>
> - **Placement**: a `--output-tracelog` for `optimize` with the columns of the timetree trace log that apply, a likelihood field in the model document, a console line per iteration at `--verbosity debug`, or more than one of these
> - **Per-partition or total**: the loop maximizes one total over all partitions, but a per-partition breakdown is what a scheme comparison reads. The indel term is part of the objective unless `--no-indels` is set, so a reported value needs its composition named
> - **Parameter count**: an information criterion needs the number of free parameters, which depends on the models of the partitions, on `--shared-params`, on `--partition-rates`, on `--gamma-categories`, and on whether the topology moved. Reporting a likelihood without that count invites wrong criteria. Decide whether TreeTime reports the count, reports a criterion itself, or reports neither

## Validation

- An `optimize` run reports a log-likelihood that equals the value its stop rule compares, checked against a reconstruction run on the returned tree and models
- The reported value of a run with one partition equals the value of the same data without partition settings
- A run that stops on a decreased likelihood shows the decrease in the trace log, and the returned iteration is the one with the highest value

## Related

- [kb/issues/M-optimize-loop-stops-on-decreased-likelihood-with-long-branches.md](M-optimize-loop-stops-on-decreased-likelihood-with-long-branches.md) - the convergence behavior that a trace log would show
- [kb/reports/auto-partitioning.md](../reports/auto-partitioning.md) - model selection for partitioned analyses, blocked on this
- [kb/reports/partition-configuration.md](../reports/partition-configuration.md) - the measured scheme comparison, which reports tree lengths only for this reason
