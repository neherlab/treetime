# Timetree --vary-rate is unimplemented

The `--vary-rate` CLI flag does not exist in v1. The clap parser rejects it ("unexpected argument"). In v0, `vary_rate` is a parameter to `TreeTime.run()` derived from `--clock-std-dev` and `--covariation`/`--confidence` flags rather than a standalone flag.

## v0 behavior

v0 `wrappers.py:454-469` derives `vary_rate` as follows:

- `--clock-std-dev` provided: `vary_rate = clock_std_dev` (if `--confidence`) or `False`
- `--confidence` and `--covariation`: `vary_rate = True`
- `--confidence` without `--covariation`: warns, sets `vary_rate = False`
- Otherwise: `vary_rate = False`

This is passed to `TreeTime.run()` at `wrappers.py:493`. When active, it calls `calc_rate_susceptibility` at `treetime.py:380-382`, which redoes time tree estimation for rates +/- one standard deviation.

## v1 status

v1 has `--clock-std-dev`, `--covariation`, and `--confidence` flags. The `vary_rate` derivation logic and `calc_rate_susceptibility` behavior are not wired.

The smoke test includes `--vary-rate` as a direct CLI flag (`dev/run-smoke-tests` line 312), which is not how v0 exposes this feature. The smoke test entry should either test the underlying flags (`--clock-std-dev`, `--confidence --covariation`) or be removed.

## Related issues

- [N-timetree-dead-cli-flags.md](N-timetree-dead-cli-flags.md) -- inventory of dead/unimplemented timetree flags
