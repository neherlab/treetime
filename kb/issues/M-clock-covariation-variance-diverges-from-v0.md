# Clock covariation variance diverges from v0

With `--covariation`, the `clock` command builds the branch variance of the root-to-tip regression in the app layer, and three parts of it differ from v0 and from v1 timetree.

- v1 `clock`: [packages/app-commands/src/commands/clock/run.rs#L58-L70](../../packages/app-commands/src/commands/clock/run.rs#L58-L70)
- v1 timetree: `fn build_covariation_clock_params` in [packages/treetime/src/timetree/params.rs#L64-L94](../../packages/treetime/src/timetree/params.rs#L64-L94)
- v0 (both commands): `ClockTree.setup_TreeRegression` in [packages/legacy/treetime/treetime/clock_tree.py#L275-L287](../../packages/legacy/treetime/treetime/clock_tree.py#L275-L287)

## v0 formula

With covariation, v0 gives each branch the variance

$$
\sigma^2_e = \left(\ell_e + [e \text{ terminal}] \cdot s^2 \cdot \frac{1}{L}\right) \cdot \frac{1}{L}
$$

where $\ell_e$ is the branch length (`clock_length` when set, else `mutation_length`), $s$ is `tip_slack` and $L$ is the sequence length (`one_mutation = 1/L`). In v1 terms this is `variance_factor = 1/L`, `variance_offset = 0`, `variance_offset_leaf = s²/L²`. v0 has no separate overdispersion factor.

## Divergences

- **Invented factor 2**: `clock` sets `overdispersion = 2.0` and uses `variance_factor = 2/L`. v0 and v1 timetree use `1/L`. The factor doubles the branch term relative to the tip term, so tips get relatively less slack than in v0
- **Sequence length ignores the alignment**: `clock` requires `--sequence-length` for `--covariation` even when `--alignment` is given, although its help text says the length is "Not required if alignment is provided" ([packages/app-commands/src/commands/clock/args.rs#L185-L187](../../packages/app-commands/src/commands/clock/args.rs#L185-L187)). v0 accepts either `--aln` or `--sequence-length` ([packages/legacy/treetime/treetime/wrappers.py#L945-L947](../../packages/legacy/treetime/treetime/wrappers.py#L945-L947)) and takes $L$ from the alignment when one is given, and so does v1 timetree
- **Variance flags ignored under covariation**: `--variance-factor`, `--variance-offset` and `--variance-offset-leaf` ([packages/app-commands/src/commands/clock/args.rs#L429-L442](../../packages/app-commands/src/commands/clock/args.rs#L429-L442)) have no v0 counterpart, and `--covariation` replaces them without a warning

## Matches v0

- `--tip-slack` is wired, and its default 3 matches the v0 CLI default ([packages/legacy/treetime/treetime/argument_parser.py#L180-L185](../../packages/legacy/treetime/treetime/argument_parser.py#L180-L185)), which v0 applies to the clock command at [packages/legacy/treetime/treetime/wrappers.py#L965](../../packages/legacy/treetime/treetime/wrappers.py#L965)
- Without covariation, the v1 defaults (factor 0, offset 0, leaf offset 1) equal v0's variance of 1 for terminal branches and 0 for internal ones

## Not verified

- With `--covariation` and without `--keep-root`, v0 runs `myTree.run(root='least-squares', max_iter=0)` before rerooting ([packages/legacy/treetime/treetime/wrappers.py#L979-L981](../../packages/legacy/treetime/treetime/wrappers.py#L979-L981)). This sets `clock_length`, which then replaces the input branch length in the variance. v1 `clock` always uses the input branch length. The size of this difference was not measured

## Fix

Compute the clock command's covariation variance with the same core function timetree uses, so both commands apply the v0 formula from one place:

- Make `fn build_covariation_clock_params` reachable from `app-commands` (it is `pub(crate)` in `treetime`) and call it from the clock command, passing the clock command's alignment (when given), `--sequence-length` and `--tip-slack`
- Delete the local `overdispersion = 2.0` construction
- Keep the clock command's `--tip-slack` default at 3 (v0 CLI default). The timetree default is tracked in [M-timetree-covariation-tip-slack-default-differs-from-v0.md](M-timetree-covariation-tip-slack-default-differs-from-v0.md); once that issue is resolved, one shared default serves both commands
- Warn when `--covariation` overrides an explicitly given `--variance-factor`, `--variance-offset` or `--variance-offset-leaf`

## Validation

- `clock --covariation --alignment=<aln>` without `--sequence-length` succeeds and uses the alignment length
- The variance parameters equal the v0 formula for a hand-computed $L$ and tip slack (oracle: [packages/legacy/treetime/treetime/clock_tree.py#L277-L285](../../packages/legacy/treetime/treetime/clock_tree.py#L277-L285))
- Smoke on the clock cases with `--covariation`: changed clock rates in those cases are expected and must be explained by the factor change

## Related

- [M-timetree-covariation-tip-slack-default-differs-from-v0.md](M-timetree-covariation-tip-slack-default-differs-from-v0.md): timetree's tip-slack default
- [M-clock-filter-regression-uses-covariation.md](M-clock-filter-regression-uses-covariation.md): the outlier filter's regression
