# The mugration rate step without a bracket keeps the messages at the start rate

When the rate step of the mugration GTR refinement finds no minimum inside its search range, it keeps the start rate together with the backward messages and node states evaluated at the start rate. The next fit counts with these messages. v0 keeps the start rate but the messages of its last cost evaluation, at 100 times the square root of the rate, so v1 diverges from v0 in runs where a bracket fails.

**Type**: Scientific divergence from v0.

**Status**: Approved on 2026-10-21.

## Rule

`fn optimize_gtr_rate()` ([packages/treetime/src/gtr/refinement.rs#L126-L234](../../packages/treetime/src/gtr/refinement.rs#L126-L234)) evaluates the negative log likelihood at $0.01\sqrt{\mu}$, $\sqrt{\mu}$ and $100\sqrt{\mu}$. Without a bracket (the cost at $\sqrt{\mu}$ is not below both ends), `RateRangeEdge::KeepStart` returns the candidate evaluated at $\sqrt{\mu}$, with the rate set to the start rate ([packages/treetime/src/gtr/refinement.rs#L242-L262](../../packages/treetime/src/gtr/refinement.rs#L242-L262)). Only mugration uses `KeepStart`; the relative-rate step of runs with several alignment files uses `RateRangeEdge::LowerCostEnd`, which already returns the messages of the rate it takes ([kb/decisions/multi-alignment-partitions.md](multi-alignment-partitions.md)).

## v0 behavior

`optimize_gtr_rate()` restores `old_mu` after a failed bracket and keeps the profiles of the last cost evaluation ([packages/legacy/treetime/treetime/treeanc.py#L1679-L1708](../../packages/legacy/treetime/treetime/treeanc.py#L1679-L1708)). scipy evaluates the three bracket points in the order $a$, $b$, $c$ before it checks the bracket (`Brent.get_bracket_info()` of scipy 1.18.1, the version of the pixi environment), so these profiles belong to $100\sqrt{\mu}$, that is to $10^4\mu$. The next `infer_gtr()` then combines backward messages at $10^4\mu$ with transition matrices at $\mu$.

## Why

Expected transition counts combine the messages of a branch with its transition matrix under one model. Messages at $10^4\mu$ with transition matrices at $\mu$ count transitions of no model, and the fit that follows is not an estimate under either rate.

## Effect

Measured on the mugration golden-master datasets (20 tips each, 6 rate steps per run): none of the rate steps finds a bracket on `zika/20` (country, also with `?` values and with leaves missing from the metadata) and on `dengue/20` (country). Compared with v0, the confidences move by up to 0.038 on the zika cases and by up to 0.14 on `dengue/20`, where one node changes its state. On `lassa/L/20`, `rsv/a/20`, `mpox/clade-ii/20`, `tb/20`, and zika with weights, every rate step finds a bracket.

The golden masters stay a comparison with v0: the capture runs v0's own backward pass at the rate that each v0 rate step keeps ([packages/treetime/src/mugration/__tests__/__fixtures__/gm_mugration_capture](../../packages/treetime/src/mugration/__tests__/__fixtures__/gm_mugration_capture)). With a bracket, that pass moves the v0 confidences by less than 1e-8, inside the tolerance of 1e-7 ([kb/algo/mugration.md](../algo/mugration.md#golden-master-tolerance)).
