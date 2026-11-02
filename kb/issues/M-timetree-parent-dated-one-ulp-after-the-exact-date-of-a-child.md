# An inferred parent can be dated one ulp after the exact date of its child

The time inference can date an inferred node one ulp (unit in the last place) after the exact date of one of its children, although the posterior of the node ends at that date. The output tree then holds a branch of negative length, $-2.27 \cdot 10^{-13}$ years at the year 2026, and the node data dates the parent after its child.

The consumers of node times inside TreeTime compare them with `fn runs_backwards()` [packages/treetime/src/timetree/inference/result.rs](../../packages/treetime/src/timetree/inference/result.rs), which treats times that `fn coordinates_coincide()` finds equal (4 ulps of the larger magnitude) as one point: polytomy resolution accepts such a pair, and the warning about branches that run backwards in time does not count it. The dates of the outputs keep the one-ulp inversion.

## Evidence

- **Production build**: in `lassa-gpc`, a Nextstrain build replayed through `augur refine` (`--keep-root --clock-rate 0.00078 --coalescent opt --max-iter 2`), one backward pass sends 3480 messages from exactly dated leaves to their parents. 2780 of them end exactly at the leaf date, 668 end one ulp after it, and 32 end one ulp before it. Each pass counted 92 branches that run backwards, each one ulp long, and `timetree.nwk` holds 92 branch lengths of `-2.27e-13`
- **Earlier counts**: every branch that the backward-branch warning counted was exactly one ulp backwards in 368 counts over the passes of lassa GPC, 38 of norovirus p22, 7 of measles N450, and 19 of the H5N1 cattle-outbreak HA build
- **Rounding dependence**: `data/flu/h3n2/500` with `--seed=2251488387996782098` met the one-ulp pair in polytomy resolution with the x86_64 glibc dist binaries of `rust` `926cdf77b` and not with the musl binary, because the operands of the shift below come from `exp` and `log` of the C library
- **Property test**: [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md) records a simulated case whose tip, dated exactly at 2024.3333333333333, lies $2 \cdot 10^{-13}$ years before its parent

## Mechanism

- **Message from an exactly dated child**: `fn send_backward_message()` in `packages/treetime/src/timetree/inference/backward_pass.rs` negates the branch-length distribution, whose grid then ends at exactly 0, and convolves the point distribution of the child date $d$ with it. `fn convolution_point_function()` in `packages/treetime-distribution/src/distribution_ops/convolve.rs` builds the shifted grid from its start, `f.x_min() + p.t()`, and the spacing
- **Grid end**: `struct Grid` in `packages/treetime-grid/src/grid.rs` stores the start, the spacing, and the point count, and computes its end as `x_min + dx * (n - 1)`. The end of the message is therefore $\mathrm{fl}(\mathrm{fl}(d - A) + A)$ with $A = dx \cdot (n - 1)$, which differs from $d$ by one ulp when the two roundings do not cancel. The hard-end law of the message is shifted exactly to $d$ (`fn shift_tail_law()`), so the grid end and its law disagree
- **Node time**: the product at the parent ends at the smallest end of the messages of its children, and the most likely time lies at that end when the branch carries no substitution. `fn committed_time()` in `packages/treetime/src/timetree/inference/forward_pass.rs` commits it, and the child keeps its exact date ([kb/decisions/timetree-uncertain-leaf-dates-are-inferred.md](../decisions/timetree-uncertain-leaf-dates-are-inferred.md))
- **Same pattern**: `Grid::from_range_dx_anchored()` with `Side::Right` computes the start from the anchored end and then recomputes the end from the start, so a grid anchored at its right end can also end one ulp away from the anchor

> [!IMPORTANT]
> **Decision required.** The fix changes how a grid represents its ends, which every distribution operation uses. Options:
>
> - **Exact grid ends** (recommended): `struct Grid` keeps both ends as given coordinates and computes only the interior points; the shift, the negation, and the anchored constructors map the ends exactly. This covers every operation that shifts, negates, or anchors a grid. It moves the end points of grids by up to one ulp wherever they differ today; how far that changes the outputs of `timetree` is not measured
> - **Local start correction**: `fn convolution_point_function()` moves the start of the shifted grid by single ulps until the recomputed end equals the shifted end. This covers only this operation and leaves the anchored constructors as they are
> - **Snap at the commit**: the forward pass commits an inferred node time that coincides with the exact date of a child at that date. This repairs the symptom in the outputs only, and belongs to the projection of node times that [M-timetree-marginal-node-times-can-violate-topology.md](M-timetree-marginal-node-times-can-violate-topology.md) leaves to a decision

## Validation

- On `lassa-gpc`, every backward message from an exactly dated leaf ends exactly at the leaf date, and `timetree.nwk` holds no negative branch length
- A unit test convolves a point at a date near 2026 with a negated branch distribution whose shifted end rounds away from the date, and finds the end of the result at the date
