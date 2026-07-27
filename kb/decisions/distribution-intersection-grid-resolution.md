# Distribution operations retain a bounded uniform intersection grid

v1 represents discretized distributions on uniform grids. Pointwise products and quotients therefore resample their positive-width support intersection onto a new uniform grid instead of preserving v0's union of operand knots.

**Type**: Intentional numerical divergence from v0.

## Grid contract

For an exact analytical support intersection `[overlap_min, overlap_max]`, define the target spacing from the concrete Function operands:

```text
dx = min(a.dx, b.dx)  # Function-Function
dx = function.dx      # Range-Function or Formula-Function
n_points = clamp(round((overlap_max - overlap_min) / dx) + 1, 2, 1_000_000)
```

The result grid includes both analytical endpoints. A disjoint intersection produces Empty; endpoint-only contact produces a Point. Function-Function and Range-Function apply this contract to multiplication and division. Formula-Function multiplication uses the Function spacing. Formula-Formula operations have no concrete input spacing and retain their separate fixed discretization.

Endpoint contact is an exact comparison of the analytical bounds. A tolerance would enlarge the intersection and could turn two narrowly disjoint supports into a Point. This preserves v0's endpoint-contact result without reproducing its grid construction: v0 collects knots inside an epsilon-expanded intersection and converts a single surviving knot into a delta distribution. [`packages/legacy/treetime/treetime/distribution.py#L82-L145`](../../packages/legacy/treetime/treetime/distribution.py#L82-L145)

### Division with explicit divisor tails

The nominal divisor grid bounds constrain division only on sides whose boundary behavior is `Error`. An explicit `Zero` or `Constant` tail declares the divisor evaluable beyond that grid boundary, so division extends that side to the dividend boundary. The two sides are independent:

```text
divisor_eval_min = divisor.x_min if divisor.left_extrap  == Error else dividend_min
divisor_eval_max = divisor.x_max if divisor.right_extrap == Error else dividend_max
result_bounds    = intersection(dividend_bounds, [divisor_eval_min, divisor_eval_max])
```

This preserves exact finite-support division for the default `Error` behavior while allowing inference messages with explicit tails to affect cavity division. The result grid still follows the spacing contract above. [kb/decisions/timetree-inference-pass-boundary-tails.md](timetree-inference-pass-boundary-tails.md) defines the meaning and representation constraints of each boundary behavior.

The one-million-point limit is the existing grid safety ceiling. It prevents unbounded allocations when a narrow input spacing is combined with a wide intersection.

Coalescent contributions do not introduce another grid. The backward pass multiplies child messages first and then evaluates the contribution pointwise on the completed message grid in negative-log space.

## Rationale

Uniform grids are a v1 representation invariant and support direct convolution without an additional nonuniform-to-uniform conversion. Using the finer concrete input spacing retains more operand resolution than selecting the larger operand length independently of overlap width. Exact analytical endpoints avoid support drift when a boundary lies between input samples.

The contract follows the distribution and coalescent design recorded in commit [`542ac860c7cfa4bab6764aee1d1b3810a09eb54f`](https://github.com/neherlab/treetime/commit/542ac860c7cfa4bab6764aee1d1b3810a09eb54f).

## Accepted limitations

This is a bounded resampling rule rather than a derived interpolation-error criterion. Rounding can make the realized spacing slightly finer or coarser than the target spacing, and the safety cap can make it substantially coarser. Sequential multiplication can also resample at each binary operation and is not guaranteed to be associative, although Function-Function multiplication is commutative because its support and spacing choices are symmetric.

Validation therefore covers the stated contract directly: exact support endpoints, point-count rounding and capping, use of the finer Function spacing, multiplication commutativity, pointwise values for known linear operands, default finite-support division, and division over explicitly declared tails. Scientific workflows that require posterior peak, normalization, or integrated-probability error bounds need a separate accuracy contract before tightening this approximation.
