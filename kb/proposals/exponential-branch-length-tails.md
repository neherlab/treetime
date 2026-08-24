# Proposal: exponential tails on branch length distributions

## Problem

Branch length distributions currently have hard grid boundaries (`Error` tails). The right side (long branches) physically decays exponentially, but the grid truncates it. After negation and convolution, this truncation limits the backward parent message's grid extent, making disjoint-support collisions more likely when child messages disagree temporally.

The `Constant` tail on backward messages provides a crude workaround -- it declares the message evaluable beyond its grid at a flat boundary value. An exponential tail would be more physically accurate: the probability of long branches decays smoothly rather than holding flat.

## Proposal

Add a new `BoundaryBehavior` variant that extrapolates the tail using the slope at the grid boundary:

$$f(t) = f(t_0) \exp\!\left(\frac{\ln f(t_1) - \ln f(t_0)}{t_1 - t_0}\,(t - t_0)\right)$$

where $t_0$ is the boundary point, $t_1$ is the adjacent interior point, and $f(t_i)$ are the corresponding grid values. For $t < t_0$ (left extrapolation) or $t > t_0$ (right extrapolation), this continues the local log-linear slope as an exponential decay.

This matches v0's `convolve_fft` behavior, which constructs slope-based tails (linear in neg-log, exponential in probability) from the outermost trusted points at [packages/legacy/treetime/treetime/node_interpolator.py#L231-L256](../../packages/legacy/treetime/treetime/node_interpolator.py#L231-L256).

## Where it applies

- **Branch length distributions (right side)**: the primary target. Long branches are unlikely but not impossible; the tail decays exponentially. The left side is typically at zero (hard boundary -- branches can't be negative). If the left side is not at zero, exponential extrapolation could apply there too, but extending the grid to zero is likely simpler.
- **Convolution results (backward messages, left side)**: after negating the branch length and convolving, the parent message inherits the branch length's right tail as its left tail. An exponential tail here would replace the current flat `Constant` approximation.
- **Forward messages (right side)**: symmetric case. The convolution with a branch length distribution could produce an exponential right tail.

## Grid size consideration

Exponential tails require discretizing the tail onto additional grid points before convolution. The grid grows by the number of tail points. For well-behaved exponential decays, the tail reaches negligible values within a few multiples of the boundary slope, so the extension is modest. The 1M point safety cap prevents unbounded growth.

Previous grid explosion issues (fixed in `94a519a3`) were caused by a units mixup between years and divergence, not by tail extension.

## Relationship to `Constant` tail multiplication fix

The `Constant` tail fix in multiplication is complementary. It ensures multiplication reads operand tails when computing support intersection, regardless of tail shape. If exponential tails are implemented:

- `BoundaryBehavior::Exponential` would need the same treatment as `Constant` in `multiplication_support_intersection` (now shared by both multiplication and division): extend the evaluable domain on that side
- The product in the extension region would use the exponential decay values instead of a flat constant, producing a more physically accurate result

## Open questions

- Should the exponential tail be a new `BoundaryBehavior` variant (e.g. `Exponential`) or a separate mechanism on `DistributionFunction` that stores the slope parameters?
- Should `GridFn::interp` handle exponential extrapolation, or should it be applied at a higher level (convolution-specific)?
- What truncation threshold should determine how far the tail extends? A fixed probability threshold (e.g. $10^{-10}$ of the boundary value) or a fixed number of grid spacings?

## Motivation

Branch length distributions have exponential tails physically. v0 already constructs slope-based (linear in neg-log, exponential in probability) tails from the outermost trusted points in `fn NodeInterpolator.convolve_fft()`. v1 currently uses flat `Constant` tails as a simpler approximation. An exponential variant would close this parity gap and improve the accuracy of backward messages in the tail region.
