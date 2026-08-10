# generate_subtree draws the mutating branch from weights inconsistent with its own rate

## v0 location

`TreeTime.generate_subtree()` (`#TreeTime`, `#generate_subtree`) [packages/legacy/treetime/treetime/treetime.py#L872-L1010](../../packages/legacy/treetime/treetime/treetime.py#L872-L1010)

## Erratum

The stochastic polytomy resolution sweep is a Gillespie loop: at each step it computes a
total event rate, draws a waiting time, chooses between a mutation event and a coalescence
event in proportion to their rates, and then chooses which branch or pair the event applies
to. The mutation rate and the branch-selection weights disagree.

The rates at [packages/legacy/treetime/treetime/treetime.py#L919-L920](../../packages/legacy/treetime/treetime/treetime.py#L919-L920):

```python
total_mut_rate = mutation_rate * total_mutations + coalescent_rate * n_branches_w_mutations
total_coalescent_rate = max(0, (len(ready_to_coalesce) - 1)) * (coalescent_rate + mutation_rate)
```

with `coalescent_rate` itself already carrying an added `mutation_rate` term
([L910](../../packages/legacy/treetime/treetime/treetime.py#L910),
[L912](../../packages/legacy/treetime/treetime/treetime.py#L912)). Writing $\mu$ for
`mutation_rate`, $\kappa$ for the coalescent per-branch merger rate, $M$ for
`total_mutations`, $R$ for `ready_to_coalesce` and $W$ for the alive branches that still
carry mutations:

$$R_{\text{mut}} = \mu M + (\kappa + \mu)\lvert W\rvert, \qquad
R_{\text{coal}} = \max(0, \lvert R\rvert - 1)(\kappa + 2\mu)$$

Having selected a mutation event with probability
$R_{\text{mut}}/(R_{\text{mut}} + R_{\text{coal}})$, the branch to mutate is drawn at
[packages/legacy/treetime/treetime/treetime.py#L944-L951](../../packages/legacy/treetime/treetime/treetime.py#L944-L951):

```python
p /= total_mut_rate * total_rate_inv
p *= total_mutations
branch_to_mutate = None
for b in branches_alive:
    p -= mutations_per_branch.get(b.name, 0)
    if p < 0:
        branch_to_mutate = b
        break
```

`p` is rescaled to be uniform on $[0, M]$ and then walked down the per-branch mutation
counts. Those counts sum to $M$, but the rate that gated the draw was
$\mu M + (\kappa + \mu)\lvert W\rvert$. Only the $\mu M$ portion has branches behind it. The
remaining $(\kappa + \mu)\lvert W\rvert$ portion corresponds to no branch, so whenever the
draw lands in it, the loop runs off the end of `branches_alive` with `p >= 0` and
`branch_to_mutate` stays `None`.

The fallback at [packages/legacy/treetime/treetime/treetime.py#L953-L960](../../packages/legacy/treetime/treetime/treetime.py#L953-L960)
then removes a mutation from `branches_alive[-1]`, an arbitrary branch that may have no
mutations left to remove, driving its count negative. v0's own comment identifies this as a
defect rather than intended behaviour:

```python
# this should never happen (but recreate previous behavior for now)
# TODO: raise
self.logger('TreeTime.generate_subtree: did not find a mutation to remove -- error in total_mutation count calculation', 2)
```

The failure is not rare. The misallocated fraction is
$(\kappa + \mu)\lvert W\rvert / R_{\text{mut}}$, which approaches 1 whenever the alive
branches carry few remaining mutations each -- the ordinary state late in a sweep, and the
normal state throughout for a polytomy whose children carry one or two mutations.

Separately, the three added $\mu$ terms have no stated derivation. The comment at
[packages/legacy/treetime/treetime/treetime.py#L916-L918](../../packages/legacy/treetime/treetime/treetime.py#L916-L918)
motivates the *structure* -- branches with mutations may only mutate, mutation-free branches
may only coalesce -- but that structure is already expressed by partitioning the branches
into $W$ and $R$. It does not explain adding a coalescent rate to the mutation channel, nor
adding $\mu$ to the coalescent channel twice.

## Correct formulation

Partitioning the branches already encodes the conditioning. The self-consistent rates are

$$R_{\text{mut}} = \mu M, \qquad R_{\text{coal}} = \max(0, \lvert R\rvert - 1)\,\kappa(t)$$

with the mutating branch drawn $\propto m_b$ (weights summing to $M$, matching
$R_{\text{mut}}$) and the coalescing pair drawn uniformly from $R$. Every unit of rate now
has an event behind it, so the selection loop cannot fall through.

## Evidence

- v0's own comment and `TODO: raise` at [packages/legacy/treetime/treetime/treetime.py#L953-L957](../../packages/legacy/treetime/treetime/treetime.py#L953-L957) label the reachable branch an error in the rate calculation
- The branch-selection weights at [L945](../../packages/legacy/treetime/treetime/treetime.py#L945) sum to `total_mutations`, not to `total_mut_rate`, while the acceptance test at [L940](../../packages/legacy/treetime/treetime/treetime.py#L940) uses `total_mut_rate`
- The coalescing pair is drawn uniformly from `ready_to_coalesce` at [L966](../../packages/legacy/treetime/treetime/treetime.py#L966), consistent with a rate $\propto \lvert R\rvert$ but not with the $(\kappa + 2\mu)$ per-pair factor at [L920](../../packages/legacy/treetime/treetime/treetime.py#L920)
- The docstring for `resolve_polytomies` at [packages/legacy/treetime/treetime/treetime.py#L665-L668](../../packages/legacy/treetime/treetime/treetime.py#L665-L668) describes the method as generating "a stochastic binary coalescent tree with mutation", i.e. a two-channel process, with no third rate contribution

## v0 impact

- A mutation is deducted from a branch chosen by position in `branches_alive` rather than by
  mutation count, biasing which lineages become eligible to coalesce first
- Branch mutation counts can go negative, after which
  `mutations_per_branch.get(b.name, 0) == 0` at [L908](../../packages/legacy/treetime/treetime/treetime.py#L908)
  is false forever and that branch can never coalesce, forcing it to remain a direct child
  of the polytomy
- The inflated rates shorten the drawn waiting times, compressing the resolved subtree toward
  the present relative to the intended coalescent

## v1 status

Implemented with the correction in [`simulate_subtree()`](../../packages/treetime/src/timetree/optimization/polytomy/sweep.rs). v1 uses the two rates above, selects the mutation channel from those rates, and selects a branch with an exact integer draw over the same mutation counts. Event-for-event v0 comparison is not valid because the corrected process has a different event distribution.
