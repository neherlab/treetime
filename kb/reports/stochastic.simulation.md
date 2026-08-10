# Stochastic Simulation and Continuous-Time Markov Processes in TreeTime

## Scope and reader orientation

TreeTime is a phylogenetic **inference** tool, not a forward simulator. It does not evolve populations or sequences forward in time by drawing random events. Yet the same mathematical object that underlies exact stochastic simulation -- a continuous-time, pure-jump Markov process built from event rates and exponential holding times -- is the substrate for three things TreeTime does compute: the coalescent prior on node times, the substitution model on branches, and the marginal reconstruction of ancestral states.

This report maps that substrate onto TreeTime v1 as implemented: where event-rate reasoning appears in the code, how the implementation relates to the established theory, and which methods from the wider stochastic-simulation literature fall outside the current design.

Three findings anchor the report:

- **F1. TreeTime contains a single event-driven simulator**: the stochastic polytomy resolver, a two-channel <a id="gloss-use-1"></a>Gillespie <sup>[1](#gloss-1)</sup> sampler. It is the only place TreeTime draws a random trajectory. Verified by source trace.
- **F2. The coalescent and the substitution model are event processes evaluated _analytically_, never simulated.** This finding has two independent halves: **F2a** (Section 3), the Kingman coalescent, which enters as a closed-form log-prior on node times; and **F2b** (Section 4), the substitution <a id="gloss-use-2"></a>CTMC <sup>[2](#gloss-2)</sup>, which enters through the matrix exponential $P(t)=e^{Qt}$ and <a id="gloss-use-3"></a>Felsenstein pruning <sup>[3](#gloss-3)</sup>. Both choices match the standard recommendation that a finite-state endpoint likelihood is computed with matrix methods, not path simulation. Verified by source trace.
- **F3. Marginal ancestral reconstruction is an <a id="gloss-use-4"></a>endpoint-conditioned CTMC <sup>[4](#gloss-4)</sup> computation** -- the two-pass sum-product that returns conditional state distributions given the observed tips. It is the deterministic counterpart of <a id="gloss-use-5"></a>stochastic mapping <sup>[5](#gloss-5)</sup>, and it relates to the <a id="gloss-use-6"></a>Doob $h$-transform <sup>[6](#gloss-6)</sup> as the framework for conditioning a Markov process on future observations. Verified by source trace.

---

## 1. Background: the shared jump-process construction

### 1.1 Doob's construction

Every continuous-time, finite-or-countable-state Markov jump process is built the same way <a id="cite-1"></a>[Doob 1945](https://doi.org/10.2307/1990339) [[1](#ref-1)]. Let the current state be $x$, and let $q(x,y)$ be the transition rate from $x$ to a distinct state $y$. Define the total exit rate

$$q(x) = \sum_{y \neq x} q(x,y).$$

If $q(x) > 0$, the process holds in $x$ for an exponentially distributed time with rate $q(x)$, and then jumps to state $y$ with probability $q(x,y)/q(x)$. Conditional on $x$, the holding time and the destination are independent. This factorization -- total rate, exponential waiting time, rate-proportional event choice, state update -- is the universal skeleton of exact event-driven simulation.

### 1.2 The Gillespie direct method

Applied to reacting systems, this skeleton is the Gillespie stochastic simulation algorithm <a id="cite-2"></a>[Gillespie 1976](<https://doi.org/10.1016/0021-9991(76)90041-3>) [[2](#ref-2)]; <a id="cite-3"></a>[Gillespie 1977](https://doi.org/10.1021/j100540a008) [[3](#ref-3)]. Each event channel $j$ has a <a id="gloss-use-7"></a>**propensity** <sup>[7](#gloss-7)</sup> $a_j(x)$ -- its instantaneous rate given the current state. One step of the _direct method_ is:

$$a_0 = \sum_j a_j(x), \qquad \tau = -\frac{\ln r_1}{a_0}, \qquad \mu = \min\Big\{ k : \sum_{j=1}^{k} a_j \ge r_2\, a_0 \Big\},$$

where $a_0$ is the total propensity, $\tau$ is the next-event waiting time, $\mu$ is the selected channel, and $r_1, r_2$ are independent uniforms on $(0,1)$. The state advances by $\mu$'s state-change, the clock advances by $\tau$, and the propensities are recomputed. The result is an _exact_ sample of the process trajectory, with no time-discretisation error.

Two properties of the exponential make this correct and also make it composable:

- <a id="gloss-use-8"></a>**Memorylessness** <sup>[8](#gloss-8)</sup>. The residual waiting time has the same exponential law regardless of elapsed time. A candidate event that is discarded because a rate changed can be redrawn under the new rate with no bias.
- **Piecewise-constant handling.** When rates change at a _deterministic boundary_ (a scheduled time, an arrival, a rate breakpoint), the correct procedure is to integrate the hazard only up to that boundary, apply the change, and continue -- never to carry a stale rate across it, and never to fire an event drawn past the boundary.

TreeTime's event-driven sampler follows both rules; its v0 predecessor violated them (Section 2).

---

## 2. F1 -- The stochastic polytomy resolver

### 2.1 Purpose

Tree builders emit multifurcations (polytomies) wherever the sequence data lack the signal to order divergences -- for example when children carry one or two mutations each. TreeTime collapses builder-inserted zero-length bifurcations back into polytomies and re-resolves them into dated binary subtrees consistent with the temporal ordering, following the TreeTime pipeline <a id="cite-4"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[4](#ref-4)]. The resolver samples one mutation-conditioned binary history per polytomy. It is a forward-in-time event process, and the only one in TreeTime.

### 2.2 The two-channel event model

The sweep moves from the most recent child of the polytomy toward the parent. A child becomes a _live_ lineage when the sweep reaches its calendar time. Two event channels compete. Let $M$ be the total remaining substitution count over live lineages, let $R$ be the set of live mutation-free lineages, let $\mu$ be the whole-alignment mutation rate, and let $\kappa(t)$ be the whole-tree per-lineage coalescent merger rate schedule (Section 3). The propensities are

$$R_{\mathrm{mut}} = \mu M, \qquad R_{\mathrm{coal}} = \max\!\big(0,\, |R| - 1\big)\,\kappa(t).$$

A mutation event is selected against a coalescent event in proportion to these two rates -- the Gillespie channel choice. A mutation event then selects a live lineage in exact proportion to its integer remaining mutation count; a coalescent event selects a pair from $R$ uniformly. A lineage may only merge after every substitution mapped to its incoming branch has been consumed. This is the direct-method factorisation of Section 1.2, with a labelled forest as the state instead of molecule counts. The rate/multiplicity source is `packages/treetime/src/timetree/optimization/polytomy/sweep.rs` [packages/treetime/src/timetree/optimization/polytomy/sweep.rs](../../packages/treetime/src/timetree/optimization/polytomy/sweep.rs), with the merger-rate schedule assembled in `packages/treetime/src/timetree/pipeline.rs` [packages/treetime/src/timetree/pipeline.rs](../../packages/treetime/src/timetree/pipeline.rs).

The first-merger waiting time has analytical mean $1/\big((k-1)\kappa\big)$ for constant $\kappa$ and $k$ live lineages, and a unit test checks the sampler against it.

### 2.3 Time-varying hazard and deterministic boundaries

The sampler draws **one unit-exponential hazard threshold** $E \sim \mathrm{Exponential}(1)$ per event and integrates it across intervals where the live-lineage set and $\kappa(t)$ are constant, i.e. it solves for the time $T$ at which

$$\int \big(R_{\mathrm{mut}} + R_{\mathrm{coal}}\big)\,dt = E$$

with the integrand held piecewise-constant. The deterministic boundaries are lineage arrivals (serial sampling of the children by date), merger-rate breakpoints (skyline steps of $\kappa$), and the parent time. Crossing an arrival or a rate breakpoint carries the _remaining_ hazard into the next interval -- valid because the exponential is memoryless. Crossing the parent ends the sweep with no event, leaving the surviving lineages as a residual polytomy.

### 2.4 The three corrections over the reference implementation

The v0 (Python) resolver was a Gillespie loop with three defects, each a violation of the boundary/consistency rules in Section 1. TreeTime v1 corrects all three; the corrections are the reason event-for-event parity with v0 is neither achievable nor a target.

- **C1. Self-consistent channel rates.** v0 added a coalescent term to the mutation channel and a doubled mutation term to the coalescent channel, so the branch-selection weights (summing to $M$) were inconsistent with the rate that gated the draw (containing extra $(\kappa+\mu)|W|$ mass with no branch behind it). The selection loop could then run off the end and drive a mutation count negative, permanently barring a lineage from coalescing. v1 uses $R_{\mathrm{mut}}=\mu M$ with weights $\propto m_b$ that sum to exactly $M$, so every unit of rate has an event behind it. See [kb/v0-errata/timetree-stochastic-resolve-rate-selection-mismatch.md](../v0-errata/timetree-stochastic-resolve-rate-selection-mismatch.md).
- **C2. No event past the parent bound.** v0 tested the parent bound only at the top of the loop, so an unbounded exponential draw could place a merger _older than the parent_, yielding a negative branch length. v1 makes the parent time an explicit hazard boundary: crossing it ends the sweep, and plan application in `packages/treetime/src/timetree/optimization/polytomy/apply.rs` [packages/treetime/src/timetree/optimization/polytomy/apply.rs](../../packages/treetime/src/timetree/optimization/polytomy/apply.rs) rejects any merger at or before the parent before the graph is mutated. See [kb/v0-errata/timetree-stochastic-resolve-event-past-parent.md](../v0-errata/timetree-stochastic-resolve-event-past-parent.md).
- **C3. No skipped interval at an arrival.** v0 detected a lineage arrival crossing but resumed from the already-advanced time $t+\mathrm{d}t$, never evaluating the interval $[t_a,\, t+\mathrm{d}t]$ at the higher post-arrival rate -- a systematic overestimate of waiting times. v1 resumes exactly at the crossed boundary $t_a$, so each interval is evaluated once under its own rate. See [kb/v0-errata/timetree-stochastic-resolve-skipped-arrival-interval.md](../v0-errata/timetree-stochastic-resolve-skipped-arrival-interval.md).

The resolver is deterministic given its seed (`--seed`), and its correctness contract (finite non-overflowing rates, no merger more recent than either child or at/after the parent, each lineage referenced once, validation before graph mutation) is enforced before any topology change. The design decision is recorded in [kb/decisions/timetree-stochastic-polytomy-resolution.md](../decisions/timetree-stochastic-polytomy-resolution.md).

**Interpretation boundary.** This process samples _plausible_ mutation-conditioned histories to give the downstream time inference a resolved topology. It is not a calibrated posterior over binary and unresolved topologies; a soft polytomy reflects missing signal, and no sampler can manufacture signal the data do not contain.

---

## 3. F2a -- The Kingman coalescent as an analytic prior

### 3.1 The rate and the closed-form contribution

The coalescent <a id="cite-5"></a>[Kingman 1982](<https://doi.org/10.1016/0304-4149(82)90011-4>) [[5](#ref-5)] supplies a prior on node times. With $k(t)$ ancestral lineages and coalescent time scale $T_c(t)$, the classical rates are the pairwise (all-pairs) form

$$\lambda(t) = \frac{k(t)\big(k(t)-1\big)}{2\,T_c(t)} = \frac{\binom{k}{2}}{T_c(t)}, \qquad \kappa(t) = \frac{k(t)-1}{2\,T_c(t)},$$

where $\lambda$ is the total merger rate and $\kappa$ the per-lineage rate. In code, both are written through an effective count $n = \max(0.5,\, k-1)$ as $\kappa = 0.5\,n/T_c$ and $\lambda = 0.5\,n(n+1)/T_c$; in the ordinary regime $k \ge 2$ this is exactly the pairwise form above, and the clamp only alters the degenerate $k \le 1$ tail. The rate primitives are in `packages/treetime/src/coalescent/integration.rs` [packages/treetime/src/coalescent/integration.rs#L73-L88](../../packages/treetime/src/coalescent/integration.rs#L73-L88).

TreeTime does **not** simulate coalescent trajectories. It stores the cumulative per-lineage hazard $H(t) = \int_t^{P} \kappa(s)\,ds$ (with $P$ the most recent event) and evaluates a closed-form negative-log contribution per node, grouped from per-edge survival terms $H(t_p)-H(t_c)$:

| Node role             | Negative-log contribution                     |
| --------------------- | --------------------------------------------- |
| Internal (incl. root) | $\big(m-1\big)\big(H(t) - \ln\lambda(t)\big)$ |
| Leaf                  | $-H(t_{\mathrm{leaf}})$                       |
| Root correction       | $+H(t_{\mathrm{root}})$                       |

where $m$ is the node's child count. These telescope into the whole-tree Kingman objective. The construction and ownership are documented in [kb/algo/coalescent-contribution-refactor.md](../algo/coalescent-contribution-refactor.md).

### 3.2 Serial sampling without a simulation loop

Heterochronous (serially sampled) data are handled without any event loop over sampling boundaries. The lineage count $k(t)$ is built directly as a piecewise-constant step function from event deltas: each leaf contributes $+1$ at its calendar time and each $m$-child internal node contributes $-(m-1)$, aggregated in ascending calendar order in `packages/treetime/src/coalescent/events.rs` [packages/treetime/src/coalescent/events.rs#L57-L63](../../packages/treetime/src/coalescent/events.rs#L57-L63). Sampling up-steps and coalescence down-steps interleave automatically. This is the analytic equivalent of admitting sampled lineages at their collection times.

### 3.3 Time-varying $T_c$: constant, optimal, and skyline

$T_c$ is carried as a distribution and appears in three forms, all reducing to the same model primitives:

- **Constant $T_c$**, with the maximum-likelihood value available in closed form. The constant-$T_c$ Kingman likelihood is $\log L = -M\ln T_c - I/T_c + \text{const}$, where $M$ is the total merger count and $I=\int \tfrac{k(k-1)}{2}\,dt$ is the pairwise-rate integral; setting $\partial_{T_c}\log L=0$ gives $T_c^{*}=I/M$ with no iterative search. See [kb/decisions/coalescent-analytic-tc-optimization.md](../decisions/coalescent-analytic-tc-optimization.md).
- **Skyline** (piecewise-constant $T_c(t)$), a step function estimated by a convex Newton solve in $z=\ln T_c$ with a smoothness penalty on adjacent log-steps. This is the frequentist analogue of the Bayesian skyline family for inferring population-size history from a genealogy <a id="cite-6"></a>[Pybus, Rambaut, and Harvey 2000](https://doi.org/10.1093/genetics/155.3.1429) [[6](#ref-6)]; <a id="cite-7"></a>[Drummond et al. 2005](https://doi.org/10.1093/molbev/msi103) [[7](#ref-7)]. See [kb/decisions/coalescent-skyline-convex-log-tc.md](../decisions/coalescent-skyline-convex-log-tc.md).

When effective population size varies continuously the coalescence hazard is time-dependent, and freezing the rate until the next event would give the wrong waiting-time law. TreeTime avoids this by integrating $H(t)$ rather than sampling. The cumulative hazard is built by midpoint quadrature over the lineage-count breakpoints in `packages/treetime/src/coalescent/integration.rs` [packages/treetime/src/coalescent/integration.rs#L30-L66](../../packages/treetime/src/coalescent/integration.rs#L30-L66).

### 3.4 Actual merger multiplicity (a corrected reference defect)

At a polytomy node with $m$ children, the merger-rate credit carries a factor $(m-1)/m$. The reference v0 code used the correct $m$ in per-node inference but silently defaulted to $m=2$ in the whole-tree likelihood and in $T_c$ optimisation, underweighting merger terms at multifurcations by 25-33% for $m=3,4$. TreeTime uses the parent's actual child count throughout, read from the graph in `packages/treetime/src/coalescent/edge_data.rs` [packages/treetime/src/coalescent/edge_data.rs#L84-L96](../../packages/treetime/src/coalescent/edge_data.rs#L84-L96). This has zero effect on fully resolved binary trees and a real effect on trees with residual polytomies. See [kb/decisions/coalescent-total-lh-actual-multiplicity.md](../decisions/coalescent-total-lh-actual-multiplicity.md) and [kb/v0-errata/coalescent-total-lh-fixed-multiplicity.md](../v0-errata/coalescent-total-lh-fixed-multiplicity.md).

### 3.5 How the coalescent is consumed

Three uses, all analytic: (i) a log-prior applied to internal- and root-node time distributions in the backward belief-propagation pass, `packages/treetime/src/timetree/inference/backward_pass.rs` [packages/treetime/src/timetree/inference/backward_pass.rs](../../packages/treetime/src/timetree/inference/backward_pass.rs); (ii) $T_c$ estimation from the timed tree; (iii) a whole-tree likelihood term reported for convergence. The prior reads a lineage count $k(t)$ that is frozen after the first coalescent-free pass, so the loop does not chase a receding target ([kb/decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md](../decisions/timetree-frozen-lineage-counts-for-coalescent-prior.md)); the same $\kappa(t)$ schedule feeds the polytomy resolver of Section 2.

---

## 4. F2b -- Substitution as an analytic CTMC, not an event simulation

### 4.1 Matrix exponential over branches

Sequence evolution is a finite-state CTMC on the nucleotide or amino-acid alphabet. Each general-time-reversible model defines a rate matrix $Q$, and the branch transition matrix is $P(t) = e^{Qt}$ <a id="cite-8"></a>[Tavaré 1986](https://openalex.org/W1593676244) [[8](#ref-8)]. Because only the endpoint likelihood is needed -- not a realised substitution history -- TreeTime computes $P(t)$ analytically and never draws individual substitution events.

$P(t)$ is obtained by eigendecomposition, reusing the reversibility structure. Detailed balance $\pi_i Q_{ij} = \pi_j Q_{ji}$ permits symmetrisation $S = \Pi^{1/2} Q\, \Pi^{-1/2}$ into a real symmetric matrix with real eigenvalues $\lambda_i$ and orthogonal eigenvectors $V$, so that

$$P(t) = \Pi^{-1/2}\, V\, \operatorname{diag}\!\big(e^{\mu\,\lambda_i t}\big)\, V^{\top}\, \Pi^{1/2},$$

where $\Pi=\operatorname{diag}(\pi)$ holds the equilibrium frequencies $\pi_i$, $\mu$ is the overall rate scale, and $\lambda_i$ are the eigenvalues (one is always $0$, the rest negative, guaranteeing convergence to $\pi$). The eigendecomposition is done once per model in `packages/treetime/src/gtr/gtr.rs` [packages/treetime/src/gtr/gtr.rs#L74-L105](../../packages/treetime/src/gtr/gtr.rs#L74-L105); each branch then costs $O(k^2)$ via `fn Gtr.expQt()` [packages/treetime/src/gtr/gtr.rs#L453-L463](../../packages/treetime/src/gtr/gtr.rs#L453-L463). The model hierarchy (JC69 through GTR) is documented in [kb/algo/gtr.md](../algo/gtr.md).

### 4.2 Substitution is computed, not simulated

A source trace confirms that TreeTime contains no event-by-event, tau-leaping, or chemical-Langevin sequence evolution. For a finite-state substitution endpoint, matrix exponentiation with pruning yields the exact likelihood at lower cost than simulating every substitution. Gillespie-style sequence simulation only becomes attractive when the _full substitution path_, sequence-context dependence, or variable-length (indel) states are required. At pandemic scale, such simulators exist -- for example `phastSim`, which uses a Gillespie algorithm with efficient search structures for large, closely related viral and bacterial phylogenies <a id="cite-10"></a>[De Maio et al. 2022](https://doi.org/10.1371/journal.pcbi.1010056) [[10](#ref-10)] -- but they solve the _simulation_ problem, not the _endpoint-likelihood_ problem that TreeTime solves. The two occurrences of "Poisson" in TreeTime are the indel term in per-edge branch-length likelihood and the Poisson branch-length distributions used to discretise time; neither simulates substitution events.

### 4.3 Substitution counts as a Poisson branch-length model: how its variance reaches node dates

The branch-length likelihood models the per-branch substitution count as Poisson -- exact for equal-exit-rate models such as JC69, an approximation for HKY/GTR, whose per-state exit rates make the true change count Markov-modulated rather than Poisson. Under this model the number of changes on a branch of a given duration is a random count, so a branch length inferred from an observed number of substitutions is a distribution, not a point. TreeTime propagates that uncertainty rather than discarding it. Each edge's substitution likelihood is converted to a time-domain distribution, and the backward and forward message passing that also carries the coalescent prior (Section 3.5) yields a marginal posterior time distribution at every node whose width encodes the propagated substitution stochasticity. The node-date confidence interval is the highest-posterior-density region of that marginal, computed in `packages/treetime/src/timetree/confidence.rs` [packages/treetime/src/timetree/confidence.rs](../../packages/treetime/src/timetree/confidence.rs). A second, independent source of date uncertainty -- the standard error of the clock-rate regression slope -- is combined with it in quadrature, but that is a regression effect rather than a jump-process effect and sits outside this report's scope. See [kb/algo/timetree.md](../algo/timetree.md).

---

## 5. F3 -- Marginal reconstruction as an endpoint-conditioned CTMC

### 5.1 Two-pass sum-product

Marginal ancestral reconstruction is Felsenstein pruning <a id="cite-9"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[9](#ref-9)] followed by an outward pass -- the sum-product (belief-propagation) algorithm on the tree. The upward (post-order) pass computes each node's partial likelihood, the probability of the subtree data given each ancestral state $X$: for an internal node $k$ with children $i,j$,

$$w_k(X) = \Big[\sum_{Y} P_{X\to Y}(t_i)\, w_i(Y)\Big]\Big[\sum_{Z} P_{X\to Z}(t_j)\, w_j(Z)\Big],$$

where $w_k(X)$ is the partial likelihood at node $k$ in state $X$, and $P_{X\to Y}(t)$ is the entry of $P(t)=e^{Qt}$. Tips enter as boundary conditions: a leaf's distribution is one-hot on its observed state (spread over states for ambiguity codes). At the root the site likelihood is $\sum_X \pi_X\, w_{\mathrm{root}}(X)$. The downward (pre-order) pass forms, for each node, the _cavity_ message -- the information from the rest of the tree excluding the node's own subtree, obtained by dividing out that subtree's upward message -- and multiplies it into the upward message to yield the node's marginal posterior. Orchestration is in `packages/treetime/src/ancestral/marginal.rs` [packages/treetime/src/ancestral/marginal.rs#L30-L41](../../packages/treetime/src/ancestral/marginal.rs#L30-L41); the recurrence core is `packages/treetime/src/partition/marginal_core.rs` [packages/treetime/src/partition/marginal_core.rs#L109-L269](../../packages/treetime/src/partition/marginal_core.rs#L109-L269). Both dense (all positions) and sparse (variable positions plus Fitch-compressed fixed blocks) backends run this same recurrence; the algorithm is documented in [kb/algo/ancestral.md](../algo/ancestral.md).

Joint reconstruction (the single most-likely assignment across all nodes, a Viterbi-style argmax) is intentionally removed: the pipeline rejects it at runtime in `packages/treetime/src/ancestral/pipeline.rs` [packages/treetime/src/ancestral/pipeline.rs#L225-L233](../../packages/treetime/src/ancestral/pipeline.rs#L225-L233). See [kb/decisions/ancestral-joint-reconstruction-removed.md](../decisions/ancestral-joint-reconstruction-removed.md).

### 5.2 The conditioning connection

The marginal posterior at an internal node is, by definition, a CTMC state distribution _conditioned on the observed tip states_ -- an endpoint-conditioned CTMC quantity. The theory of conditioning a Markov process on future observations is the Doob $h$-transform: for a process conditioned to reach state $z$ at time $T$, with $h(t,x)=P\big(X(T)=z \mid X(t)=x\big)$, the conditioned off-diagonal jump rates become

$$q_h(t; x, y) = q(x,y)\,\frac{h(t,y)}{h(t,x)},$$

which become time-inhomogeneous whenever $h$ depends on the remaining time to $T$, even when the original process is homogeneous <a id="cite-14"></a>[Corstanje, van der Meulen, and Schauer 2023](https://doi.org/10.1080/17442508.2022.2150081) [[14](#ref-14)]. This is why a constant-rate waiting time cannot simply be reused inside a conditioned bridge.

TreeTime computes the conditional _distributions_ (sum over paths) rather than sampling conditioned _paths_. The path-sampling counterpart is phylogenetic stochastic mapping -- drawing a random substitution history given both branch endpoints -- for which the standard methods are Nielsen's mutation mapping <a id="cite-11"></a>[Nielsen 2002](https://doi.org/10.1080/10635150290102393) [[11](#ref-11)], the endpoint-conditioned CTMC simulation methods surveyed by <a id="cite-12"></a>[Hobolth and Stone 2009](https://doi.org/10.1214/09-AOAS247) [[12](#ref-12)], and the uniformization-based mapping of <a id="cite-13"></a>[Irvahn and Minin 2014](https://doi.org/10.1089/cmb.2014.0062) [[13](#ref-13)]. TreeTime implements the deterministic (marginal-distribution) side of this problem, not the sampling side. Conditioning-on-endpoints is the shared framework; the $h$-transform is one way to describe it and does not by itself imply that any particular algorithm simulates the $h$-transformed process.

### 5.3 Expected event counts without sampling

The quantities that stochastic mapping estimates by drawing conditioned paths -- the expected number of transitions of each type along a branch, and the expected time spent in each state -- TreeTime approximates from the branch endpoint joint posterior, without sampling a single path. GTR inference consumes these endpoint-marginal statistics as its sufficient statistics: an endpoint substitution-count matrix and midpoint-rule per-state dwell times, accumulated in `packages/treetime/src/gtr/infer_gtr/common.rs` [packages/treetime/src/gtr/infer_gtr/common.rs](../../packages/treetime/src/gtr/infer_gtr/common.rs). Writing $n_{ij}$ for the branch endpoint joint mass $P(\text{parent}=j,\ \text{child}=i)$ (used as an approximate substitution count, exact only when at most one substitution falls on the branch) and $T_k$ for the approximate time in state $k$,

$$T_k \mathrel{+}= \tfrac{1}{2}\,\ell_e\,\big(P(\text{parent}=k) + P(\text{child}=k)\big),$$

summed over every edge $e$ of branch length $\ell_e$ and every site, where $P(\text{parent}=k)$ and $P(\text{child}=k)$ are the marginal posteriors at the edge endpoints. These are endpoint-marginal approximations to the path summaries produced by the stochastic-mapping methods of Section 5.2 -- $T_k$ is the trapezoidal rule for the conditional occupation integral $\int_0^{\ell_e} P\big(X(s)=k \mid \text{endpoints}\big)\,ds$, and $n_{ij}$ replaces the endpoint-conditioned expected transition count with endpoint-difference mass; both are exact only in the single-substitution-per-branch limit and carry an endpoint-approximation bias that grows with branch length and per-site rate. The GTR coordinate descent then re-estimates the exchangeabilities, equilibrium frequencies, and overall rate from $n_{ij}$ and $T_k$. So approximations to the _outputs_ of endpoint-conditioned mapping (transition counts, dwell times) are present in TreeTime even though the _sampling procedure_ -- and the exact expectations it targets -- are not. See [kb/algo/gtr.md](../algo/gtr.md).

### 5.4 Discrete-trait reconstruction (mugration)

The `mugration` command reconstructs discrete metadata -- sampling country, host species, lineage -- by treating each state as a character evolving under a GTR-like transition matrix, run through the identical two-pass pruning of Section 5.1. A leaf with a known trait enters as a one-hot boundary distribution; a missing trait (`"?"`) enters as a uniform prior that is marginalized during message passing. The alphabet changes from four nucleotides to an arbitrary discrete state set; the endpoint-conditioned CTMC framework does not. This is the implemented realization of reconstructing discrete host or geographic-state histories on a phylogeny, and it shares `packages/treetime/src/partition/marginal_core.rs` [packages/treetime/src/partition/marginal_core.rs](../../packages/treetime/src/partition/marginal_core.rs) with nucleotide reconstruction through the discrete partition in `packages/treetime/src/partition/marginal_discrete.rs` [packages/treetime/src/partition/marginal_discrete.rs](../../packages/treetime/src/partition/marginal_discrete.rs). See [kb/algo/mugration.md](../algo/mugration.md).

### 5.5 Method note: eigendecomposition, not uniformization

For the endpoint problem, the stochastic-mapping literature often reaches for uniformization -- representing the CTMC as a Poisson-distributed number of steps of an embedded discrete-time chain -- because it samples conditioned paths without any matrix exponential (the methods of Section 5.2). TreeTime takes the complementary route: it forms $P(t)=e^{Qt}$ once per model by eigendecomposition of the symmetrized reversible generator (Section 4.1) and reuses that factorization across all branches, then obtains marginals by pruning. The choice follows from the target. TreeTime needs endpoint likelihoods and marginal distributions, for which an amortized per-model eigendecomposition is a direct fit; uniformization's advantages -- avoiding matrix exponentials and remaining robust for stiff or non-diagonalizable generators -- apply to the path-sampling problem TreeTime does not solve. The eigendecomposition and per-branch $P(t)$ live in `fn Gtr.eig_single_site()` [packages/treetime/src/gtr/gtr.rs#L74-L105](../../packages/treetime/src/gtr/gtr.rs#L74-L105) and `fn Gtr.expQt()` [packages/treetime/src/gtr/gtr.rs#L453-L463](../../packages/treetime/src/gtr/gtr.rs#L453-L463).

---

## 6. Directions for future work

The following extensions are not implemented in v1. They range from a v0 feature not yet ported (reassortment-aware dating, Section 6.1) through models representable in the current data model but not inferred (recombination, Section 6.1) to methods the project has never modeled (birth-death priors, Section 6.4; scalable ARG inference, Section 6.5; simulation-based neural inference, Section 6.7). Each subsection states the basis in the codebase and the corresponding prior art.

The ordering reflects demand rather than incremental distance from the current code: the directions that turn expert-only, BEAST2-dominated capabilities into fast, pipeline-integrated defaults are listed first, because this is the pattern TreeTime has already applied twice -- to molecular-clock dating and to geographic reconstruction via mugration.

### 6.1 Reassortment and recombination

**Reassortment (v0 parity gap).** TreeTime v0 has an `arg` command that consumes `TreeKnit` segment trees and their maximally-compatible-clade (MCC) assignments <a id="cite-15"></a>[Barrat-Charlaix et al. 2022](https://doi.org/10.1371/journal.pcbi.1010394) [[15](#ref-15)], masks each branch to the combined or segment-specific alignment according to whether parent and child share an MCC, and infers reassortment-aware time trees. That command is not yet ported to v1 (see [kb/features/arg.md](../features/arg.md)). The v0 implementation is tree-based: it runs two separate segment trees with per-branch alignment masking and does not require multi-parent nodes. Porting it to v1 reuses the existing partition and timetree infrastructure with no graph changes. Full Bayesian reassortment-network inference, in which the network itself is the inferred object, is a separate, larger target <a id="cite-16"></a>[Muller et al. 2020](https://doi.org/10.1073/pnas.1918304117) [[16](#ref-16)].

**Recombination (representable, not inferred).** TreeTime's directed graph admits multiple-parent nodes and first-class edges, a choice made to host reticulate ancestry (see [kb/decisions/graph-based-phylogenetic-representation.md](../decisions/graph-based-phylogenetic-representation.md)). The inference algorithms require a tree: `fn get_exactly_one_parent()` [packages/treetime-graph/src/graph_traverse.rs#L84-L87](../../packages/treetime-graph/src/graph_traverse.rs#L84-L87) and `fn get_exactly_one_parent_edge()` [packages/treetime-graph/src/graph_traverse.rs#L160-L161](../../packages/treetime-graph/src/graph_traverse.rs#L160-L161) reject multi-parent nodes, and the marginal core assumes one parent edge [packages/treetime/src/partition/marginal_core.rs#L165-L167](../../packages/treetime/src/partition/marginal_core.rs#L165-L167). A recombination extension requires two additions: genomic-interval annotations on edges, for which the ARG and tree-sequence data models of `msprime` <a id="cite-17"></a>[Baumdicker et al. 2022](https://doi.org/10.1093/genetics/iyab229) [[17](#ref-17)] and <a id="cite-18"></a>[Wong et al. 2024](https://doi.org/10.1093/genetics/iyae100) [[18](#ref-18)] define the structure, and a likelihood on the network, for which the coalescent with recombination of <a id="cite-19"></a>[Hudson 1983](<https://doi.org/10.1016/0040-5809(83)90013-8>) [[19](#ref-19)] specifies the event process. Bayesian inference on bacterial ARGs is implemented in `Bacter` <a id="cite-20"></a>[Vaughan et al. 2017](https://doi.org/10.1534/genetics.116.193425) [[20](#ref-20)].

For bacteria the relevant process is gene-conversion-like tract import rather than crossover. The field-specific tools are `ClonalOrigin` <a id="cite-21"></a>[Didelot et al. 2010](https://doi.org/10.1534/genetics.110.120121) [[21](#ref-21)], `SimBac` <a id="cite-22"></a>[Brown et al. 2016](https://doi.org/10.1099/mgen.0.000044) [[22](#ref-22)], and the Bacterial Sequential Markov Coalescent <a id="cite-23"></a>[De Maio and Wilson 2017](https://doi.org/10.1534/genetics.116.198796) [[23](#ref-23)].

### 6.2 Birth-death-sampling tree prior and effective reproduction number

The coalescent is TreeTime's only node-time prior. The birth-death-sampling process is the main alternative, modeling transmission, becoming-noninfectious, and sampling directly at the lineage level; its skyline form estimates piecewise-constant epidemiological rates and has been applied to HIV and hepatitis C <a id="cite-24"></a>[Stadler et al. 2013](https://doi.org/10.1073/pnas.1207965110) [[24](#ref-24)]. A birth-death-sampling contribution would enter as an alternative to the coalescent term of Section 3, reusing the same node-time belief propagation. Its principal output -- the effective reproduction number $R_e(t)$ -- is currently available only through BEAST2 BDSKY, which requires MCMC. A fast ML birth-death-skyline emitting $R_e(t)$ alongside the dated tree would bring this quantity to the same level of accessibility as node dates.

A calibration requirement applies: birth-death rate functions can have entire congruence classes with identical tree likelihoods <a id="cite-25"></a>[Louca and Pennell 2020](https://doi.org/10.1038/s41586-020-2176-1) [[25](#ref-25)], so any $R_e(t)$ output must ship identifiability-aware intervals, not bare point estimates.

### 6.3 Structured coalescent for transmission

The mugration command (Section 5.4) treats geography or host as a discrete character evolving on a fixed tree. This is a mutational approximation to migration: it does not let population structure shape the genealogy. The structured coalescent is the model that does, with lineages coalescing within demes and migrating between them, and hosts acting as demes for transmission reconstruction, as implemented in `SCOTTI` <a id="cite-26"></a>[De Maio et al. 2016](https://doi.org/10.1371/journal.pcbi.1005130) [[26](#ref-26)]. Because mugration already runs discrete-state marginal reconstruction through the shared core, it is the point from which a structured-coalescent prior would be introduced. The structured coalescent couples state and genealogy in a way the current discrete-character model does not: lineage-specific coalescent rates depend on the deme, and migration changes the deme assignment backward in time. An analytic or linearized approximation, rather than full MCMC over the structured process, would preserve TreeTime's ML character.

### 6.4 Multiple-merger coalescents

The coalescent prior and the polytomy resolver both assume pairwise mergers. A $\Lambda$-coalescent allows one event to merge $k \ge 2$ of the $b$ active lineages with total rate $a_k(b)=\binom{b}{k}\,\lambda_{b,k}$ <a id="cite-27"></a>[Pitman 1999](https://doi.org/10.1214/aop/1022874819) [[27](#ref-27)]; the Beta-coalescent family of $\Lambda$ measures is the most studied parametric case <a id="cite-28"></a>[Schweinsberg 2003](<https://doi.org/10.1016/S0304-4149(03)00028-0>) [[28](#ref-28)]; and a $\Xi$-coalescent further permits several groups to merge in one event <a id="cite-29"></a>[Schweinsberg 2000](https://doi.org/10.1214/ejp.v5-68) [[29](#ref-29)]. Eldon and Wakeley introduced an explicit model of reproductive sweepstakes that produces multiple-merger genealogies in marine organisms, with direct relevance to pathogen superspreading events <a id="cite-30"></a>[Eldon and Wakeley 2006](https://doi.org/10.1534/genetics.105.052175) [[30](#ref-30)].

The motivation is specific to pathogens. Superspreading, transmission bottlenecks, and rapid clonal expansion produce skewed offspring distributions that violate the Kingman assumptions; Kingman-based inference then mis-estimates demography and selection, as shown for influenza A <a id="cite-31"></a>[Sackman et al. 2019](https://doi.org/10.1534/genetics.118.301684) [[31](#ref-31)], and multiple-merger genealogies fit _Mycobacterium tuberculosis_ outbreaks better than Kingman genealogies <a id="cite-32"></a>[Menardo et al. 2021](https://doi.org/10.1093/molbev/msaa179) [[32](#ref-32)].

A $\Lambda$-coalescent requires a deeper change to the coalescent objective than a rate-table swap. Under Kingman, the per-lineage rate $\kappa = (k-1)/(2T_c)$ satisfies $k \cdot \kappa = \lambda = \binom{k}{2}/T_c$, so the cumulative hazard $H(t)=\int \kappa\,ds$ is a per-edge sufficient statistic, and an $m$-furcation decomposes into $(m-1)$ pairwise mergers. This is the structure `packages/treetime/src/coalescent/integration.rs` and the node contribution $(m-1)(H - \ln\lambda)$ in `packages/treetime/src/coalescent/coalescent.rs` encode. Under a general $\Lambda$, the total rate $\lambda_b = \sum_{k \ge 2}\binom{b}{k}\lambda_{b,k}$ is not a sum of independent per-lineage hazards, so survival must be computed per interval as $\exp(-\lambda_b \Delta t)$, and an $m$-merger is a single event with density $\lambda_{b,m}$, not $(m-1)$ pairwise terms. Both the edge-factored hazard decomposition and the $(m-1)$ multiplicity break, and the closed-form skyline $T_c^* = I/M$ of Section 3.3 is no longer available. The $(\alpha, T_c)$ parameter pair under a Beta-coalescent is weakly identifiable from a single tree, so the recommended approach is to profile $\alpha$ on a grid or fix it from external evidence, not to estimate it jointly.

### 6.5 Scalable ARG inference and external-genealogy dating

The landscape of ancestral recombination graph inference has expanded rapidly since 2019. For TreeTime's purposes three developments are relevant.

First, **scalable heuristic ARG construction** -- `tsinfer` <a id="cite-33"></a>[Kelleher et al. 2019](https://doi.org/10.1038/s41588-019-0483-y) [[33](#ref-33)], `Relate` <a id="cite-34"></a>[Speidel et al. 2019](https://doi.org/10.1038/s41588-019-0484-x) [[34](#ref-34)], and `ARG-Needle` <a id="cite-35"></a>[Zhang et al. 2023](https://doi.org/10.1038/s41588-023-01379-x) [[35](#ref-35)] -- now produces point-estimate genealogies at biobank scale (10^5--10^6 genomes). These are primarily human-genomics tools, but the two-step pattern they define -- construct an ARG externally, then date/annotate it separately -- is exactly the pattern TreeTime already follows for reassortment (v0 `arg` + TreeKnit) and for UShER mutation-annotated trees (v1 already ingests UShER MAT via `packages/treetime-io/src/usher_mat.rs`). An official `tskit` Rust crate (v0.16.3, `tskit-dev/tskit-rust`) exists, so ingesting a tree-sequence into TreeTime's `Graph` is a tractable I/O step.

Second, **Bayesian ARG sampling with uncertainty** -- `SINGER` <a id="cite-36"></a>[Deng, Nielsen, and Song 2025](https://doi.org/10.1038/s41588-025-02317-9) [[36](#ref-36)] accelerated ARG posterior sampling by two orders of magnitude over ARGweaver, making full posterior exploration feasible for hundreds of genomes. This enables genuine uncertainty quantification over the genealogy, not just a point estimate.

Third, **deep-learning methods operating on ARGs** -- SIA uses LSTM networks on ARG-derived features to infer selection coefficients and allele-frequency trajectories <a id="cite-37"></a>[Hejase et al. 2022](https://doi.org/10.1093/molbev/msab332) [[37](#ref-37)], and PhyloDeep uses simulation-based deep learning on tree shapes to estimate epidemiological parameters <a id="cite-38"></a>[Voznica et al. 2022](https://doi.org/10.1038/s41467-022-31511-0) [[38](#ref-38)]. These are amortized estimators: expensive simulation happens once during training, and each subsequent inference is a forward pass.

A correctness gate applies to any external-genealogy dating in TreeTime: independent per-local-tree dating double-counts shared ARG nodes. Correct dating requires ARG-node-centric time variables that aggregate edge likelihoods across the local trees in which each node participates -- a genuine re-architecture of the message-passing system, not a simple iteration over local trees.

### 6.6 Topology uncertainty

TreeTime's node-date confidence (Section 4.3) conditions on a single fixed input tree and ignores topology uncertainty. For pathogen trees with low per-branch mutation counts, topology uncertainty can dominate branch-length uncertainty. The existing rate-susceptibility machinery (`packages/treetime/src/timetree/confidence.rs`) already re-runs `run_timetree` at perturbed clock rates and aggregates per-node dates by `GraphNodeKey` -- the same pattern generalized to an input set of alternative topologies (bootstrap or near-optimal ML trees) yields per-clade date intervals that incorporate topological variation. The output is a topological sensitivity band, reported alongside clade support values; it is not a calibrated posterior interval, and should not be labeled as one. For recombining organisms, posterior tree samples from SINGER (Section 6.5) would be the appropriate input; for the non-recombining single-segment case that dominates TreeTime's current use, bootstrap or UFBoot trees are the natural source.

### 6.7 Simulation-based amortized inference

For models with no tractable tree likelihood -- many-deme structured coalescent, reassortment networks, selection, combined multiple-merger and population structure -- the analytic likelihood approach that serves TreeTime for Kingman $T_c$ is not available. Simulation-based inference (SBI) trains a neural estimator on a large set of simulated genealogies and uses the trained network for fast amortized parameter estimation <a id="cite-39"></a>[Cranmer, Brehmer, and Louppe 2020](https://doi.org/10.1073/pnas.1912789117) [[39](#ref-39)]. Pathogen-scale trees (10^2--10^3 tips) are the regime where PhyloDeep <a id="cite-38">[Voznica et al. 2022](https://doi.org/10.1038/s41467-022-31511-0)</a> [[38](#ref-38)] and SIA <a id="cite-37">[Hejase et al. 2022](https://doi.org/10.1093/molbev/msab332)</a> [[37](#ref-37)] already operate.

An SBI estimator for TreeTime parameters must never replace the analytic path for quantities TreeTime already computes exactly (clock rate, Kingman $T_c$ via the closed-form skyline). It applies only where no analytic likelihood exists. The simulator must include a sampling process (pathogen tips are preferentially sampled; ignoring this produces confidently wrong estimates on real data). In-distribution calibration -- the standard SBI quality check -- cannot detect structural non-identifiability such as the birth-death congruence classes of <a id="cite-25a"></a>[Louca and Pennell 2020](https://doi.org/10.1038/s41586-020-2176-1) [[25](#ref-25)]; an identifiability analysis of the target model before training is required.

Forward pathogen simulation (`VGsim` <a id="cite-40"></a>[Shchur et al. 2022](https://doi.org/10.1371/journal.pcbi.1010409) [[40](#ref-40)], `phastSim` Section 4.2) is separate from inference and shares the direct-method sampler of Section 2; its relevance to TreeTime is as a training-data source for the simulation-based estimators of Section 6.7.

---

## 7. Summary: one construction, three uses

A single construction -- total hazard, exponential waiting time, rate-proportional event choice, state update, with correct handling of deterministic boundaries -- appears in TreeTime in one _simulated_ form and two _analytic_ forms:

- **Simulated**: the polytomy resolver runs it forward as a labelled-forest event process (Section 2).
- **Analytic, coalescent**: the same rate structure is integrated to a closed-form log-prior on node times (Section 3).
- **Analytic, substitution**: the same CTMC structure is exponentiated per branch and summed by pruning (Sections 4-5).

The design consequence, visible in the code layout, is that the reusable parts (hazard integration, waiting-time sampling under piecewise-constant rates, deterministic-boundary handling, reproducible random streams) are isolated in the polytomy sweep, while the coalescent and substitution paths reuse only the _rate definitions_ and integrate them analytically. This is why the coalescent's $\kappa(t)$ schedule can be shared by the resolver without either path simulating the other's process.

## 8. Gaps and open questions

- **Skyline hazard quadrature.** The cumulative coalescent hazard $H(t)$ is integrated by midpoint quadrature over lineage-count breakpoints. A skyline $T_c$ discontinuity that falls strictly inside a lineage-count interval is approximated by the interval midpoint rather than represented as its own breakpoint. Whether this is exact enough against the approved coalescent target is an open verification item; it is a candidate correctness question, not an established defect, and should be checked against the coalescent specification and the skyline quadrature-contract issues already tracked in [kb/algo/coalescent-contribution-refactor.md](../algo/coalescent-contribution-refactor.md).
- **Coalescent rate in the resolver.** The resolver's coalescent channel uses $R_{\mathrm{coal}}=(|R|-1)\,\kappa(t)$, where $|R|$ is the count of mutation-free lineages _within the polytomy_ but $\kappa(t)=(k_{\mathrm{global}}-1)/(2T_c)$ is the _whole-tree_ per-lineage schedule shared with Section 3. The effective total is therefore $(|R|-1)(k_{\mathrm{global}}-1)/(2T_c)$ -- a local count times a global per-lineage rate, neither the local all-pairs $\binom{|R|}{2}/T_c$ nor a purely local per-lineage law, though consistent with the documented first-merger mean $1/((|R|-1)\kappa)$. This is the implemented and tested behaviour; whether it is the intended coalescent law for resolution is a modelling question for the maintainers.

---

## Glossary

1. <a id="gloss-1"></a> **Gillespie direct method.** Exact event-driven sampler for a continuous-time Markov jump process: compute the total event rate, draw an exponential waiting time, choose an event in proportion to its rate, update the state, repeat ([Gillespie 1977](https://doi.org/10.1021/j100540a008) [[3](#ref-3)]). [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **CTMC (continuous-time Markov chain).** A Markov process on a discrete state space that holds in each state for an exponential time and then jumps; defined by a rate matrix $Q$ whose branch transition matrix is $P(t)=e^{Qt}$. [↩](#gloss-use-2)
3. <a id="gloss-3"></a> **Felsenstein pruning.** Post-order dynamic program that computes the likelihood of tip data given a tree and substitution model by propagating per-state partial likelihoods from leaves to root ([Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[9](#ref-9)]); equivalent to the sum-product algorithm on a tree. [↩](#gloss-use-3)
4. <a id="gloss-4"></a> **Endpoint-conditioned CTMC.** A CTMC whose law is conditioned on the states at both ends of an interval (e.g. observed tips); its effective jump rates become time-inhomogeneous. [↩](#gloss-use-4)
5. <a id="gloss-5"></a> **Stochastic mapping.** Sampling a full substitution history along a branch conditional on its endpoint states -- the sampling counterpart of marginal reconstruction ([Nielsen 2002](https://doi.org/10.1080/10635150290102393) [[11](#ref-11)]). [↩](#gloss-use-5)
6. <a id="gloss-6"></a> **Doob $h$-transform.** A construction that conditions a Markov process on a future event using a (space-time) harmonic function $h$, reweighting the jump rates by $h(t,y)/h(t,x)$ ([Corstanje et al. 2023](https://doi.org/10.1080/17442508.2022.2150081) [[14](#ref-14)]). [↩](#gloss-use-6)
7. <a id="gloss-7"></a> **Propensity.** The instantaneous rate $a_j(x)$ of event channel $j$ given the current state $x$; the probability of one firing of $j$ in $[t,t+\mathrm{d}t)$ is $a_j(x)\,\mathrm{d}t + o(\mathrm{d}t)$. [↩](#gloss-use-7)
8. <a id="gloss-8"></a> **Memorylessness.** The defining property of the exponential distribution: the residual waiting time is independent of elapsed time, which makes discarding-and-redrawing an event under a changed rate unbiased. [↩](#gloss-use-8)

## References

1. <a id="ref-1"></a> Doob, Joseph L. 1945. "Markoff Chains--Denumerable Case." _Transactions of the American Mathematical Society_ 58(3):455-473. https://doi.org/10.2307/1990339 [↩](#cite-1)
2. <a id="ref-2"></a> Gillespie, Daniel T. 1976. "A General Method for Numerically Simulating the Stochastic Time Evolution of Coupled Chemical Reactions." _Journal of Computational Physics_ 22(4):403-434. https://doi.org/10.1016/0021-9991(76)90041-3 [↩](#cite-2)
3. <a id="ref-3"></a> Gillespie, Daniel T. 1977. "Exact Stochastic Simulation of Coupled Chemical Reactions." _The Journal of Physical Chemistry_ 81(25):2340-2361. https://doi.org/10.1021/j100540a008 [↩](#cite-3)
4. <a id="ref-4"></a> Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-4)
5. <a id="ref-5"></a> Kingman, John F. C. 1982. "The Coalescent." _Stochastic Processes and their Applications_ 13(3):235-248. https://doi.org/10.1016/0304-4149(82)90011-4 [↩](#cite-5)
6. <a id="ref-6"></a> Pybus, Oliver G., Andrew Rambaut, and Paul H. Harvey. 2000. "An Integrated Framework for the Inference of Viral Population History from Reconstructed Genealogies." _Genetics_ 155(3):1429-1437. https://doi.org/10.1093/genetics/155.3.1429 [↩](#cite-6)
7. <a id="ref-7"></a> Drummond, Alexei J., Andrew Rambaut, Beth Shapiro, and Oliver G. Pybus. 2005. "Bayesian Coalescent Inference of Past Population Dynamics from Molecular Sequences." _Molecular Biology and Evolution_ 22(5):1185-1192. https://doi.org/10.1093/molbev/msi103 [↩](#cite-7)
8. <a id="ref-8"></a> Tavaré, Simon. 1986. "Some Probabilistic and Statistical Problems in the Analysis of DNA Sequences." _Lectures on Mathematics in the Life Sciences_ 17:57-86. https://openalex.org/W1593676244 [↩](#cite-8)
9. <a id="ref-9"></a> Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-9)
10. <a id="ref-10"></a> De Maio, Nicola, William Boulton, Lukas Weilguny, Conor R. Walker, Yatish Turakhia, Russell Corbett-Detig, and Nick Goldman. 2022. "phastSim: Efficient Simulation of Sequence Evolution for Pandemic-Scale Datasets." _PLOS Computational Biology_ 18(4):e1010056. https://doi.org/10.1371/journal.pcbi.1010056 [↩](#cite-10)
11. <a id="ref-11"></a> Nielsen, Rasmus. 2002. "Mapping Mutations on Phylogenies." _Systematic Biology_ 51(5):729-739. https://doi.org/10.1080/10635150290102393 [↩](#cite-11)
12. <a id="ref-12"></a> Hobolth, Asger, and Eric A. Stone. 2009. "Simulation from Endpoint-Conditioned, Continuous-Time Markov Chains on a Finite State Space, with Applications to Molecular Evolution." _The Annals of Applied Statistics_ 3(3):1204-1231. https://doi.org/10.1214/09-AOAS247 [↩](#cite-12)
13. <a id="ref-13"></a> Irvahn, Jan, and Vladimir N. Minin. 2014. "Phylogenetic Stochastic Mapping without Matrix Exponentiation." _Journal of Computational Biology_ 21(9):676-690. https://doi.org/10.1089/cmb.2014.0062 [↩](#cite-13)
14. <a id="ref-14"></a> Corstanje, Marc, Frank van der Meulen, and Moritz Schauer. 2023. "Conditioning Continuous-Time Markov Processes by Guiding." _Stochastics_ 95(6):963-996. https://doi.org/10.1080/17442508.2022.2150081 [↩](#cite-14)
15. <a id="ref-15"></a> Barrat-Charlaix, Pierre, Timothy G. Vaughan, and Richard A. Neher. 2022. "TreeKnit: Inferring Ancestral Reassortment Graphs of Influenza Viruses." _PLOS Computational Biology_ 18(8):e1010394. https://doi.org/10.1371/journal.pcbi.1010394 [↩](#cite-15)
16. <a id="ref-16"></a> Muller, Nicola F., Ugne Stolz, Gytis Dudas, Tanja Stadler, and Timothy G. Vaughan. 2020. "Bayesian Inference of Reassortment Networks Reveals Fitness Benefits of Reassortment in Human Influenza Viruses." _Proceedings of the National Academy of Sciences_ 117(29):17104-17111. https://doi.org/10.1073/pnas.1918304117 [↩](#cite-16)
17. <a id="ref-17"></a> Baumdicker, Franz, Gertjan Bisschop, Matthew Goldstein, et al. 2022. "Efficient Ancestry and Mutation Simulation with msprime 1.0." _Genetics_ 220(3):iyab229. https://doi.org/10.1093/genetics/iyab229 [↩](#cite-17)
18. <a id="ref-18"></a> Wong, Yan, Anastasia Ignatieva, Jere Koskela, Gregor Gorjanc, Anthony W. Wohns, and Jerome Kelleher. 2024. "A General and Efficient Representation of Ancestral Recombination Graphs." _Genetics_ 228(1):iyae100. https://doi.org/10.1093/genetics/iyae100 [↩](#cite-18)
19. <a id="ref-19"></a> Hudson, Richard R. 1983. "Properties of a Neutral Allele Model with Intragenic Recombination." _Theoretical Population Biology_ 23(2):183-201. https://doi.org/10.1016/0040-5809(83)90013-8 [↩](#cite-19)
20. <a id="ref-20"></a> Vaughan, Timothy G., David Welch, Alexei J. Drummond, Patrick J. Biggs, Tessy George, and Nigel P. French. 2017. "Inferring Ancestral Recombination Graphs from Bacterial Genomic Data." _Genetics_ 205(2):857-870. https://doi.org/10.1534/genetics.116.193425 [↩](#cite-20)
21. <a id="ref-21"></a> Didelot, Xavier, Daniel Lawson, Aaron Darling, and Daniel Falush. 2010. "Inference of Homologous Recombination in Bacteria Using Whole-Genome Sequences." _Genetics_ 186(4):1435-1449. https://doi.org/10.1534/genetics.110.120121 [↩](#cite-21)
22. <a id="ref-22"></a> Brown, Thomas, Xavier Didelot, Daniel J. Wilson, and Nicola De Maio. 2016. "SimBac: Simulation of Whole Bacterial Genomes with Homologous Recombination." _Microbial Genomics_ 2(1):e000044. https://doi.org/10.1099/mgen.0.000044 [↩](#cite-22)
23. <a id="ref-23"></a> De Maio, Nicola, and Daniel J. Wilson. 2017. "The Bacterial Sequential Markov Coalescent." _Genetics_ 206(1):333-343. https://doi.org/10.1534/genetics.116.198796 [↩](#cite-23)
24. <a id="ref-24"></a> Stadler, Tanja, Denise Kuhnert, Sebastian Bonhoeffer, and Alexei J. Drummond. 2013. "Birth-Death Skyline Plot Reveals Temporal Changes of Epidemic Spread in HIV and Hepatitis C Virus (HCV)." _Proceedings of the National Academy of Sciences_ 110(1):228-233. https://doi.org/10.1073/pnas.1207965110 [↩](#cite-24)
25. <a id="ref-25"></a> Louca, Stilianos, and Matthew W. Pennell. 2020. "Extant Timetrees Are Consistent with a Myriad of Diversification Histories." _Nature_ 580:502-505. https://doi.org/10.1038/s41586-020-2176-1 [↩¹](#cite-25) [↩²](#cite-25a)
26. <a id="ref-26"></a> De Maio, Nicola, Chieh-Hsi Wu, and Daniel J. Wilson. 2016. "SCOTTI: Efficient Reconstruction of Transmission within Outbreaks with the Structured Coalescent." _PLOS Computational Biology_ 12(9):e1005130. https://doi.org/10.1371/journal.pcbi.1005130 [↩](#cite-26)
27. <a id="ref-27"></a> Pitman, Jim. 1999. "Coalescents with Multiple Collisions." _The Annals of Probability_ 27(4):1870-1902. https://doi.org/10.1214/aop/1022874819 [↩](#cite-27)
28. <a id="ref-28"></a> Schweinsberg, Jason. 2003. "Coalescent Processes Obtained from Supercritical Galton-Watson Processes." _Stochastic Processes and their Applications_ 106(1):107-139. https://doi.org/10.1016/S0304-4149(03)00028-0 [↩](#cite-28)
29. <a id="ref-29"></a> Schweinsberg, Jason. 2000. "Coalescents with Simultaneous Multiple Collisions." _Electronic Journal of Probability_ 5:1-50. https://doi.org/10.1214/ejp.v5-68 [↩](#cite-29)
30. <a id="ref-30"></a> Eldon, Bjarki, and John Wakeley. 2006. "Coalescent Processes When the Distribution of Offspring Number Among Individuals Is Highly Skewed." _Genetics_ 172(4):2621-2633. https://doi.org/10.1534/genetics.105.052175 [↩](#cite-30)
31. <a id="ref-31"></a> Sackman, Andrew M., Rebecca B. Harris, and Jeffrey D. Jensen. 2019. "Inferring Demography and Selection in Organisms Characterized by Skewed Offspring Distributions." _Genetics_ 211(3):1019-1028. https://doi.org/10.1534/genetics.118.301684 [↩](#cite-31)
32. <a id="ref-32"></a> Menardo, Fabrizio, Sebastian Gagneux, and Fabian Freund. 2021. "Multiple Merger Genealogies in Outbreaks of _Mycobacterium tuberculosis_." _Molecular Biology and Evolution_ 38(1):290-306. https://doi.org/10.1093/molbev/msaa179 [↩](#cite-32)
33. <a id="ref-33"></a> Kelleher, Jerome, Yan Wong, Anthony W. Wohns, Chaimaa Fadil, Patrick K. Albers, and Gil McVean. 2019. "Inferring Whole-Genome Histories in Large Population Datasets." _Nature Genetics_ 51(9):1330-1338. https://doi.org/10.1038/s41588-019-0483-y [↩](#cite-33)
34. <a id="ref-34"></a> Speidel, Leo, Marie Forest, Sinan Shi, and Simon R. Myers. 2019. "A Method for Genome-Wide Genealogy Estimation for Thousands of Samples." _Nature Genetics_ 51(9):1321-1329. https://doi.org/10.1038/s41588-019-0484-x [↩](#cite-34)
35. <a id="ref-35"></a> Zhang, Brian C., Arjun Biddanda, Pier F. Gunnarsson, Fergus Cooper, and Pier Francesco Palamara. 2023. "Biobank-Scale Inference of Ancestral Recombination Graphs Enables Genealogical Analysis of Complex Traits." _Nature Genetics_ 55:768-776. https://doi.org/10.1038/s41588-023-01379-x [↩](#cite-35)
36. <a id="ref-36"></a> Deng, Yun, Rasmus Nielsen, and Yun S. Song. 2025. "Robust and Accurate Bayesian Inference of Genome-Wide Genealogies for Hundreds of Genomes." _Nature Genetics_. https://doi.org/10.1038/s41588-025-02317-9 [↩](#cite-36)
37. <a id="ref-37"></a> Hejase, Hussein A., Ziyi Mo, Leonardo Campagna, and Adam Siepel. 2022. "A Deep-Learning Approach for Inference of Selective Sweeps from the Ancestral Recombination Graph." _Molecular Biology and Evolution_ 39(1):msab332. https://doi.org/10.1093/molbev/msab332 [↩](#cite-37)
38. <a id="ref-38"></a> Voznica, Jakub, Alix Zhukova, Veronika Boskova, Emma Saulnier, Frederic Lemoine, Marie Moslonka-Lefebvre, and Olivier Gascuel. 2022. "Deep Learning from Phylogenies to Uncover the Epidemiological Dynamics of Outbreaks." _Nature Communications_ 13:3896. https://doi.org/10.1038/s41467-022-31511-0 [↩](#cite-38)
39. <a id="ref-39"></a> Cranmer, Kyle, Johann Brehmer, and Gilles Louppe. 2020. "The Frontier of Simulation-Based Inference." _Proceedings of the National Academy of Sciences_ 117(48):30055-30062. https://doi.org/10.1073/pnas.1912789117 [↩](#cite-39)
40. <a id="ref-40"></a> Shchur, Vladimir, Vadim Spirin, Dmitry Sirotkin, Evgeni Burovski, Nicola De Maio, and Russell Corbett-Detig. 2022. "VGsim: Scalable Viral Genealogy Simulator for Global Pandemic." _PLOS Computational Biology_ 18(8):e1010409. https://doi.org/10.1371/journal.pcbi.1010409 [↩](#cite-40)
