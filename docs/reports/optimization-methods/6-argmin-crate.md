# Chapter 6: The argmin crate

[Back to index](_index.md) | Previous: [Chapter 5: Supporting optimizations](5-supporting-optimizations.md) | Next: [Chapter 7: Audit](7-audit.md)

TreeTime v1 uses [argmin](https://docs.rs/argmin/0.10.0/argmin/) as its optimization framework. This chapter documents argmin's architecture, how to write custom solvers, and practical patterns for this project.

## Solver inventory

| Solver                | Used in v1? | Potential use                                       |
| :-------------------- | :---------: | :-------------------------------------------------- |
| `BrentOpt`            | E1, E4, E6  | P5 (grid fallback replacement)                      |
| `BrentRoot`           |     No      | P2 (HPD root-finding)                               |
| `GoldenSectionSearch` |     E2      | Already used                                        |
| `NelderMead`          |     E5      | Already used for skyline                            |
| `Newton`              |     No      | P1 alt (requires ndarray vectors, not ideal for 1D) |
| `LBFGS`               |     No      | P8 (skyline upgrade)                                |
| `BFGS`                |     No      | Not needed (LBFGS is the limited-memory variant)    |
| `SteepestDescent`     |     No      | Not competitive                                     |
| `SimulatedAnnealing`  |     No      | Not needed (problems are well-behaved)              |
| `ParticleSwarm`       |     No      | Not needed                                          |
| `GaussNewton`         |     No      | Not applicable (not least-squares)                  |

## Core traits

For 1D problems, use `f64` directly as the parameter type:

```rust
impl CostFunction for MyProblem {
    type Param = f64;
    type Output = f64;
    fn cost(&self, x: &Self::Param) -> Result<Self::Output, Error> { ... }
}

impl Gradient for MyProblem {
    type Param = f64;
    type Gradient = f64;
    fn gradient(&self, x: &Self::Param) -> Result<Self::Gradient, Error> { ... }
}

impl Hessian for MyProblem {
    type Param = f64;
    type Hessian = f64;
    fn hessian(&self, x: &Self::Param) -> Result<Self::Hessian, Error> { ... }
}
```

`BrentOpt` and `BrentRoot` need only `CostFunction`. `LBFGS` needs `CostFunction + Gradient`. `Newton` needs all three.

`CostFunction::cost(&self)` is immutable by design. Use `RefCell` for mutable internal state.

## Custom solver pattern (for P1)

The `Solver` trait enables wrapping TreeTime's Newton-Raphson as an argmin solver:

```rust
pub struct NewtonBranchLength {
    max_iter: u64,
    relative_tolerance: f64,
}

impl Solver<BranchLengthProblem, IterState<f64, (), (), (), (), f64>>
    for NewtonBranchLength
{
    const NAME: &'static str = "Newton1D";

    fn next_iter(
        &mut self,
        problem: &mut Problem<BranchLengthProblem>,
        state: IterState<f64, (), (), (), (), f64>,
    ) -> Result<(IterState<f64, (), (), (), (), f64>, Option<KV>), Error> {
        let t = state.get_param().unwrap();
        let cost = problem.cost(t)?;
        let grad = problem.gradient(t)?;
        let hess = problem.hessian(t)?;
        // Newton step with clamping...
        let new_t = t - (grad / hess).clamp(-1.0, *t);
        Ok((state.param(new_t).cost(/* new cost */), None))
    }
}
```

Reference implementation: `crates/argmin/src/solver/newton/newton_method.rs` in the [argmin repo](https://github.com/argmin-rs/argmin).

## LBFGS configuration (for P8)

```rust
let linesearch = MoreThuenteLineSearch::new().with_c(1e-4, 0.9)?;
let solver = LBFGS::new(linesearch, 7)    // 7 correction pairs
    .with_tolerance_grad(1e-8)?
    .with_tolerance_cost(1e-12)?;
let result = Executor::new(problem, solver)
    .configure(|state| state.param(init_param).max_iters(200))
    .run()?;
```

Memory parameter `m` (number of correction pairs): typical values 5-20, 7 is common. For the skyline's 10 dimensions, this is negligible overhead.

## BrentRoot vs BrentOpt

| Aspect         | `BrentRoot`                                | `BrentOpt`                               |
| :------------- | :----------------------------------------- | :--------------------------------------- |
| Purpose        | Find x where f(x) = 0                      | Find x that minimizes f(x)               |
| Algorithm      | Bisection + secant + IQI                   | Parabolic interpolation + golden-section |
| Bracket        | Must bracket the root (f(min)\*f(max) < 0) | Must bracket the minimum                 |
| Required trait | `CostFunction`                             | `CostFunction`                           |

Both are derivative-free and work with `f64` parameters.

## Observer system

The `Observe<I>` trait has three optional methods: `observe_init`, `observe_iter`, `observe_final`. Attach via `.add_observer(observer, ObserverMode::Every(10))`.

v1 has 5 near-identical observer implementations. A shared observer would eliminate this duplication (see [audit proposal P3](7-audit.md)):

```rust
pub struct PeriodicDebugObserver {
    pub log_every: u64,
    pub label: &'static str,
}
```

## Practical notes

- ndarray feature flags must match exactly: `argmin-math = { features = ["ndarray_v0_16"] }` for ndarray 0.16
- Timer broken by default in argmin 0.11.0 - must call `.timer(true)` explicitly
- `ParticleSwarm` panics with scalar `Param = f64`
- `NelderMead` can hang in edge cases
- Disable Ctrl+C handling (`.ctrlc(false)`) for nested executors

## Alternatives

| Crate               | Algorithms | Constraints     | 1D solvers                  | Build deps       |
| :------------------ | :--------- | :-------------- | :-------------------------- | :--------------- |
| argmin              | ~20        | External only   | BrentOpt, BrentRoot, Golden | None (pure Rust) |
| nlopt               | ~43        | Full (eq+ineq)  | None                        | cmake + C        |
| cobyla              | 1          | Inequality      | No                          | libc only        |
| optimization-engine | 3          | Projection sets | No                          | None             |

argmin is the correct choice for this project: pure Rust, good ndarray integration, the `Solver` trait enables custom algorithms, and maintenance is active.
