# Rerooting code is spread over the clock and reroot modules

## Problem

Choosing and applying a new root takes code from four places in `packages/treetime/src/`:

- [clock/reroot.rs](../../packages/treetime/src/clock/reroot.rs): the driver (`reroot_at_best_root`,
  `reroot_undated`), the choice among the reroot methods (`select_root`) and `RerootParams`;
- [clock/find_best_root/](../../packages/treetime/src/clock/find_best_root/): the least-squares root search (cost
  function, exact per-edge optimum) and the method types `RerootSpec`, `RerootMethod` and `RootObjective`;
- [reroot/](../../packages/treetime/src/reroot/): placement on an edge and the topology change (`apply_root_at_edge`,
  `RerootTopologyParams`);
- [timetree/optimization/reroot.rs](../../packages/treetime/src/timetree/optimization/reroot.rs): the dated wrapper
  that also carries the branch model through the reroot.

`optimize` imports `clock::reroot::reroot_undated` to reroot trees without dates, and so depends on the clock
module only for rerooting. The clock fit (`estimate_clock_model_with_reroot_policy`) calls the driver itself, so the
driver cannot move to the top-level `reroot/` without a dependency cycle between `clock` and `reroot`.

> [!IMPORTANT]
> **Decision required.** Where the rerooting code lives:
>
> 1. `clock/reroot/`: the driver and `find_best_root` move into one clock submodule. This is a pure move; `optimize`
>    still imports the undated reroot from `clock`.
> 2. Top-level `reroot/`: the driver and the root search join placement and topology in `reroot/`, which then
>    depends on `clock` for the regression sums (`ClockState`, `clock_regression`). The clock fit stops calling the
>    driver: its callers reroot first and then fit with the root kept. This changes the clock fit's signature
>    (`keep_root` and `reroot_params` go away) and all its callers.

## Validation

The rerooting code sits under one module and does not form a dependency cycle. The workspace tests pass unchanged
apart from their import paths, and timetree and clock outputs are byte-identical before and after the move.
