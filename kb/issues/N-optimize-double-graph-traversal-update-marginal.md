# Optimize marginal updates traverse mixed representations separately

Optimize stores dense and sparse marginal partitions in different concrete vectors. Calls to `update_marginal()` therefore traverse the graph once for each representation at several pipeline boundaries [packages/treetime/src/optimize/pipeline.rs#L104-L107](../../packages/treetime/src/optimize/pipeline.rs#L104-L107) [packages/treetime/src/optimize/pipeline.rs#L150-L152](../../packages/treetime/src/optimize/pipeline.rs#L150-L152). Numerical-failure and worsened-likelihood recovery repeat that pair [packages/treetime/src/optimize/run_loop.rs#L108-L113](../../packages/treetime/src/optimize/run_loop.rs#L108-L113) [packages/treetime/src/optimize/run_loop.rs#L133-L137](../../packages/treetime/src/optimize/run_loop.rs#L133-L137).

Each call performs a backward and forward graph pass. Source inspection establishes the repeated traversal, but does not establish that traversal and synchronization overhead materially affect complete-command runtime relative to per-partition computation.

## Evidence required

- Profile mixed-representation runs and separate graph-traversal overhead from marginal kernel work.
- Preserve dense and sparse numerical results and failure atomicity when evaluating consolidation options.
- Create an implementation ticket only after measurements identify a meaningful bottleneck and the chosen representation preserves static dispatch or supplies an equally explicit type boundary.
