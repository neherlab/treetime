# Clock filter deadlock blocks all default timetree runs

`apply_outlier_bad_branches()` self-deadlocks on a `parking_lot::RwLock`.
`iter_depth_first_postorder_forward` constructs `GraphNodeBackward` which
write-locks the current node's payload via `node.payload().write_arc()`. The
callback re-acquires a write lock on the same payload via
`graph.get_node(node.key).read_arc().payload().write_arc()`. `parking_lot`
write locks are non-reentrant.

All timetree runs with the default `--clock-filter=3` hang after the "Outlier
filtering" log message. All threads block on futex.

## Root cause

[`clock_filter.rs#L86-L108`](../../packages/treetime/src/commands/timetree/optimization/clock_filter.rs#L86-L108):

```rust
graph.iter_depth_first_postorder_forward(|node| {
    // node.payload is already WRITE-locked by GraphNodeBackward::new
    graph
      .get_node(node.key)       // gets Arc<RwLock<Node<N>>>
      .expect("node must exist")
      .read_arc()               // read-locks Node<N> (OK, re-entrant)
      .payload()                // gets Arc<RwLock<N>>
      .write_arc()              // DEADLOCK: write-locks N, same lock as GraphNodeBackward::payload
      .bad_branch = all_children_bad;
});
```

## Fix direction

Use `node.payload` and `node.children` directly, matching the pattern in
[`date_constraints.rs#L34-L69`](../../packages/treetime/src/commands/timetree/optimization/date_constraints.rs#L34-L69)
and
[`assign_dates.rs#L16-L31`](../../packages/treetime/src/commands/timetree/optimization/assign_dates.rs#L16-L31):

```rust
graph.iter_depth_first_postorder_forward(|mut node| {
    let all_children_bad = node.children.iter()
        .all(|(child_payload, _)| child_payload.read_arc().bad_branch);
    node.payload.bad_branch = all_children_bad;
});
```

## Workaround

`--clock-filter=0` bypasses the outlier filtering code path entirely.

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- timetree \
  --tree=data/flu/h3n2/20/tree.nwk \
  --dates=data/flu/h3n2/20/metadata.tsv \
  --outdir=tmp/repro-deadlock data/flu/h3n2/20/aln.fasta.xz
# Hangs after "Outlier filtering" log. Must kill with Ctrl+C.
```
