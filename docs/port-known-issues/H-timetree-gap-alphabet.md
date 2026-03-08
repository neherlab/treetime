# Gap character not handled in alphabet

The alphabet used in marginal dense timetree inference does not handle the gap character `-`. Alignments containing gaps trigger an error during sequence processing. The ebola_20 dataset reproduces this, and the [marginal dense runner test](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L32) has ebola_20 commented out. Marginal sparse and poisson paths are unaffected or share the dataset exclusion with the [zero-length branch panic](H-timetree-zero-length-panic.md).
