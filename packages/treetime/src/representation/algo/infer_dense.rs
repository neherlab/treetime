pub fn infer_dense() -> bool {
  // TODO: include a heuristic when to use the dense or sparse representation.
  // Generally, long branches in the tree --> dense, short branches --> sparse.
  // For small datasets, it doesn't really matter, but for large ones the sparse is more memory efficient when
  // branches are short
  false
}
