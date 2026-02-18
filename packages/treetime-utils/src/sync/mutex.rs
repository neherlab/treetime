use eyre::Report;
use parking_lot::{Mutex, RwLock};
use std::sync::Arc;

use crate::make_internal_report;

/// Extract inner value from Arc<RwLock<T>>
pub fn unwrap_arc_rwlock<T>(arc_rwlock: Arc<RwLock<T>>) -> Result<T, Report> {
  let inner = Arc::try_unwrap(arc_rwlock)
    .map_err(|val| make_internal_report!("Failed to unwrap Arc<T>: it still has references"))
    .map(|rwlock| rwlock.into_inner())?;
  Ok(inner)
}

/// Extract error from a parallel error collector.
///
/// Used with parallel traversals that store the first error in `Arc<Mutex<Option<Report>>>`.
/// Handles both cases: unique ownership (try_unwrap succeeds) and shared references.
pub fn extract_parallel_error(error: Arc<Mutex<Option<Report>>>) -> Result<(), Report> {
  match Arc::try_unwrap(error) {
    Ok(mutex) => {
      if let Some(e) = mutex.into_inner() {
        return Err(e);
      }
    },
    Err(arc) => {
      if let Some(e) = arc.lock().take() {
        return Err(e);
      }
    },
  }
  Ok(())
}
