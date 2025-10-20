use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;

use crate::make_internal_report;

/// Extract inner value from Arc<RwLock<T>>
pub fn unwrap_arc_rwlock<T>(arc_rwlock: Arc<RwLock<T>>) -> Result<T, Report> {
  let inner = Arc::try_unwrap(arc_rwlock)
    .map_err(|val| make_internal_report!("Failed to unwrap Arc<T>: it still has references"))
    .map(|rwlock| rwlock.into_inner())?;
  Ok(inner)
}
