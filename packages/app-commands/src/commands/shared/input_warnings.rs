use itertools::Itertools;
use treetime::progress::{LogSink, RunWarning, RunWarningKind};

const NAMES_SHOWN: usize = 3;

pub fn warn_duplicate_names(
  log: &dyn LogSink,
  kind: RunWarningKind,
  subject: &str,
  consequence: &str,
  names: &[String],
) {
  if names.is_empty() {
    return;
  }
  let names = names.iter().sorted().dedup().cloned().collect_vec();
  log.warning(&RunWarning {
    kind,
    message: format!("{subject} {}. {consequence}", name_list(&names)),
    names,
  });
}

pub fn name_list(names: &[String]) -> String {
  let shown = names.iter().take(NAMES_SHOWN).join(", ");
  if names.len() > NAMES_SHOWN {
    format!("{shown}, and {} more", names.len() - NAMES_SHOWN)
  } else {
    shown
  }
}
