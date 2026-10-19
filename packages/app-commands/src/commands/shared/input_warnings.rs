use itertools::Itertools;
use treetime::progress::{LogSink, RunWarning, RunWarningKind};
use treetime_utils::fmt::string::name_list;

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
