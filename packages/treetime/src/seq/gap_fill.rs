use serde::{Deserialize, Serialize};
use treetime_primitives::{AsciiChar, Seq};

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]
pub enum GapFill {
  #[default]
  OnlyTerminal,
  All,
  None,
}

pub fn apply_gap_fill(seq: &mut Seq, mode: GapFill, gap: AsciiChar, unknown: AsciiChar) {
  match mode {
    GapFill::None => {},
    GapFill::All => {
      for ch in seq.as_mut_slice() {
        if *ch == gap {
          *ch = unknown;
        }
      }
    },
    GapFill::OnlyTerminal => {
      let slice = seq.as_mut_slice();
      if slice.is_empty() {
        return;
      }

      let first_non_gap = slice.iter().position(|&ch| ch != gap);
      let Some(first) = first_non_gap else {
        // All-gap sequence: fill everything with unknown
        for ch in slice.iter_mut() {
          *ch = unknown;
        }
        return;
      };

      let last = slice.iter().rposition(|&ch| ch != gap).unwrap_or(first);

      for ch in &mut slice[..first] {
        *ch = unknown;
      }
      for ch in &mut slice[last + 1..] {
        *ch = unknown;
      }
    },
  }
}
