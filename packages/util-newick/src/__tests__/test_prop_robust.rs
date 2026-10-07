#[cfg(test)]
mod tests {
  use crate::dialect::NewickDialect;
  use crate::nexus::read::nexus_from_str;
  use crate::read::options::{NewickReadOptions, ReadMode};
  use crate::read::stream::{newick_from_str, newick_trees};
  use generators::{gen_mutated_fixture, gen_text};
  use helpers::runs_without_panic;
  use proptest::prelude::*;

  proptest! {
    #![proptest_config(ProptestConfig { cases: 256, ..ProptestConfig::default() })]

    #[test]
    fn test_prop_robust_random_text_reads_without_panic(text in gen_text()) {
      prop_assert!(runs_without_panic(text), "the reader panicked");
    }

    #[test]
    fn test_prop_robust_mutated_fixture_reads_without_panic(text in gen_mutated_fixture()) {
      prop_assert!(runs_without_panic(text), "the reader panicked");
    }
  }

  mod helpers {
    use super::{NewickDialect, NewickReadOptions, ReadMode, newick_from_str, newick_trees, nexus_from_str};
    use std::thread;

    pub(super) fn runs_without_panic(text: String) -> bool {
      thread::Builder::new()
        .stack_size(2 << 20)
        .spawn(move || {
          for dialect in NewickDialect::pairs() {
            for mode in [ReadMode::Strict, ReadMode::Tolerant] {
              let options = NewickReadOptions {
                dialect,
                mode,
                ..NewickReadOptions::default()
              };
              drop(newick_from_str(&text, &options));
              newick_trees(text.as_bytes(), options.clone()).for_each(drop);
              drop(nexus_from_str(&format!("#NEXUS\nBegin Trees;\n{text}"), &options));
              drop(nexus_from_str(&text, &options));
            }
          }
        })
        .unwrap()
        .join()
        .is_ok()
    }
  }

  mod generators {
    use proptest::collection::vec;
    use proptest::prelude::*;

    const FIXTURES: [&str; 10] = [
      "[&R]((A[&rate=1.5,s=\"x,y\",a={1,{2}}]:[&r=1]1[c],B)[&posterior=0.9]:1,C);",
      "((A[&&NHX:S=human:B=90:C=1.2.3]:1,B):1,C[&&NHX:Ev=1>2>dup]);",
      "[&U]((A:1[&B TK02Brlens 0.1],B[&E ibr 2: 0.1]):1,C);",
      "(A,B,((C,(Y)x#H1)c,(x##H1,D)d)e)f;",
      "[&U][&W 0.5]((A:1:90,(B)#H1:::0.3),(#H1:::0.7,C)):0.1;",
      "('it''s':0.1,'a b'[note [nested]]:2e-3,(C,D)80.5/95:1)root;",
      "#NEXUS\nBegin Taxa;\nDimensions ntax=2;\nTaxLabels A B;\nEnd;\nBegin Trees;\nTranslate 1 A, 2 B;\nTree t [&lnP=1] = [&R] (1:1,2:2);\nEnd;\n",
      "(A[5\" tall],B[say \"hi\"]);",
      "\u{feff}(A:.5,B:+1)C",
      "(A:1,B:2);(C:3,D:4);",
    ];

    pub(super) fn gen_text() -> impl Strategy<Value = String> {
      "[(),;:\\[\\]{}'\"&#=>!*A-Za-z0-9 ._\\-\n\t\u{feff}\u{e9}]{0,4096}"
    }

    pub(super) fn gen_mutated_fixture() -> impl Strategy<Value = String> {
      let edit = (
        any::<usize>(),
        0_u8..3,
        prop_oneof![
          Just('('),
          Just(')'),
          Just(','),
          Just(';'),
          Just(':'),
          Just('['),
          Just(']'),
          Just('{'),
          Just('}'),
          Just('\''),
          Just('"'),
          Just('&'),
          Just('#'),
          Just('='),
          Just('1'),
        ],
      );
      (0..FIXTURES.len(), vec(edit, 1..8)).prop_map(|(fixture, edits)| {
        let mut chars: Vec<char> = FIXTURES[fixture].chars().collect();
        for (position, operation, character) in edits {
          let at = position % (chars.len() + 1);
          match operation {
            0 => chars.insert(at, character),
            1 if at < chars.len() => {
              chars.remove(at);
            },
            _ if at < chars.len() => chars[at] = character,
            _ => chars.push(character),
          }
        }
        chars.into_iter().collect()
      })
    }
  }
}
