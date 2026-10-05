#[cfg(test)]
mod tests {
  use helpers::exports;
  use pretty_assertions::assert_eq;
  use treetime_utils::vec_of_owned;

  const INDEX_D_TS: &str = include_str!("../../index.d.ts");

  #[test]
  fn test_exports_of_the_addon_are_the_allowlist() {
    assert_eq!(
      vec_of_owned![
        "class Backend: constructor(), fetch(request, scope, onReply), rejectRequest(seq, message)",
        "class PortExchange: abort()",
        "function appStartup()",
        "interface AppStartup",
        "interface PortHeader",
        "type PortMessage",
        "type PortReply",
        "interface PortRequest",
      ],
      exports(INDEX_D_TS)
    );
  }

  #[test]
  fn test_exports_parse_classes_functions_and_types() {
    let source = "export declare class A {\n  constructor()\n  run(x: number, f: ((a: B) => void)): C\n}\n\n\
                  export declare function make(options?: { a: number }): A\n\nexport type T =\n  | { kind: 'a' }\n";
    assert_eq!(
      vec_of_owned!["class A: constructor(), run(x, f)", "function make(options)", "type T"],
      exports(source)
    );
  }

  mod helpers {
    pub(super) fn exports(source: &str) -> Vec<String> {
      let mut found = vec![];
      let mut lines = source.lines();
      while let Some(line) = lines.next() {
        if let Some(rest) = line.strip_prefix("export declare class ") {
          let name = rest.split_whitespace().next().unwrap_or_default();
          let members = lines
            .by_ref()
            .take_while(|member| *member != "}")
            .map(|member| signature(member.trim()))
            .collect::<Vec<_>>()
            .join(", ");
          found.push(format!("class {name}: {members}"));
        } else if let Some(rest) = line.strip_prefix("export declare function ") {
          found.push(format!("function {}", signature(rest)));
        } else if let Some(rest) = line.strip_prefix("export interface ") {
          found.push(format!("interface {}", first_word(rest)));
        } else if let Some(rest) = line.strip_prefix("export type ") {
          found.push(format!("type {}", first_word(rest)));
        }
      }
      found
    }

    fn signature(declaration: &str) -> String {
      let Some((name, parameters)) = declaration.split_once('(') else {
        return declaration.to_owned();
      };
      let parameters = top_level_items(parameters)
        .iter()
        .map(|parameter| {
          parameter
            .split(':')
            .next()
            .unwrap_or_default()
            .trim()
            .trim_end_matches('?')
            .to_owned()
        })
        .filter(|parameter| !parameter.is_empty())
        .collect::<Vec<_>>()
        .join(", ");
      format!("{name}({parameters})")
    }

    fn top_level_items(parameters: &str) -> Vec<String> {
      let mut items = vec![String::new()];
      let mut depth = 0_usize;
      let mut previous = ' ';
      for c in parameters.chars() {
        let arrow = previous == '=' && c == '>';
        previous = c;
        match c {
          '(' | '{' | '[' | '<' => depth += 1,
          ')' if depth == 0 => break,
          '>' if arrow => {},
          ')' | '}' | ']' | '>' => depth = depth.saturating_sub(1),
          ',' if depth == 0 => {
            items.push(String::new());
            continue;
          },
          _ => {},
        }
        if let Some(item) = items.last_mut() {
          item.push(c);
        }
      }
      items
    }

    fn first_word(text: &str) -> &str {
      text
        .split(|c: char| !c.is_alphanumeric() && c != '_')
        .next()
        .unwrap_or_default()
    }
  }
}
