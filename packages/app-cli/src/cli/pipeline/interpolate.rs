use eyre::Report;
use itertools::Itertools;
use minijinja::{Environment, UndefinedBehavior};
use serde::Deserialize;
use serde_json::{Map, Value};
use std::collections::{BTreeSet, HashSet};
use treetime_utils::{make_error, make_report};

/// Namespaces available to a template, as they appear at the root of a `{{ ... }}` expression.
pub const NAMESPACES: [&str; 3] = ["vars", "env", "steps"];

/// Renders `{{ ... }}` templates embedded in config string leaves.
///
/// Minijinja is used in substitution-only mode: `{% %}` statement blocks are never documented or
/// used, and undefined names are a hard error (`UndefinedBehavior::Strict`) so a typo surfaces
/// instead of silently rendering empty. Two rendering modes preserve types (the C-tension): a leaf
/// that is exactly one expression (`"{{ vars.rate }}"`) is evaluated to a typed value (number, bool,
/// or string), while an expression embedded in surrounding text (`"{{ vars.dir }}/tree.nwk"`) renders
/// to a string.
pub struct Interpolator {
  env: Environment<'static>,
}

impl Default for Interpolator {
  fn default() -> Self {
    let mut env = Environment::new();
    env.set_undefined_behavior(UndefinedBehavior::Strict);
    Self { env }
  }
}

impl Interpolator {
  /// Render one string leaf against `context`.
  ///
  /// A leaf with no `{{` is returned verbatim. A whole-value expression is evaluated to a typed
  /// value; any other template renders to a string.
  pub fn interpolate_str(&self, leaf: &str, context: &Value) -> Result<Value, Report> {
    if !leaf.contains("{{") {
      return Ok(Value::String(leaf.to_owned()));
    }

    if let Some(expr) = whole_value_expr(leaf) {
      let compiled = self
        .env
        .compile_expression(expr)
        .map_err(|err| minijinja_error("expression", leaf, &err))?;
      let value = compiled
        .eval(context)
        .map_err(|err| minijinja_error("expression", leaf, &err))?;
      let value = Value::deserialize(value).map_err(|err| make_report!("evaluating `{leaf}`: {err}"))?;
      Ok(value)
    } else {
      let rendered = self
        .env
        .render_str(leaf, context)
        .map_err(|err| minijinja_error("template", leaf, &err))?;
      Ok(Value::String(rendered))
    }
  }

  /// Render every string leaf inside `value` against `context`, preserving structure and types.
  pub fn interpolate_value(&self, value: &Value, context: &Value) -> Result<Value, Report> {
    map_string_leaves(value, &mut |leaf| self.interpolate_str(leaf, context))
  }

  /// Dotted references a leaf uses, e.g. `{"vars.data", "steps.tt.outputs.nwk"}`.
  ///
  /// Empty when the leaf has no template. `nested = true` yields the full dotted path so the caller
  /// can validate each reference against the namespace it names.
  pub fn references(&self, leaf: &str) -> Result<HashSet<String>, Report> {
    if !leaf.contains("{{") {
      return Ok(HashSet::new());
    }
    let template = self
      .env
      .template_from_str(leaf)
      .map_err(|err| minijinja_error("template", leaf, &err))?;
    Ok(template.undeclared_variables(true))
  }
}

/// Resolve the `vars` block in dependency order, then return the concrete values.
///
/// A var value may reference earlier-resolved vars and `env` (`workdir: "{{ vars.base }}/run"`),
/// resolved in repeated passes until a fixed point. A var may not reference `steps` (steps are
/// produced later). A set of vars that never becomes resolvable is a cycle (or an unknown reference)
/// and is reported rather than looped forever.
pub fn resolve_vars(
  interp: &Interpolator,
  raw: &Map<String, Value>,
  env: &Value,
) -> Result<Map<String, Value>, Report> {
  let mut resolved: Map<String, Value> = Map::new();
  let mut pending: Vec<(String, Value)> = raw.iter().map(|(name, value)| (name.clone(), value.clone())).collect();

  while !pending.is_empty() {
    let mut progressed = false;
    let mut next = Vec::new();
    for (name, value) in pending {
      let deps = var_dependencies(interp, &name, &value)?;
      if deps.iter().all(|dep| resolved.contains_key(dep)) {
        let context = template_context(&resolved, env);
        let rendered = interp
          .interpolate_value(&value, &context)
          .map_err(|err| eyre::eyre!("resolving var `{name}`: {err}"))?;
        resolved.insert(name, rendered);
        progressed = true;
      } else {
        next.push((name, value));
      }
    }
    if !progressed {
      let unresolved = next.iter().map(|(name, _)| format!("`{name}`")).sorted().join(", ");
      return make_error!("cannot resolve pipeline `vars` ({unresolved}): circular or unknown reference among them");
    }
    pending = next;
  }

  Ok(resolved)
}

/// Wrap the `vars` mapping and the process environment as a `{ vars, env }` template context.
pub(crate) fn template_context(vars: &Map<String, Value>, env: &Value) -> Value {
  let mut context = Map::new();
  context.insert("vars".to_owned(), Value::Object(vars.clone()));
  context.insert("env".to_owned(), env.clone());
  Value::Object(context)
}

/// Rebuild `value`, replacing every string leaf with the result of `f`, preserving structure.
///
/// Arrays and objects recurse; non-string scalars are cloned untouched. The one traversal is shared
/// by typed interpolation and by the `{{ steps.* }}` substitution pass, which differ only in `f`.
pub(crate) fn map_string_leaves(
  value: &Value,
  f: &mut impl FnMut(&str) -> Result<Value, Report>,
) -> Result<Value, Report> {
  match value {
    Value::String(leaf) => f(leaf),
    Value::Array(items) => items
      .iter()
      .map(|item| map_string_leaves(item, f))
      .collect::<Result<Vec<_>, _>>()
      .map(Value::Array),
    Value::Object(map) => map
      .iter()
      .map(|(key, val)| Ok((key.clone(), map_string_leaves(val, f)?)))
      .collect::<Result<Map<_, _>, Report>>()
      .map(Value::Object),
    other => Ok(other.clone()),
  }
}

/// Names of other vars that `value` depends on, rejecting a `steps` reference from within `vars`.
fn var_dependencies(interp: &Interpolator, name: &str, value: &Value) -> Result<BTreeSet<String>, Report> {
  let mut deps = BTreeSet::new();
  collect_var_dependencies(interp, name, value, &mut deps)?;
  Ok(deps)
}

fn collect_var_dependencies(
  interp: &Interpolator,
  name: &str,
  value: &Value,
  deps: &mut BTreeSet<String>,
) -> Result<(), Report> {
  match value {
    Value::String(leaf) => {
      for reference in interp.references(leaf)? {
        let mut segments = reference.split('.');
        match segments.next() {
          Some("vars") => {
            if let Some(dep) = segments.next() {
              deps.insert(dep.to_owned());
            }
          },
          Some("steps") => {
            return make_error!("pipeline var `{name}` references `steps`; vars may only use `vars` and `env`");
          },
          _ => {},
        }
      }
    },
    Value::Array(items) => {
      for item in items {
        collect_var_dependencies(interp, name, item, deps)?;
      }
    },
    Value::Object(map) => {
      for val in map.values() {
        collect_var_dependencies(interp, name, val, deps)?;
      }
    },
    _ => {},
  }
  Ok(())
}

/// Convert a minijinja error into an eyre report naming the offending leaf.
fn minijinja_error(kind: &str, leaf: &str, err: &minijinja::Error) -> Report {
  match err.range() {
    Some(range) => {
      let snippet = leaf.get(range).unwrap_or(leaf);
      make_report!("invalid {kind} in `{leaf}`: {err} (near `{snippet}`)")
    },
    None => make_report!("invalid {kind} in `{leaf}`: {err}"),
  }
}

/// Extract the inner expression of a whole-value leaf `"{{ EXPR }}"`, or `None` when the leaf mixes
/// a template with surrounding text or contains more than one expression.
fn whole_value_expr(leaf: &str) -> Option<&str> {
  let trimmed = leaf.trim();
  let inner = trimmed.strip_prefix("{{")?.strip_suffix("}}")?;
  if inner.contains("{{") || inner.contains("}}") {
    return None;
  }
  Some(inner.trim())
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use treetime_utils::assert_error;

  fn empty_env() -> Value {
    json!({})
  }

  // A whole-value numeric reference evaluates to a JSON number, not the string "0.003", so a typed
  // command field like `clock_rate` deserializes correctly.
  #[test]
  fn test_interpolate_whole_value_ref_preserves_number() {
    let interp = Interpolator::default();
    let context = json!({ "vars": { "rate": 0.003 }, "env": {} });
    let value = interp.interpolate_str("{{ vars.rate }}", &context).unwrap();
    assert_eq!(json!(0.003), value);
  }

  // A whole-value boolean reference stays a JSON boolean.
  #[test]
  fn test_interpolate_whole_value_ref_preserves_bool() {
    let interp = Interpolator::default();
    let context = json!({ "vars": { "flag": true }, "env": {} });
    let value = interp.interpolate_str("{{ vars.flag }}", &context).unwrap();
    assert_eq!(json!(true), value);
  }

  // An embedded reference renders to a string with the surrounding text intact.
  #[test]
  fn test_interpolate_embedded_ref_renders_string() {
    let interp = Interpolator::default();
    let context = json!({ "vars": { "data": "data/flu" }, "env": {} });
    let value = interp.interpolate_str("{{ vars.data }}/tree.nwk", &context).unwrap();
    assert_eq!(json!("data/flu/tree.nwk"), value);
  }

  // A leaf with no template is returned verbatim.
  #[test]
  fn test_interpolate_plain_leaf_unchanged() {
    let interp = Interpolator::default();
    let value = interp.interpolate_str("data/flu/tree.nwk", &empty_env()).unwrap();
    assert_eq!(json!("data/flu/tree.nwk"), value);
  }

  // Under strict undefined behavior, referencing a name absent from the context is an error. This is
  // how backward-only step references are enforced: a step not yet in `steps` is undefined.
  #[test]
  fn test_interpolate_undefined_reference_errors() {
    let interp = Interpolator::default();
    let context = json!({ "vars": {}, "env": {}, "steps": {} });
    let result = interp.interpolate_str("{{ steps.later.outputs.nwk }}", &context);
    assert!(result.is_err(), "undefined step reference must error under strict mode");
  }

  // Vars resolve in dependency order: `workdir` uses `base` regardless of declaration order.
  #[test]
  fn test_resolve_vars_dependency_order() {
    let interp = Interpolator::default();
    let raw = json!({ "workdir": "{{ vars.base }}/run", "base": "tmp" })
      .as_object()
      .unwrap()
      .clone();
    let resolved = resolve_vars(&interp, &raw, &empty_env()).unwrap();
    assert_eq!(json!("tmp/run"), resolved["workdir"]);
    assert_eq!(json!("tmp"), resolved["base"]);
  }

  // A var may reference `env`.
  #[test]
  fn test_resolve_vars_uses_env() {
    let interp = Interpolator::default();
    let raw = json!({ "home": "{{ env.HOME }}/data" }).as_object().unwrap().clone();
    let env = json!({ "HOME": "/users/me" });
    let resolved = resolve_vars(&interp, &raw, &env).unwrap();
    assert_eq!(json!("/users/me/data"), resolved["home"]);
  }

  // A cycle among vars is reported rather than looped forever.
  #[test]
  fn test_resolve_vars_cycle_errors() {
    let interp = Interpolator::default();
    let raw = json!({ "a": "{{ vars.b }}", "b": "{{ vars.a }}" })
      .as_object()
      .unwrap()
      .clone();
    let result = resolve_vars(&interp, &raw, &empty_env());
    assert_error!(
      result,
      "cannot resolve pipeline `vars` (`a`, `b`): circular or unknown reference among them"
    );
  }

  // A var that references `steps` is rejected: steps are produced later than vars.
  #[test]
  fn test_resolve_vars_rejects_steps_reference() {
    let interp = Interpolator::default();
    let raw = json!({ "x": "{{ steps.tt.outputs.nwk }}" })
      .as_object()
      .unwrap()
      .clone();
    let result = resolve_vars(&interp, &raw, &empty_env());
    assert_error!(
      result,
      "pipeline var `x` references `steps`; vars may only use `vars` and `env`"
    );
  }
}
