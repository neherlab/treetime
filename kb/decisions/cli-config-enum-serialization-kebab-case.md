# Kebab-case serialization for all serialized enums

Every enum in the workspace that derives deser serializes and deserializes its variants in `kebab-case`, via an unconditional `#[deser(rename_all = "kebab-case")]` on each enum; enums that also derive `JsonSchema` repeat it as `#[schemars(rename_all = "kebab-case")]` for the schema. This matches the values clap's `ValueEnum` derive already accepts on the command line, whose default rename is also `kebab-case`.

## Context

The `--config <file>` argument accepts a command's full argument object as JSON or YAML and merges it under the precedence explicit CLI flag > config file > default. Deserialization of that object goes through deser, while the same fields parse through clap on the command line. Before this decision the two disagreed: the serializer emitted and expected PascalCase variant names (`Marginal`, `LeastSquares`) while clap accepted kebab-case (`marginal`, `least-squares`). A config file therefore needed different casing than the equivalent CLI flag. Using `kebab-case` everywhere gives one representation for both paths.

The redundant per-enum `#[cfg_attr(feature = "clap", value(rename_all = "kebab-case"))]` attributes are removed, because `kebab-case` is clap's `ValueEnum` default.

## Acronym variants

The `kebab-case` of deser inserts a separator before every uppercase letter (`apply_to_variant()` in `deser-derive` 0.9.0), so an acronym run splits: `JC69` becomes `j-c69`. clap derives command-line values through the `heck` crate, which keeps an acronym run together: `JC69` becomes `jc69`. The two conventions agree for ordinary PascalCase (`LeastSquares` -> `least-squares` on both sides) and diverge only on consecutive capitals.

`GtrModelName` [packages/treetime/src/gtr/get_gtr.rs](../../packages/treetime/src/gtr/get_gtr.rs) is the only command-line enum with acronym variants (`JC69`, `HKY85`, `TN93`). Each carries an explicit `#[deser(rename = "jc69")]` (and the same `#[schemars(...)]`), `"hky85"`, `"tn93"` so the configuration emits and accepts exactly the spelling clap uses, preserving CLI and config parity. Enums without acronym variants rely on `rename_all` alone.

## Impact

Serialized enum strings in output JSON are recased. The GTR model output `model_name` changes from `"JC69"` to `"jc69"` and `model_type` stays lowercase; internal graph and partition JSON tags such as `Estimated`, `Point`, and `Fitch` become `estimated`, `point`, and `fitch`. `GtrModelName` and `AlphabetName` reach augur-facing output, so downstream consumers see the lowercase spelling. The change is representation only; no numeric value moves.
