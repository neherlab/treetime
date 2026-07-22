# PhyloXML boolean properties are silently coerced

The PhyloXML reader parses recognized boolean properties (`treetime:bad_branch`, `treetime:date_inferred`) through equality with the string `"true"`. The valid XML Schema value `"1"` becomes false, while every invalid lexical value also becomes false without an error.

> [!NOTE]
> The tree-output refactor removed the former `tree_ir` PhyloXML reader that assigned these booleans via `prop.value == "true"`. The current reader lives in [`packages/treetime-io/src/phyloxml.rs`](../../packages/treetime-io/src/phyloxml.rs). Whether recognized boolean properties are now parsed per the XSD lexical space (accepting `1`/`0`, rejecting other text) and whether `bad_branch`/`date_inferred` still round-trip is **not yet confirmed**; the required contract below stands regardless.

## Evidence

The former property loop assigned both `treetime:bad_branch` and `treetime:date_inferred` using `prop.value == "true"`, so `"1"` and every invalid value silently became false.

XML Schema `boolean` permits exactly `true`, `false`, `1`, and `0`; its value space contains only true and false. [XML Schema Part 2: boolean](https://www.w3.org/TR/xmlschema-2/#boolean)

## Required contract

- Parse `true` and `1` as true.
- Parse `false` and `0` as false.
- Reject every other lexical value instead of coercing it to false.
- Include the PhyloXML property `ref` and original value in the error.

## Potential solutions

- O1. Parse the four XML Schema lexical forms through one shared fallible helper at the property boundary.
- O2. Deserialize recognized property values into a schema-aware boolean type before clade conversion. This centralizes XML lexical validation but requires the property parser to retain contextual `ref` information for errors.

## Recommendation

Use O1: one fallible XML Schema boolean parser for both recognized boolean properties. This matches the external format contract and prevents invalid text from becoming valid TreeIR state.

## Validation

- Parameterized cases for all four legal lexical forms on `REF_BAD_BRANCH` and `REF_DATE_INFERRED`.
- Representative invalid values, including case variants, surrounding whitespace, an empty value, and unrelated text.
- Whole-document tests proving legal values round-trip and invalid values return contextual errors.

