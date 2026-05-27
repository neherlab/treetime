# Augur node data JSON format

Implementation reference for producing augur-compatible node data JSON from treetime v1. The node data JSON is the intermediate format produced by augur commands (`augur ancestral`, `augur refine`, `augur traits`). Multiple node data files are merged by `augur export v2` to build the final auspice dataset.

## Format overview

A flat JSON object with metadata at the top level and a `nodes` dict mapping node names to per-node property objects. Written with `indent=2`, `sort_keys=true` (Python `json.dump` defaults in augur). No formal JSON schema for the full format. Only `annotations` has a sub-schema ([`schema-annotations.json`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/data/schema-annotations.json), draft-06). Validation is programmatic: [`node_data_file.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/util_support/node_data_file.py#L62-L114) checks `nodes` is a dict and validates `annotations`.

Multiple files are deep-merged by [`node_data_reader.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/util_support/node_data_reader.py): dict values merge recursively, non-dict values are overwritten by the later file, mixing dict/non-dict for the same key raises an error. `generated_by` is filtered out during merging.

Nextstrain workflows run a sequence of augur commands on the same tree, each adding different per-node fields to a separate node data JSON file. A typical pipeline:

1. `augur refine` -> `branch_lengths.json`: dates (`numdate`), branch lengths, clock model
2. `augur ancestral` -> `nt_muts.json`: nucleotide mutations (`muts`), sequences, optionally AA mutations
3. `augur traits` -> `traits.json`: discrete trait assignments (`region`, `country`), confidence, entropy
4. `augur clades` -> `clades.json`: clade membership based on mutation signatures
5. Custom scripts -> `epiweeks.json`, `recency.json`, etc.: derived metadata fields

These are different analyses, not the same command over different data. Each produces non-overlapping per-node fields. `augur export v2 --node-data branch_lengths.json nt_muts.json traits.json clades.json ...` deep-merges all files into one combined dict (dict values merge recursively, non-dict values overwrite), then builds the auspice JSON from the merged result. The format is an open extension point: any `{"nodes": {name: {key: value}}}` file can participate in the merge.

Neither augur nor any workflow calls treetime CLI to produce node data JSON. Augur calls treetime v0 as a Python library, reads BioPython tree node attributes, and constructs the JSON itself. Treetime v0 has no awareness of this format.

## `generated_by` (all commands)

Added by [`write_augur_json()`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/utils.py#L123-L137):

```python
data["generated_by"] = {"program": "augur", "version": get_augur_version()}
```

For treetime v1: `{"program": "treetime", "version": "<version>"}`. Filtered out during merging by augur. Not consumed downstream.

## `augur ancestral` node data JSON

Corresponds to treetime v1 `ancestral` command.

### Transformation layer

Augur calls `TreeAnc.infer_ancestral_sequences()`, then builds JSON from tree node attributes. The transformation functions are in [`ancestral.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/ancestral.py).

### Complete structure

```json
{
  "generated_by": { "program": "augur", "version": "25.0.0" },
  "annotations": {
    "nuc": { "start": 1, "end": 50, "strand": "+", "type": "source" }
  },
  "reference": { "nuc": "AAAACAAAAATGCTCTGTGGGTAAAAAAAAAAAAACTACTTGACCATAAA" },
  "mask": "00000000000000000000000000000000000000000000000000",
  "nodes": {
    "node_root": {
      "muts": ["A5C", "C14T"],
      "sequence": "AAAACAAAAATGCTCTGTGGGTAAAAAAAAAAAAACTACTTGACCATAAA"
    },
    "sample_A": {
      "muts": ["G33C", "C39T"],
      "sequence": "AAAACAGAAATGCTCTGCGGGTAAAAAAAAAACAACTATTTGTCCATAAA"
    }
  }
}
```

With AA reconstruction (`--genes`), per-node objects gain `aa_muts` and the root gains `aa_sequences`. Top-level `annotations` and `reference` gain per-gene entries:

```json
{
  "annotations": {
    "nuc": { "start": 1, "end": 10794, "strand": "+", "type": "source" },
    "ENV": { "start": 937, "end": 2413, "strand": "+", "type": "CDS" }
  },
  "reference": { "nuc": "ATGC...", "ENV": "MDFL..." },
  "nodes": {
    "root": {
      "muts": [],
      "sequence": "ATGC...",
      "aa_muts": { "ENV": [] },
      "aa_sequences": { "ENV": "MDFL..." }
    },
    "sample_A": {
      "muts": ["A123T"],
      "sequence": "ATGC...",
      "aa_muts": { "ENV": ["D614G", "N501Y"] }
    }
  }
}
```

### Field specifications

#### `annotations` (top-level)

Object mapping sequence/CDS names to coordinate objects. Validated against [`schema-annotations.json`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/data/schema-annotations.json).

`annotations.nuc`: always present. `start` = 1, `end` = alignment length, `strand` = `"+"`, `type` = `"source"`. Built from `len(root_seq)`.

`annotations.<gene>`: one per gene when AA reconstruction is done. Built by [`genome_features_to_auspice_annotation()`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/utils.py#L582-L632) from BioPython `SeqFeature` objects. SimpleLocation: `start` = `int(feat.location.start) + 1` (BioPython zero-origin to 1-based), `end` = `int(feat.location.end)` (already 1-based inclusive in GFF convention). CompoundLocation: `segments` array of `{start, end}` instead of top-level start/end. `strand` = `{+1: "+", -1: "-"}[feat.location.strand]`. `type` = `feat.type` (e.g. `"CDS"`). Optional `seqid` if `ref_seq_name` provided.

Consumed by `augur export v2` verbatim as `meta.genome_annotations`.

#### `reference` (top-level)

Object mapping `"nuc"` and gene names to sequence strings.

`reference.nuc`: if `--root-sequence` provided, that file's sequence. Otherwise the TreeTime-inferred root sequence. Augur code (`ancestral.py:298-301`):

```python
if root_sequence:
    root_seq = str(root_sequence.seq)
else:
    root_seq = tt.sequence(T.root, as_string=True, reconstructed=True)
```

`reference.<gene>`: AA sequence for each gene. From `--aa-root-sequence`, translated from nuc ref, or TreeTime-inferred AA root.

Consumed by `augur export v2` as `root_sequence` in auspice JSON (inline or sidecar file).

#### `mask` (top-level, optional)

String of `"0"` and `"1"` characters, length = alignment length. `"1"` = masked (every tip is ambiguous at this position). Built by [`create_mask()`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/ancestral.py#L139-L183):

```python
def create_mask(tt, is_vcf):
    if is_vcf:
        mask = np.zeros(tt.data.full_length, dtype=bool)
        # VCF: mask non-variable positions
    else:
        mask = np.ones(tt.data.full_length, dtype=bool)
        for leaf in tt.tree.get_terminals():
            for pos, state in enumerate(tt.sequence(leaf, reconstructed=False)):
                if state not in [tt.gtr.ambiguous, '-']:
                    mask[pos] = False
    return mask
```

The boolean numpy array is converted to a string in `_to_ancestral_json()`:

```python
j['mask'] = "".join(['1' if x else '0' for x in anc["mutations"]["mask"]])
```

`True` (masked) becomes `"1"`, `False` (not masked) becomes `"0"`.

Not consumed by `augur export v2`. Used by augur's VCF reconstruction.

#### `nodes.<name>.muts`

Array of nucleotide mutation strings. Format: `<from><1-based-pos><to>` (e.g. `"A7G"`).

Built by [`collect_mutations()`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/ancestral.py#L185-L229). TreeTime stores mutations as `(ancestral_char, 0-based-position, derived_char)` tuples on `node.mutations`. Augur converts:

```python
# all nodes: skip masked positions (mask[pos]=True means masked)
data[n.name] = [a + str(int(pos) + 1) + d
                for a, pos, d in n.mutations if not mask[int(pos)]]

# root WITH --root-sequence: overwrite with char-by-char diff against reference
if reference_sequence:
    data[tt.tree.root.name] = []
    for pos, (root_state, tree_state) in enumerate(
            zip(reference_sequence, tt.sequence(tt.tree.root, reconstructed=infer_ambiguous, as_string=True))):
        if mask[pos]:
            continue
        if root_state != tree_state:
            data[tt.tree.root.name].append(f"{root_state}{pos+1}{tree_state}")
```

Root WITHOUT `--root-sequence`: `n.mutations` returns `[]` for root (no parent to compare against).

Non-root: mutations relative to parent.

Consumed by `augur export v2` as `branch_attrs.mutations.nuc`. Empty lists excluded.

#### `nodes.<name>.sequence`

Full reconstructed nucleotide sequence string. Present when `full_sequences=True` (FASTA input, not VCF).

Built by [`collect_sequences()`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/ancestral.py#L231-L269):

```python
for n in tt.tree.find_clades():
    tmp = tt.sequence(n, reconstructed=infer_ambiguous, as_string=False)  # mutable numpy array
    if ref_mask is not None:   # infer_ambiguous=True AND reference provided
        tmp[mask] = ref_mask   # masked positions get reference base
    else:
        tmp[mask] = tt.gtr.ambiguous  # masked positions get 'N'
    sequences[n.name] = "".join(tmp)
```

Not consumed by `augur export v2` (excluded by `node_data_prop_is_normal_trait()`).

#### `nodes.<name>.aa_muts`

Object mapping gene name to array of AA mutation strings. Present when AA reconstruction done.

Format identical to `muts` but AA characters. Each gene has its own list. Empty gene lists included in node data (unlike auspice where they are excluded).

Consumed by `augur export v2` as `branch_attrs.mutations.<gene>`. Empty lists excluded.

#### `nodes.<name>.aa_sequences`

Object mapping gene name to AA sequence string. Present only on the root node when AA reconstruction done.

Not consumed by `augur export v2`.

### Field presence matrix

| Field                  | Nuc (FASTA) | Nuc (VCF) | With AA   |
| ---------------------- | ----------- | --------- | --------- |
| `annotations.nuc`      | yes         | yes       | yes       |
| `annotations.<gene>`   | no          | no        | yes       |
| `reference.nuc`        | yes         | yes       | yes       |
| `reference.<gene>`     | no          | no        | yes       |
| `mask`                 | yes         | yes       | yes       |
| `nodes.*.muts`         | yes         | yes       | yes       |
| `nodes.*.sequence`     | yes         | no        | yes       |
| `nodes.*.aa_muts`      | no          | no        | yes       |
| `nodes.*.aa_sequences` | no          | no        | root only |

## `augur refine` node data JSON

Corresponds to treetime v1 `timetree` command.

### Transformation layer

Augur calls `TreeTime.run()`, then extracts node attributes via `collect_node_data()`:

```python
def collect_node_data(T, attributes):
    data = {}
    for n in T.find_clades():
        data[n.name] = {
            attr: getattr(n, attr)
            for attr in attributes
            if getattr(n, attr, None) is not None
        }
    return data
```

The `attributes` list starts as `['branch_length', 'confidence']` and grows conditionally:

```python
# always:
attributes = ['branch_length', 'confidence']
# timetree adds:
attributes.extend(['numdate', 'clock_length', 'mutation_length', 'raw_date', 'date', 'date_inferred'])
# with --date-confidence:
attributes.append('num_date_confidence')
```

### Complete structure

```json
{
  "generated_by": { "program": "augur", "version": "25.0.0" },
  "alignment": "path/to/alignment.fasta",
  "input_tree": "path/to/tree.nwk",
  "clock": {
    "rate": 0.003,
    "intercept": -2.21,
    "rtt_Tmrca": 2012.3,
    "cov": [
      [1e-8, 0],
      [0, 0.5]
    ],
    "rate_std": 1e-4
  },
  "nodes": {
    "<name>": {
      "branch_length": 0.00119745,
      "confidence": 0.998,
      "numdate": 2013.999,
      "clock_length": 0.001,
      "mutation_length": 0.002,
      "raw_date": "2013-12-XX",
      "date": "2013-12-31",
      "date_inferred": true,
      "num_date_confidence": [2013.5, 2014.5]
    }
  }
}
```

### Per-node field specifications

#### `branch_length`

Float. BioPython node attribute from TreeTime. Substitutions/site by default. With `--divergence-units=mutations`: integer mutation count (augur replaces the float with a count of non-ambiguous, non-masked mutations from `node.mutations`).

#### `confidence`

Float. BioPython tree node bootstrap/support value from the input Newick/Nexus tree. Not computed by augur or TreeTime. Augur just reads it via `getattr(n, 'confidence')`.

Consumed by `augur export v2` as a normal trait (coloring `node_attrs.confidence.value`).

#### `numdate`

Float. TreeTime node attribute. Decimal year (e.g. `2013.999`). Set by TreeTime during time tree inference.

Consumed by `augur export v2` as `node_attrs.num_date.value`.

#### `clock_length`

Float. TreeTime node attribute. Expected branch length from clock model.

Not consumed by `augur export v2` as a trait (excluded by `node_data_prop_is_normal_trait()`). Used for divergence calculation (`mutation_length` preferred over `branch_length`).

#### `mutation_length`

Float. TreeTime node attribute. Observed mutations normalized to substitutions/site. With `--divergence-units=mutations`: integer count.

Not consumed as trait. Used for cumulative `div` computation.

#### `raw_date` and `date`

Strings. Set differently:

`raw_date`: set by augur (not TreeTime) from the metadata TSV, tips only:

```python
for n in T.get_terminals():
    if n.name in metadata.index and METADATA_DATE_COLUMN in metadata.columns:
        n.raw_date = metadata.at[n.name, METADATA_DATE_COLUMN]
```

Internal nodes do not have this attribute (skipped by `getattr(n, attr, None) is not None`).

`date`: set by TreeTime during inference via `datestring_from_numeric(node.numdate)` which converts decimal year to `"YYYY-MM-DD"` format.

Not consumed as traits (excluded by `node_data_prop_is_normal_trait()`). In export: `raw_date` -> `node_attrs.num_date.raw_value`, `date_inferred` -> `node_attrs.num_date.inferred` (tips only).

#### `date_inferred`

Bool. Computed by augur (not TreeTime) after `tt.run()` returns:

```python
for node in T.find_clades():
    inferred = not isinstance(node.raw_date_constraint, float)
    setattr(node, 'date_inferred', inferred)
```

`raw_date_constraint` is set by TreeTime during date parsing: exact numeric date produces a `float`, date range produces a tuple/array, absent date produces `None`. So `date_inferred=False` for known dates, `True` for inferred or absent.

#### `num_date_confidence`

2-element float array. 90% HPD bounds from TreeTime's marginal posterior:

```python
n.num_date_confidence = list(tt.get_max_posterior_region(n, 0.9))
```

Consumed by `augur export v2` as `node_attrs.num_date.confidence`.

### Top-level fields

#### `clock`

Present for timetree mode. Built from TreeTime's `date2dist` regression object:

```python
node_data['clock'] = {
    'rate': tt.date2dist.clock_rate,
    'intercept': tt.date2dist.intercept,
    'rtt_Tmrca': -tt.date2dist.intercept / tt.date2dist.clock_rate
}
if hasattr(tt.date2dist, "cov") and tt.date2dist.cov is not None:
    node_data["clock"]["cov"] = tt.date2dist.cov        # 2x2 numpy array -> list of lists
    node_data["clock"]["rate_std"] = np.sqrt(tt.date2dist.cov[0, 0])
```

Not consumed by `augur export v2`.

#### `skyline` (optional)

4-element array `[x_values, lower_conf, y_values, upper_conf]` from `tt.merger_model.skyline_inferred()`. Present with `--coalescent=skyline`. Not consumed by `augur export v2`.

#### `alignment` and `input_tree`

Paths to input files. Set from CLI args. Can be `null` (when no alignment).

### Divergence units

When `--divergence-units=mutations`, augur replaces the float `branch_length` and `mutation_length` with integer mutation counts. The count excludes ambiguous characters and fully-ambiguous positions:

```python
n_muts = len([
    position
    for ancestral, position, derived in node.mutations
    if are_sequence_states_different(ancestral, derived) and position not in ambiguous_positions
])
```

Where `are_sequence_states_different` returns False if either state is `'-'` or ambiguous (`'N'`/`'X'`), and uses the profile map to check non-overlapping states.

### Downstream consumption

| Node data field       | Auspice JSON target              | Notes                                       |
| --------------------- | -------------------------------- | ------------------------------------------- |
| `numdate`             | `node_attrs.num_date.value`      | `format_number()` applied                   |
| `num_date_confidence` | `node_attrs.num_date.confidence` | Each element through `format_number()`      |
| `date_inferred`       | `node_attrs.num_date.inferred`   | Tips only                                   |
| `raw_date`            | `node_attrs.num_date.raw_value`  | Tips only                                   |
| `branch_length`       | Cumulative `node_attrs.div`      | `mutation_length` preferred if present      |
| `mutation_length`     | Cumulative `node_attrs.div`      | Preferred over `branch_length`              |
| `confidence`          | `node_attrs.confidence.value`    | Normal trait (not excluded)                 |
| `clock_length`        | Excluded                         | `node_data_prop_is_normal_trait()` excludes |
| `date`                | Excluded                         | Same                                        |
| `clock`               | Not consumed                     | Top-level metadata                          |
| `skyline`             | Not consumed                     | Top-level metadata                          |

## `augur traits` node data JSON

Corresponds to treetime v1 `mugration` command.

### Transformation layer

For each `--columns` column, augur calls `mugration_inference()` which wraps `reconstruct_discrete_traits()` (TreeTime's `TreeAnc` mugration). After inference, augur reads attributes from tree nodes.

### Alphabet mapping

TreeTime uses single-letter internal codes. The mapping is:

```python
unique_states = sorted(unique_states)
reverse_alphabet = {state: chr(65+i) for i, state in enumerate(unique_states) if state != missing_data}
letter_to_state = {v: k for k, v in reverse_alphabet.items()}
# missing char gets chr(65 + n_states), one past the real alphabet
missing_char = chr(65 + n_states)
reverse_alphabet[missing_data] = missing_char
letter_to_state[missing_char] = missing_data
```

Example: 3 states "North America", "Oceania", "South America" -> internal `'A'`, `'B'`, `'C'`, missing `'D'`.

### Complete structure

```json
{
  "generated_by": { "program": "augur", "version": "25.0.0" },
  "models": {
    "region": {
      "rate": 3.737,
      "alphabet": ["North America", "Oceania", "South America", "?"],
      "equilibrium_probabilities": [0.349, 0.305, 0.346],
      "transition_matrix": [
        [0.0, 1.221, 1.921],
        [1.221, 0.0, 1.308],
        [1.921, 1.308, 0.0]
      ]
    }
  },
  "nodes": {
    "<name>": {
      "region": "South America",
      "region_confidence": { "South America": 0.95, "North America": 0.05 },
      "region_entropy": 0.286
    }
  },
  "branches": {
    "<name>": {
      "labels": { "region": "North America \u2192 South America" }
    }
  }
}
```

### Per-node field specifications

#### `<column>`

String. The inferred discrete state. Read from TreeTime's compressed sequence:

```python
node.__setattr__(field, letter_to_state[node.cseq[0]])
```

`node.cseq[0]` is the single-character internal code (mugration uses length-1 compressed sequences). `letter_to_state` maps to the display name.

#### `<column>_confidence`

Object mapping state name to probability. Present with `--confidence`. Built from `node.marginal_profile[0]` (posterior over real states, excluding missing):

```python
TINY = 1e-12
pdis = node.marginal_profile[0]  # shape (n_states,)
marginal = [(letter_to_state[tt.gtr.alphabet[i]], pdis[i]) for i in range(len(tt.gtr.alphabet))]
marginal.sort(key=lambda x: x[1], reverse=True)
marginal = [(a, b) for a, b in marginal if b > 0.001]
conf = {a: b for a, b in marginal}
```

Key: `tt.gtr.alphabet` is the internal letter array (size n_states, excludes missing). Values are numpy float64, serialized at full precision.

#### `<column>_entropy`

Float. Shannon entropy of the posterior. Present with `--confidence`:

```python
S = -np.sum(pdis * np.log(pdis + TINY))  # TINY = 1e-12
```

Computed over real states only (excludes missing).

### Top-level fields

#### `models.<column>`

Present when GTR model available (absent when 0 or 1 unique states). Built from TreeTime's GTR object:

```python
models[column]['rate'] = gtr.mu                                         # float: overall rate
models[column]['alphabet'] = [alphabet[k] for k in sorted(alphabet.keys())]  # state names (includes "?")
models[column]['equilibrium_probabilities'] = list(gtr.Pi)              # shape (n_states,), excludes missing
models[column]['transition_matrix'] = [list(x) for x in gtr.W]         # shape (n_states, n_states), excludes missing
```

Sizing asymmetry: `alphabet` has n_states+1 elements (includes `"?"` for missing data), `equilibrium_probabilities` and `transition_matrix` have n_states (exclude missing). The missing char sorts last in `alphabet` because it gets `chr(65 + n_states)`.

Not consumed by `augur export v2`.

#### `branches.<name>.labels`

Optional. Present when `--branch-labels` is used. Maps node name to `{"labels": {"<key>": "<label>"}}`.

Label format: `"<parent_state> <ARROW> <child_state>"` where `<ARROW>` is U+2192 RIGHTWARDS ARROW (the literal Unicode character, not ASCII `->` or `=>`). Produced by Python `f"{parent_state} \u2192 {child_state}"`. Root gets just the state name (no arrow). Duplicate labels disambiguated with counter suffix (space + number).

With `--branch-confidence`: below-threshold transitions labeled `"uncertain"`.

### Downstream consumption

| Node data field       | Auspice JSON target              |
| --------------------- | -------------------------------- |
| `<column>`            | `node_attrs.<column>.value`      |
| `<column>_confidence` | `node_attrs.<column>.confidence` |
| `<column>_entropy`    | `node_attrs.<column>.entropy`    |
| `models`              | Not consumed                     |
| `branches.*.labels`   | `branch_attrs.labels`            |

## Mutation string format

All mutation strings: `<ancestral_char><1-based-position><derived_char>`.

Nucleotide characters: `[ATCGNYRWSKMDVHB-]` (IUPAC + gap). Amino acid: `[A-Z*-]` (standard AA + stop + gap).

TreeTime v0 stores 0-based positions internally; augur converts to 1-based: `f"{a}{pos+1}{d}"`.

Treetime v1 `Sub::Display` (`packages/treetime/src/seq/mutation.rs:156-159`) already produces 1-based format: `write!(f, "{}{}{}", self.reff, self.pos + 1, self.qry)`. Direct `to_string()` is augur-compatible.

## Annotations schema

[`schema-annotations.json`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/data/schema-annotations.json) (`$id: https://nextstrain.org/schemas/augur/annotations`, draft-06).

The `nuc` entry: `start` must equal 1, `strand` must equal `"+"`. CDS entries (pattern `^(?!nuc$)[a-zA-Z0-9*_.()-]+$`): require `strand`, either simple `start`/`end` or compound `segments` array. Optional CDS properties: `gene`, `color`, `display_name`, `description`, `seqid`. `additionalProperties: true` on entries.

## Serialization format

Augur uses `json.dump` with `indent=2`, `sort_keys=True` (unless OrderedDict), and `AugurJSONEncoder` that converts numpy types to Python primitives. For treetime v1 to produce 1:1 matching output: `serde_json` with `indent=2` and sorted keys. Treetime v1's `json_write_file()` with `JsonPretty(true)` uses `indent=2` by default (via `serde_json::to_writer_pretty`), and serde's default for `BTreeMap` produces sorted keys.

## Rust type design

The format has a fixed envelope (`generated_by`, `nodes`) with command-specific contents in both the top level and per-node objects. Types live in a dedicated `packages/util-augur-node-data-json/` crate, following the `util-*` convention used by `util-phyloxml` and `util-usher-mat`. All types use the `AugurNodeDataJson` prefix. The design uses a generic container parameterized over command-specific types:

```rust
#[derive(Serialize, Deserialize)]
pub struct AugurNodeDataJson<M, N> {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub generated_by: Option<AugurNodeDataJsonGeneratedBy>,
    pub nodes: BTreeMap<String, N>,
    #[serde(flatten)]
    pub metadata: M,
}
```

Each command defines typed metadata and per-node structs:

```rust
pub type AugurNodeDataJsonAncestral = AugurNodeDataJson<AugurNodeDataJsonAncestralMeta, AugurNodeDataJsonAncestralNode>;
pub type AugurNodeDataJsonRefine = AugurNodeDataJson<AugurNodeDataJsonRefineMeta, AugurNodeDataJsonRefineNode>;
pub type AugurNodeDataJsonTraits = AugurNodeDataJson<AugurNodeDataJsonTraitsMeta, AugurNodeDataJsonTraitsNode>;
```

All structs carry a `#[serde(flatten)] pub other: BTreeMap<String, Value>` catch-all for forward compatibility. On write, `other` is empty. On read, unknown fields land in `other` and survive round-trip. This follows the same pattern used in `packages/treetime-io/src/auspice_types.rs`.

Future extension path: promote fields from `other` into typed fields one at a time by adding `Option<T>` or `#[serde(default)]` fields to the struct. Serde deserializes matching keys into the typed field instead of `other`. No format change, no breaking changes.

For generic read/merge (no typed access), `AugurNodeDataJson<BTreeMap<String, Value>, BTreeMap<String, Value>>` accepts any node data JSON. Merging replicates augur's [`deep_add_or_update()`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/util_support/node_data.py#L16-L44): dict+dict merges recursively, non-dict overwrites, dict/non-dict mismatch errors. `BTreeMap<String, Value>` covers all JSON extension cases (`serde_json::Value` handles strings, numbers, bools, null, nested objects, arrays).

## Treetime v1 data sources

### For `ancestral` command

| Node data field        | v1 source                                                       | Status    |
| ---------------------- | --------------------------------------------------------------- | --------- |
| `annotations.nuc.end`  | `partition.length` or `get_sequence_length()`                   | Available |
| `reference.nuc`        | `PartitionMarginalSparse.root_sequence` or dense root MAP       | Available |
| `mask`                 | Not implemented                                                 | Missing   |
| `nodes.*.muts`         | `PartitionBranchOps::edge_subs()` -> `Vec<Sub>`, `Sub::Display` | Available |
| `nodes.*.sequence`     | `SparseNodePartition.seq.sequence` or dense MAP state           | Available |
| `nodes.*.aa_muts`      | Not implemented (v1 has no AA reconstruction)                   | Missing   |
| `nodes.*.aa_sequences` | Not implemented                                                 | Missing   |

Partition access: inference results (mutations, sequences) live in the partition, not on graph payloads. Per-edge mutations are in `SparseEdgePartition.ml_subs`, per-node sequences in `SparseNodePartition.seq.sequence`, root sequence in `partition.root_sequence`, alignment length in `partition.length`. The partition is local to `run_ancestral_reconstruction()` but the existing output code (FASTA via `ancestral_reconstruction_marginal` callback at `run.rs:154`, Newick mutations via `MutationCommentProvider` at `run.rs:165`) already writes while the partition is in scope. Node data JSON writing goes alongside these existing writes.

### For `timetree` command

| Node data field       | v1 source                                    | Status                      |
| --------------------- | -------------------------------------------- | --------------------------- |
| `branch_length`       | `EdgeTimetree.base.branch_length`            | Available                   |
| `confidence`          | Input tree bootstrap value                   | Available                   |
| `numdate`             | `NodeTimetree.time`                          | Available                   |
| `clock_length`        | Derivable from clock model                   | Available                   |
| `mutation_length`     | From partition `edge_subs().len()` / seq_len | Available (needs partition) |
| `raw_date`            | Original date string from metadata           | Missing (not preserved)     |
| `date`                | Resolved date string                         | Missing (not preserved)     |
| `date_inferred`       | Derivable from date constraint type          | Available                   |
| `num_date_confidence` | `NodeConfidenceInterval.lower/upper`         | Available                   |
| `clock.rate`          | `ClockModel.clock_rate()`                    | Available                   |
| `clock.intercept`     | `ClockModel.intercept()`                     | Available                   |
| `clock.rtt_Tmrca`     | `-intercept / clock_rate`                    | Available                   |
| `clock.cov`           | `ClockModel.stats.cov`                       | Available                   |
| `clock.rate_std`      | `sqrt(cov[0][0])`                            | Available                   |

Partition access: timetree's `write_outputs()` (`commands/timetree/run.rs:529`) already receives the partition as `&[Arc<RwLock<dyn PartitionTimetreeAll>>]` and uses it for `MutationCommentProvider`. Node data JSON writing goes alongside the existing `write_graph_files_with`, `write_clock_model`, and `write_auspice_json` calls. Dates and clock data come from `NodeTimetree` payloads and `ClockModel`. Original date strings (`raw_date`, `date`) are not preserved in v1 (tracked: [../issues/M-dates-raw-string-not-preserved.md](../issues/M-dates-raw-string-not-preserved.md)).

### For `mugration` command

| Node data field                             | v1 source                                                        | Status                       |
| ------------------------------------------- | ---------------------------------------------------------------- | ---------------------------- |
| `<column>`                                  | `MugrationResult.traits.assignments`                             | Available                    |
| `<column>_confidence`                       | `MugrationResult.confidence.rows[].profile` + `partition.states` | Available (needs formatting) |
| `<column>_entropy`                          | Compute from confidence: `-sum(p * ln(p + 1e-12))`               | Available                    |
| `models.<column>.rate`                      | `partition.gtr().mu`                                             | Available                    |
| `models.<column>.alphabet`                  | `partition.states` + missing data marker                         | Available                    |
| `models.<column>.equilibrium_probabilities` | `partition.gtr().Pi` diagonal                                    | Available                    |
| `models.<column>.transition_matrix`         | `partition.gtr().W`                                              | Available                    |
| `branches.*.labels`                         | Compare parent/child trait assignments                           | Derivable                    |

`MugrationResult` carries the partition. Best-positioned command for node data JSON.

### Existing utilities

- `Sub::Display` (`packages/treetime/src/seq/mutation.rs:156-159`): produces `"A7G"` format directly
- `MutationCommentProvider` (`packages/treetime/src/partition/traits.rs:119-144`): retrieves sorted per-edge subs
- `json_write_file()` / `JsonPretty` (`packages/treetime-utils/src/io/json.rs`)
- `AuspiceGenomeAnnotations` and related types (`packages/treetime-io/src/auspice_types.rs`)

## Open work items

### Date metadata preservation

`read_date()` in `packages/treetime-io/src/dates_csv.rs` discards the original date string and collapses distinct input types into `DateOrRange`. Node data JSON requires `raw_date` (original string), `date` (resolved `"YYYY-MM-DD"`), and `date_inferred` (exact vs ambiguous).

Design: replace `DateOrRange` with `DateConstraint` carrying a `DateValue` enum (`Exact(DateExact)`, `Uncertain(DateRange)`, `Range(DateRange)`) and the raw input string. `date_inferred` derives from the variant. `date` derives from `year_fraction_to_date()` + formatting.

Tracked: [../issues/M-dates-raw-string-not-preserved.md](../issues/M-dates-raw-string-not-preserved.md), [../tickets/dates-preserve-raw-string-and-input-type.md](../tickets/dates-preserve-raw-string-and-input-type.md).

### Root profile sampling

v1 uses `argmax_first()` for sequence assignment at all nodes. v0 supports stochastic sampling at the root (`sample_from_profile='root'`), which augur always passes. This affects root mutations at positions with flat posteriors. Orthogonal to joint reconstruction (not related to joint removal).

Tracked: [../issues/M-ancestral-root-sampling-not-implemented.md](../issues/M-ancestral-root-sampling-not-implemented.md), [../tickets/ancestral-implement-root-profile-sampling.md](../tickets/ancestral-implement-root-profile-sampling.md).

### Mask computation

`create_mask()` identifies positions where all tips are ambiguous. Needed for two reasons: (1) mutation filtering - `collect_mutations()` skips mutations at masked positions to avoid reporting spurious inferred bases where no tip provides signal, (2) output field - `"mask"` string in ancestral node data JSON for VCF reconstruction.

Tracked: [../tickets/ancestral-implement-mask-computation.md](../tickets/ancestral-implement-mask-computation.md).

### AA reconstruction

Node data JSON from `augur ancestral --genes` includes per-gene AA mutations, root AA sequences, and gene annotations. This is not a separate algorithm: augur loops over genes, calling the same `run_ancestral()` with `alphabet='aa'` on pre-translated AA alignments. v1's reconstruction pipeline is alphabet-agnostic and supports AA alphabet + JTT92. The gap is input/output plumbing (gene annotation parsing, per-gene loop, CLI wiring). TreeTime v0 standalone does not do per-gene AA reconstruction either - this is augur-side orchestration.

Initial node data JSON implementation is nuc-only. Tracked: [../proposals/node-data-json-aa-reconstruction.md](../proposals/node-data-json-aa-reconstruction.md).

## Implementation notes

### Confidence formatting for mugration

The `<column>_confidence` dict requires: state-to-probability map, sorted descending by probability, filtered at > 0.001. In v1, `MugrationConfidenceOutput` already has the full profile. Formatting:

```
for each state i:
  if profile[i] > 0.001:
    emit (states[i], profile[i])
sort by probability descending
```

### Entropy computation for mugration

`-sum(p * ln(p + 1e-12))` over real states (exclude missing). The `1e-12` additive prevents `log(0)`.

### Alphabet sizing asymmetry for mugration models

`models.<column>.alphabet` includes the missing data marker (`"?"`) at the end. `equilibrium_probabilities` and `transition_matrix` exclude it. This is intentional: the missing state has no equilibrium frequency or transition rate in the model.

### Branch label arrow character

Unicode U+2192 `\u2192` (rightwards arrow), not ASCII `->`.

## Sources

All augur links pinned to commit `024292af6daf` (2026-05-20).

- augur ancestral: [`augur/ancestral.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/ancestral.py) - `run()`, `_to_ancestral_json()`, `collect_mutations()`, `collect_sequences()`, `create_mask()`, `reconstruct_translations()`
- augur refine: [`augur/refine.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/refine.py) - `refine()`, `collect_node_data()`, `run()`
- augur traits: [`augur/traits.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/traits.py) - `mugration_inference()`, `run()`, `BranchLabeller`
- augur JSON writing: [`augur/utils.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/utils.py#L123-L137) `write_augur_json()`, [`augur/io/json.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/io/json.py#L94-L166) `write_json()`, `AugurJSONEncoder`
- augur validation: [`augur/util_support/node_data_file.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/util_support/node_data_file.py), [`augur/util_support/node_data_reader.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/util_support/node_data_reader.py)
- augur deep merge: [`augur/util_support/node_data.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/util_support/node_data.py#L16-L44)
- augur export consumption: [`augur/export_v2.py`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/export_v2.py) - [`node_data_prop_is_normal_trait()`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/export_v2.py#L865-L898), `create_branch_mutations()`, `parse_node_data_and_metadata()`
- annotations schema: [`augur/data/schema-annotations.json`](https://github.com/nextstrain/augur/blob/024292af6daf/augur/data/schema-annotations.json)
- format docs: [`docs/usage/json_format.md`](https://github.com/nextstrain/augur/blob/024292af6daf/docs/usage/json_format.md)
- augur test fixtures: [`tests/functional/ancestral/data/`](https://github.com/nextstrain/augur/tree/024292af6daf/tests/functional/ancestral/data), [`tests/functional/refine/data/`](https://github.com/nextstrain/augur/tree/024292af6daf/tests/functional/refine/data), [`tests/functional/traits/`](https://github.com/nextstrain/augur/tree/024292af6daf/tests/functional/traits)
- treetime v1 Sub type: `packages/treetime/src/seq/mutation.rs`
- treetime v1 mutation access: `packages/treetime/src/partition/traits.rs`
- treetime v1 ancestral command: `packages/treetime/src/commands/ancestral/run.rs`
- treetime v1 payloads: `packages/treetime/src/payload/ancestral.rs`, `packages/treetime/src/payload/timetree.rs`
- workflow repos: [`nextstrain/seasonal-flu`](https://github.com/nextstrain/seasonal-flu), [`nextstrain/ncov`](https://github.com/nextstrain/ncov)
