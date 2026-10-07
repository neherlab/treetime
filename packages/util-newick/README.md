# util-newick

A reader and writer for phylogenetic trees and networks in Newick and its dialects, and for the tree blocks of NEXUS files. A [pest](https://pest.rs) grammar defines the syntax of every dialect. The Rust code builds the nesting of parentheses and BEAST arrays with explicit stacks, maps grammar tokens to the data model, and checks rules that span several nodes. Trees of any depth read and write on a small thread stack.

## Dialects

A dialect is a pair of two independent settings: the tree structure (`NewickStructure`) and the annotation convention (`NewickAnnotations`). The structure decides how labels and colon fields read; the annotation convention decides how comments read. Every structure combines with every annotation convention, so the crate reads and writes 12 dialects. Their syntax and the tools that use them are described in [kb/reports/newick-annotation-dialects.md](../../kb/reports/newick-annotation-dialects.md).

| Structure | Covers                                                                                | Labels          | Colon fields                        |
| --------- | ------------------------------------------------------------------------------------- | --------------- | ----------------------------------- |
| `Classic` | Standard Newick (Felsenstein 1986)                                                    | `label`         | branch length                       |
| `ENewick` | Extended Newick hybrid nodes `x#H1`, acceptor `x##LGT1` (Dendroscope, SplitsTree)     | `network_label` | branch length                       |
| `Rich`    | Rich Newick (PhyloNet): eNewick, three colon fields, `[&R]`/`[&U]`, tree weight `[&W` | `network_label` | branch length, support, probability |

| Annotations | Covers                                                                    | Node and branch comments (in `src/newick.pest`)                 |
| ----------- | ------------------------------------------------------------------------- | --------------------------------------------------------------- |
| `Plain`     | Comments as text                                                          | `plain_comment`                                                 |
| `Beast`     | BEAST 1 and 2, FigTree, TreeAnnotator, MrBayes consensus, IQ-TREE `[&..]` | `beast_comment`, then `malformed_annotation`, `plain_comment`   |
| `Nhx`       | New Hampshire eXtended `[&&NHX:..]` (Forester, ETE)                       | `nhx_comment`, then `malformed_annotation`, `plain_comment`     |
| `MrBayes`   | MrBayes sampling comments `[&E ..]`, `[&B ..]`, `[&N ..]`                 | `mrbayes_comment`, then `malformed_annotation`, `plain_comment` |

`NewickDialect { structure, annotations }` holds one pair. Constants name the common pairs: `CLASSIC` (classic, plain), `BEAST` (classic, beast), `NHX` (classic, nhx), `MRBAYES` (classic, mrbayes), `ENEWICK` (enewick, plain), `RICH` (rich, plain) and `ENEWICK_BEAST` (enewick, beast), the format of BEAST2 network packages such as CoalRe, Bacter and SpeciesNetwork, and of TreeKnit. `NewickDialect::pairs()` lists all 12 pairs. A dialect reads from and writes to the text `<structure>,<annotations>` (`classic,beast`, `enewick,beast`) through `FromStr` and `Display`.

Examples:

```text
classic,plain    (A:0.1,B:0.2,(C:0.3,D:0.4)95:0.5)root;
classic,beast    [&R]((A[&rate=1.5,s="x,y"]:[&length=1]0.1,B:0.2)[&posterior=0.9]:0.5,C:0.7);
classic,mrbayes  [&U]((1:0.1[&B TK02Brlens 0.1],2:0.1):0.2,3:0.3);
classic,nhx      ((A:0.1[&&NHX:S=human:E=1.1.1.1],B:0.2[&&NHX:S=mouse]):0.05[&&NHX:D=Y:B=100],C:0.3);
enewick,plain    (A,B,((C,(Y)x#H1)c,(x#H1,D)d)e)f;
enewick,beast    ((A[&segments={0,1}]:1,(B:1)#H0[&segments={0}]:0.5):1,(#H0[&segments={1}]:0.7,C:1):1);
rich,plain       [&U][&W 0.5]((A:1:90,(B:1)#H1:1::0.3):1,(#H1:1::0.7,C:1):1);
```

The grammar spells out the composite rules of every pair (`<structure>_<annotations>_tail`, `_items`, `_preamble`, `_strict`, `_tolerant`), because pest rules take no parameters; one pest parser per pair compiled 7 times slower. A test renders these rules from one template for every pair and checks that the grammar contains each of them.

Rules shared by the dialects:

- **Labels**: an unquoted label is any run of characters other than whitespace and `( ) [ ] , ; :`, and it does not start with `'`. A label in single quotes may contain any character; `''` stands for one quote. A quoted label is always a name, also when it looks like a number (`'95'`)
- **Support labels**: a label made of numbers separated by `/` (`95`, `0.950`, `80.5/95`) is the support of the branch above an internal node. On a leaf, every label is a name
- **Hybrid tags**: only the `ENewick` and `Rich` structures split a tag `#`, `##`, an optional type (`H`, `LGT`, `R`) and an index off the end of a label, so `A#x#H1` is the name `A#x` of the hybrid node `#H1`. The tag is read from the label alone; a single occurrence of `EPI_ISL#402124` is a hybrid node with these structures. With `Classic`, `#` is part of the name
- **Colon fields**: `Rich` reads up to three fields, branch length, support and inheritance probability (`x#H1:2.25::0.2`); the other structures read one field, the branch length. A field may be empty (`A::0.9`)
- **Numbers**: an optional sign, digits with an optional decimal point or a leading point (`.5`), and an optional exponent
- **Comments**: comments may appear before and after a label and before and after the value of each colon field. A plain comment may nest brackets; quotes in it are plain text. The annotation conventions `Beast`, `Nhx` and `MrBayes` reserve comments that start with `[&`: a comment that starts with `[&` and does not follow the convention's syntax is a malformed annotation. `Plain` reads every node and branch comment as text, with every structure
- **Tree comments**: before the first token of a tree, `[&R]` and `[&U]` set the rooting of the tree with `Beast` or `MrBayes` annotations or with `Rich` structure, and `[&W 0.5]` sets the tree weight with `Rich` structure. Every other comment before the first token belongs to the root node

BEAST annotations (`beast_comment`) are comma-separated pairs `key=value` or a bare `key`, which reads as `TRUE`. A key is bare or in double or single quotes. A value is a string in double quotes (`""` is a quote) or single quotes, a color `#RRGGBB`, `TRUE` or `FALSE`, a number, an array `{..}` of values that may nest, or a bare string without `, [ ] { }` that does not start with a quote.

NHX annotations (`nhx_comment`) are `:`-separated tags `key=value` or a bare `key`, which reads as `true`. A value has parts separated by `>`. The standard tags have types: `B` a decimal number, `T` and `W` integers, `D` one of `T`, `F`, `Y`, `N`, `?`, `C` a color `red.green.blue`. Every other tag holds text, or a list of text parts when it has more than one part. A number in a tag of type text reads back as text. Forester and ETE write the NHX comment after the branch length (`ADH2:0.1[&&NHX:S=human]`), so its tags are comments of the edge above the node; `NewickEdgeData::annotations()` returns them.

MrBayes sampling comments (`mrbayes_comment`) are `[&E name ..]`, `[&B name ..]` and `[&N name ..]`: a kind letter, a name and values separated by whitespace.

## Reading

```rust
use util_newick::{NewickDialect, NewickReadOptions, ReadMode, newick_from_str, newick_trees};

let options = NewickReadOptions {
  dialect: NewickDialect::ENEWICK_BEAST,
  mode: ReadMode::Tolerant,
  ..NewickReadOptions::default()
};
let tree = newick_from_str("((A[&rate=1.5]:1,(B)#H1[&gamma=0.3]:1),(#H1[&gamma=0.7],C));", &options)?;

for tree in newick_trees(std::io::stdin(), options) {
  let tree = tree?;
}
```

| Option                  | Default   | Meaning                                                                                                                        |
| ----------------------- | --------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `dialect`               | `CLASSIC` | The one dialect to read                                                                                                        |
| `mode`                  | `Strict`  | `Tolerant` also accepts a byte order mark, a missing final `;`, and malformed annotations                                      |
| `internal_label`        | `Auto`    | `Auto`: a support label of an internal node is support. `Name`: always a name. `Support`: every internal label must be support |
| `underscores_as_spaces` | `false`   | Read `_` in unquoted labels and NEXUS words as a space                                                                         |

The reader reads each tree once, in the given dialect. The same text means different things in different dialects (`EPI_ISL#402124` is a sample name with `Classic` structure and a hybrid node with `ENewick`; `[&rate=1]` is an annotation with `Beast` and a comment with `Plain`), so the caller states the dialect and the reader does not guess it.

Strict mode is an error on any malformed annotation or NHX tag value that does not fit its type. Tolerant mode keeps them as text and adds a `NewickWarning` with the line and column. A tree that reads in strict mode reads to the same graph, without warnings, in tolerant mode.

`newick_from_str()` and `newick_from_reader()` read exactly one tree and fail when another tree follows. `newick_trees()` reads tree by tree from any `Read`: it finds the `;` of each tree with a scanner rule, which stops before a comment, quoted label or BEAST string that is not complete yet and resumes there once more input has arrived. The scanner reads comments as the annotation convention does: with `Plain`, a comment ends at its first unnested `]`, so `(A,B)[&a="];(C,D)[x"];` is two trees; the other conventions read a double-quoted BEAST string in a `[&` comment as part of the comment. Comments between two trees belong to the root of the next tree; comments after the last tree are ignored. Input that is not valid UTF-8 is an error at the first invalid byte.

Errors (`NewickError`) carry the kind, the byte offset, the line and the column of the problem.

## Data model

`NewickGraph` holds the nodes, the edges, the root, the rooting (`[&R]`, `[&U]`), the tree weight, and the edge above the root (`root_edge`, with the root branch length, support and comments). Its fields are private: `NewickGraph::new()`, `add_child()` and `add_edge()` keep nodes and edges consistent, `add_edge()` refuses a second parent of a node that is not a hybrid node, and `validate()` checks the whole graph for cycles, for nodes that the root does not reach, and for comments of hybrid occurrences in the wrong place. `preorder()` and `postorder()` visit each node once without recursion.

- **Nodes** (`NewickNodeData`): a name, an optional hybrid tag (`NewickHybrid`) and the comments of the node with their position (`LabelSide::BeforeLabel`, `AfterLabel`). The copies of a hybrid node merge into one node with one parent per copy
- **Edges** (`NewickEdgeData`): branch length, support values with their source (`SupportSource::Label` or `Field`), inheritance probability, the acceptor flag of `##`, the comments of each colon field with their position (`EdgeField::Length`, `Support`, `Probability`; `ValueSide::BeforeValue`, `AfterValue`), and the occurrence comments. Support belongs to the edge: the support label of an internal node is stored on the edge above it
- **Occurrence comments**: a hybrid node is written once under each of its parents, and each copy may carry its own comments, such as CoalRe's `segments` or the SpeciesNetwork inheritance probability `gamma` of that parent branch. The comments before and after the label of each copy are stored on the edge into that copy (`NewickEdgeData::occurrence_comments()`, on the root edge for a hybrid root), and the writer writes them back at the same copy. A hybrid node itself has no comments
- **Comments** (`NewickComment`): `Beast` and `Nhx` keep every pair in input order, repeated keys included; `annotations()` iterates the pairs of a node or an edge (occurrence comments first, then field comments), and `annotation(key)` returns the last value of a key. `MrBayesMcmc` keeps the kind, the name and the values; `Plain` keeps the text between the brackets
- **Values** (`NewickValue`): `Boolean`, `Number`, `NumberText`, `String`, `Color` and `Array`. A number whose text is the shortest text of its value reads as `Number`; any other number text (`2020.50`, `0123`, `1e5`, integers above 2^53) reads as `NumberText` and keeps its text. Nested arrays are freed without recursion

Equality compares every field; `==` ignores the order of children and `eq_ordered()` does not. The reader produces one form for data that has several spellings, and a graph built in code reads back equal only in that form:

- A node without a name and a hybrid tag has its comments `AfterLabel`
- A colon field without a value has its comments `BeforeValue`
- `NumberText` holds a text that is not the shortest text of its value

## Writing

```rust
use util_newick::{NewickDialect, NewickWriteOptions, newick_to_string};

let text = newick_to_string(&tree.graph, &NewickWriteOptions::new(NewickDialect::ENEWICK_BEAST))?;
```

| Option               | Default       | Meaning                                                                                     |
| -------------------- | ------------- | ------------------------------------------------------------------------------------------- |
| `dialect`            | `CLASSIC`     | The one dialect to write                                                                    |
| `branch_annotations` | `Recorded`    | `BeforeLength` or `AfterLength` move the annotations of the branch-length field to one side |
| `support`            | `Source`      | Write support where it was read, or as a `Label`, a Rich `Field`, or an `Annotation(key)`   |
| `root_edge`          | `true`        | Write the edge above the root                                                               |
| `quoting`            | `WhenNeeded`  | `Always` quotes every name                                                                  |
| `spaces`             | `Quote`       | `Underscore` writes spaces in names as `_`                                                  |
| `numbers`            | shortest text | `NumberFormat`: significant digits, decimal digits, and `point_zero` (`1.0` instead of `1`) |
| `indent`             | `None`        | Write one node per line, indented by the given number of spaces per level                   |

With the default options a write followed by a read in the same dialect gives an equal graph (property tests check this for all 12 dialects, networks, occurrence comments and every comment kind). Names are quoted exactly when the grammar would read them differently unquoted: a name that is not a run of label characters, that contains `'`, that ends like a hybrid tag, or, on an internal node, that reads as a support label. BEAST keys are quoted when they are not bare keys, and BEAST strings are always written in double quotes, so `"TRUE"` stays a string. A NHX tag or value that NHX cannot hold, and a number that is not finite, is an error.

Writer errors, and errors of `add_edge()`, `add_child()` and `validate()`, are `NewickWriteError`, with one message that starts with the context of the step that failed: "When writing Newick: When writing the branch above node 1 ('A'): Newick cannot represent the number inf". A dialect name that `FromStr` cannot read is a `ParseDialectError` ([kb/decisions/util-newick-typed-errors.md](../../kb/decisions/util-newick-typed-errors.md)).

`write_newick_trees()` writes one tree per line.

### Data a dialect cannot hold

`conversion(dialect, DataKind)` is the one table that decides what the writer does with data the chosen dialect cannot hold. Each row follows the setting that owns the data:

| Data                                                 | Held by                                            | Otherwise |
| ---------------------------------------------------- | -------------------------------------------------- | --------- |
| BEAST annotations                                    | `Beast` annotations                                | dropped   |
| NHX annotations                                      | `Nhx` annotations                                  | dropped   |
| MrBayes sampling comments                            | `MrBayes` annotations                              | dropped   |
| Plain comments that start with `&`                   | `Plain` annotations                                | error     |
| Plain comments `[&R]`, `[&U]` before the first token | dialects without a rooting comment                 | error     |
| Plain comments `[&W x]` before the first token       | structures without a tree weight                   | error     |
| Rooting `[&R]`, `[&U]`                               | `Beast` or `MrBayes` annotations, `Rich` structure | dropped   |
| Tree weight                                          | `Rich` structure                                   | dropped   |
| Support from a colon field                           | `Rich` structure                                   | dropped   |
| Inheritance probability                              | `Rich` structure                                   | dropped   |
| Comments of the support and probability fields       | `Rich` structure                                   | dropped   |
| Hybrid nodes                                         | `ENewick`, `Rich` structures                       | error     |

A plain comment before the first token is the comment of a root without children; written as `[&R]` with `Rich` structure, it would read back as the rooting of the tree. The options `root_edge: false`, digit limits and `spaces: Underscore` are conversions as well. Support written as `Annotation(key)` reads back as an annotation, not as support. Data that an option cannot express is an error that names the node: support on a leaf or on a named node with `support: Label`, more than one support value in a colon field.

## NEXUS

`nexus_from_str()`, `nexus_from_reader()` and `nexus_trees()` read NEXUS files command by command (`src/nexus.pest`), in the dialect of `NewickReadOptions::dialect`:

- **Trees blocks**: `Tree` and `UTree` commands with an optional `*`, quoted names, and comments between the name and `=` (`tree STATE_0 [&lnP=-3195.24] = ...`), read with the annotation convention of the dialect. `Translate` maps leaf labels to taxon names, per block. `PROPERTIES rooted=yes` sets the rooting of the following trees of the block that have no rooting comment
- **Taxa blocks**: `TaxLabels` and `Dimensions NTax`. When a taxon set is declared, a leaf label must be a taxon label, a `Translate` key or a taxon number; otherwise it is an error in strict mode and a warning in tolerant mode
- **Other commands and blocks** are skipped and returned in `NexusFile::skipped` with their block and line
- **Command ends**: a command ends at the first `;` outside quoted words and comments. A comment that starts with `[&` may hold a double-quoted BEAST string with `]` or `;` in every dialect, because the comment rule is global to the NEXUS grammar and NEXUS tree files come from BEAST-family tools
- Strict mode is an error on a tree or `Translate` command that does not follow its syntax, a command outside a block, and a block without its end. Tolerant mode also reads unquoted names with spaces (`tree my tree = ...`, `Translate 1 Homo sapiens`)

`nexus_to_string()` writes a `Taxa` block with the leaf names of all trees, an optional `Translate` table with numbered labels (`NexusWriteOptions::translate`), and the trees in the dialect of `NexusWriteOptions::newick`. Words are quoted when they contain NEXUS punctuation.

## Tests and benchmark

`src/__tests__/` tests every grammar rule, the reader and writer of each dialect, reading across buffer boundaries, NEXUS, and the conversion table. Property tests round-trip generated graphs in all 12 dialects and feed random and mutated input to the reader in every dialect on a 2 MiB stack. The fixtures modeled on BEAST, TreeAnnotator, MrBayes, IQ-TREE, FastTree, RAxML-NG, Forester, Dendroscope, SplitsTree and PhyloNet output are compared with the readings of DendroPy 5.0.8, ETE3 3.1.3 and PhyloNet 3.8.5; cases named `spec_` take their expectation from the dialect definitions, where these tools reject or lose part of the input, such as the eNewick networks with BEAST annotations of TreeKnit.

`benches/read_dialects.rs` measures reading an annotated BEAST tree with `BEAST` and an annotated network with `ENEWICK_BEAST`.
