# Newick annotation dialects: specification, grammar, and tool interoperability

Research report covering all Newick-family annotation formats relevant to implementing a unified reader/writer. Compiled from primary sources: original papers, official format specifications, and source code of 9 phylogenetic tool repositories.

## Background

The <a id="gloss-use-3"></a>Newick <sup>[3](#gloss-3)</sup> format <a id="cite-1"></a>[Felsenstein 1986](https://phylipweb.github.io/phylip/newicktree.html) [[1](#ref-1)] encodes phylogenetic tree topology, branch lengths, and node labels in a single parenthetical string. It has no metadata mechanism. Over four decades, tools added incompatible annotation conventions by repurposing the comment syntax `[...]` inherited from Newick and <a id="gloss-use-1"></a>NEXUS <sup>[1](#gloss-1)</sup> <a id="cite-2a"></a>[Maddison, Swofford, and Maddison 1997](https://doi.org/10.1093/sysbio/46.4.590) [[2](#ref-2)]. Each convention defines different delimiters, value types, and escaping rules.

The NEXUS specification provides no per-node or per-branch annotation mechanism. Comments in NEXUS are formally discardable: a spec-compliant processor may strip, reorder, or drop them. The <a id="gloss-use-2"></a>NeXML <sup>[2](#gloss-2)</sup> paper <a id="cite-3a"></a>[Vos et al. 2012](https://doi.org/10.1093/sysbio/sys025) [[3](#ref-3)] states this explicitly:

> "the NEXUS standard (which integrates Newick) does not specify that comments in an input stream must be retained in an output stream, much less retained at their original position. Thus, a valid NEXUS processor may scramble a valid NHX tree."

The only tree-level metadata NEXUS formally defines are the `[&R]` (rooted) and `[&U]` (unrooted) command comments and the `PROPERTIES rooted=yes` command in the TREES block [[src](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/ncl/nxstreesblock.cpp)]. No authoritative NEXUS revision exists after 1997. The NESCent working group concluded the format cannot be fixed without breaking backward compatibility and built NeXML instead <a id="cite-3b"></a>[Vos et al. 2012](https://doi.org/10.1093/sysbio/sys025) [[3](#ref-3)].

Every annotation scheme described below is layered on top of Newick/NEXUS by individual tools, encoded as structured comments inside the Newick string. NEXUS is just a container; the `.nexus` extension provides no annotation advantage over `.nwk`.

## Base format: standard Newick

Adopted 1986 at Newick's restaurant in Dover, New Hampshire. Generalizes notation developed by Christopher Meacham in 1984 for PHYLIP <a id="cite-1b"></a>[Felsenstein 1986](https://phylipweb.github.io/phylip/newicktree.html) [[1](#ref-1)]. No formal specification document exists; the grammar below is reconstructed from the [Felsenstein description](https://phylipweb.github.io/phylip/newicktree.html) [[1](#ref-1)].

### Grammar

```
Tree      -> Subtree ';'
Subtree   -> Leaf | Internal
Leaf      -> Name
Internal  -> '(' BranchSet ')' Name
BranchSet -> Branch | Branch ',' BranchSet
Branch    -> Subtree Length
Name      -> empty | string
Length    -> empty | ':' number
```

### Label rules

- Unquoted labels: may not contain `( ) [ ] , ; :` or whitespace. Underscores map to spaces on read
- Quoted labels: enclosed in single quotes `'...'`. Internal single quotes are doubled: `'it''s'` represents `it's`
- A label parseable as a number is ambiguous (could be a name or a weight). The `bio` crate and many parsers interpret numeric labels as names; some tools interpret them as bootstrap values

### Comments

`[...]` brackets denote comments. The base Newick grammar treats them as discardable whitespace. Comments starting with `&` are "command comments" in some tools but this is not part of the base format.

### What standard Newick carries

- Tree topology (parenthetical nesting)
- Branch lengths (`:number` after each subtree)
- Node labels (string before `:` or before `)`)

### What standard Newick does not carry

- Per-node or per-branch metadata (no key-value mechanism)
- Rooted vs unrooted distinction (no marker; convention is 2 children = rooted, 3 = unrooted)
- Node vs branch attribute distinction (no separate label positions for node and edge data)
- Any structured annotation beyond topology, labels, and branch lengths

### Examples

```
(,,(,));                                          no names
(A,B,(C,D));                                      leaf names
(A,B,(C,D)E)F;                                    all names
(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);                 distances and leaf names
(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;               distances and all names
```

## Dialect 1: BEAST/FigTree `[&key=value]`

Origin: <a id="cite-4"></a>[BEAST1 (Suchard et al. 2018)](https://doi.org/10.1093/ve/vey016) [[4](#ref-4)] [[repo](https://github.com/beast-dev/beast-mcmc)], adopted by <a id="cite-6"></a>[MrBayes (Ronquist and Huelsenbeck 2003)](https://doi.org/10.1093/bioinformatics/btg180) [[6](#ref-6)] [[repo](https://github.com/NBISweden/MrBayes)], FigTree [[repo](https://github.com/rambaut/figtree)], <a id="cite-7a"></a>[IQ-TREE (Nguyen et al. 2015)](https://doi.org/10.1093/molbev/msu300) [[7](#ref-7)] [[repo](https://github.com/Cibiv/IQ-TREE)]. De facto standard in Bayesian phylogenetics. Formalized in <a id="cite-5"></a>[BEAST2 (Bouckaert et al. 2014)](https://doi.org/10.1371/journal.pcbi.1003537) [[5](#ref-5)] [[repo](https://github.com/CompEvol/beast2)] via an ANTLR4 grammar.

### Grammar

Opener `[&`, closer `]`. Key-value pairs separated by commas. Format inside the brackets:

```
meta         = '[&' attrib (',' attrib)* ']'
attrib       = key '=' value | key           // bare key = boolean TRUE
key          = quoted_string | unquoted_key
value        = array | color | boolean | number | quoted_string | bare_string
array        = '{' value (',' value)* '}'
quoted_string = '"' [^"]* '"' | "'" [^']* "'"
color        = '#' [0-9A-Fa-f]{6}
boolean      = 'TRUE' | 'FALSE'              // case-insensitive
unquoted_key = [^,=\s]+
bare_string  = [^,\]]+
```

Source: BEAST1 regex at [src/dr/evolution/io/NexusImporter.java](https://github.com/beast-dev/beast-mcmc/blob/248263332d365fb97b1df0f5c2d4e28debdaa804/src/dr/evolution/io/NexusImporter.java) [[src](https://github.com/beast-dev/beast-mcmc/blob/248263332d365fb97b1df0f5c2d4e28debdaa804/src/dr/evolution/io/NexusImporter.java)]:

```java
Pattern.compile("(\"[^\"]*\"+|[^,=\\s]+)\\s*(=\\s*(\\{[^=]*\\}|\"[^\"]*\"+|[^,]+))?");
```

BEAST2 formalizes this in an ANTLR4 grammar at [NewickParser.g4](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/tree/treeparser/NewickParser.g4) [[src](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/tree/treeparser/NewickParser.g4)]:

```antlr
post: label? nodeMeta=meta? (':' lengthMeta=meta? length=number)? ;
meta: '[&' attrib (ACOMMA attrib)* ']' ;
attrib: attribKey=ASTRING '=' attribValue ;
attribValue: attribNumber | ASTRING | vector ;
vector: '{' attribValue (ACOMMA attribValue)* '}' ;
```

### Value types

BEAST1 `fn parseValue()` [[src](https://github.com/beast-dev/beast-mcmc/blob/248263332d365fb97b1df0f5c2d4e28debdaa804/src/dr/evolution/io/NexusImporter.java)] tries types in this order: `{...}` array, `#RRGGBB` color, `TRUE`/`FALSE` boolean, integer, double, quoted string, plain string. Nested arrays supported in BEAST1 (`{{1,2},{3,4}}`), not in BEAST2.

### Node vs branch annotation placement

Three conventions exist in practice:

```
BEAST2:      name[&NODE_ATTRS]:[&BRANCH_ATTRS]length     // branch meta between : and length
MrBayes:     name[&NODE_ATTRS]:length[&BRANCH_ATTRS]     // branch meta after length
BEAST1:      name[&ATTRS]:length[&ATTRS]                 // no distinction, both on node object
```

BEAST2 is the canonical grammar (ANTLR `nodeMeta` before `:`, `lengthMeta` after `:` before `length`). MrBayes writes node attrs before `:` and branch attrs (lengths, HPD intervals) after length [[src](https://github.com/NBISweden/MrBayes/blob/bb09fffbf9967cba62a4e8af2d487f443e9438d9/src/sumpt.c)]. BEAST1 does not distinguish at the parser level - both go to the same `FlexibleNode` object.

### Comma-in-value escaping

BEAST uses `{v1,v2}` for arrays (commas inside braces are not key-value separators). IQ-TREE uses double-quote escaping: `[&gCF="33.33",gDF1/gDF2="0/33.33"]` <a id="cite-7b"></a>[Nguyen et al. 2015](https://doi.org/10.1093/molbev/msu300) [[7](#ref-7)]. No standard escaping mechanism exists - this is a known compatibility problem.

### MrBayes MCMC format

During MCMC sampling (not consensus output), MrBayes uses a different space-separated format: `[&E paramname count: ...]`, `[&B paramname value]`, `[&N paramname value]` [[src](https://github.com/NBISweden/MrBayes/blob/bb09fffbf9967cba62a4e8af2d487f443e9438d9/src/utils.c)]. This is incompatible with BEAST-style parsing. The consensus output (`.con.tre`) switches to BEAST-style `[&prob=...,height_mean=...]` for FigTree compatibility.

### FigTree NEXUS block

FigTree [[repo](https://github.com/rambaut/figtree)] adds a custom `begin figtree;` block to NEXUS files for display settings (branch colors, line widths, label formatting). This uses the NEXUS extensible block mechanism (custom blocks are allowed by the spec and ignored by other parsers). Per-node scientific data (posteriors, heights, rates) is carried in `[&...]` annotations inside the Newick string, not in this block [[src](https://github.com/rambaut/figtree/blob/24c51adc9b4a61760b828d99511cf3449a7ba9d3/src/figtree/application/FigTreeNexusExporter.java)].

### Rooting markers

`[&R]` (rooted) and `[&U]` (unrooted) appear as a prefix before the tree string. `[&W value]` encodes tree weight/probability. These are parsed specially, not as key-value pairs.

## Dialect 2: NHX `[&&NHX:key=value:...]`

Origin: Christian Zmasek (ATV/Forester [[repo](https://github.com/cmzmasek/forester)]). NHX v2.0 is the final version (2008); Zmasek recommends phyloXML instead. <a id="cite-8"></a>[Zmasek and Eddy 2001](https://doi.org/10.1093/bioinformatics/17.4.383) [[8](#ref-8)].

### Grammar

```
nhx_comment  = '[&&NHX' (':' tag)* ']'
tag          = key '=' value
key          = predefined_tag_name         // S, T, B, D, GN, AC, E, etc.
value        = string                      // all values are strings
```

Opener: `[&&NHX` (double ampersand distinguishes from BEAST). Closer: `]`. Separator between tags: `:` (colon, not comma). Multi-part values use `>` as internal separator (avoids collision with Newick `,`).

Source: Official NHX v2.0 spec at phylosoft.org/NHX (Zmasek 2008). Parser at [NHXParser.java](https://github.com/cmzmasek/forester/blob/843bf03568c99429994b370e3451b2716b624e8d/forester/java/src/org/forester/io/parsers/nhx/NHXParser.java) [[src](https://github.com/cmzmasek/forester/blob/843bf03568c99429994b370e3451b2716b624e8d/forester/java/src/org/forester/io/parsers/nhx/NHXParser.java)].

### Standard tags

| Tag  | Type                      | Description                                |
| ---- | ------------------------- | ------------------------------------------ |
| `S`  | string                    | Species name                               |
| `T`  | integer                   | Taxonomy ID                                |
| `B`  | decimal                   | Bootstrap confidence for parent branch     |
| `D`  | `T`, `F`, or `?`          | Duplication event                          |
| `GN` | string                    | Gene name                                  |
| `AC` | string                    | Sequence accession                         |
| `E`  | string                    | EC number                                  |
| `Ev` | `int>int>int>str>str`     | Compound event descriptor                  |
| `DS` | `int>int>int>dbl>str>...` | Domain structure                           |
| `W`  | integer                   | Branch width                               |
| `C`  | `rrr.ggg.bbb`             | Branch color (RGB integers, dot-separated) |
| `XN` | string                    | Custom node data                           |
| `XB` | string                    | Custom branch data                         |

Source: [NHXtags.java](https://github.com/cmzmasek/forester/blob/843bf03568c99429994b370e3451b2716b624e8d/forester/java/src/org/forester/io/parsers/nhx/NHXtags.java) [[src](https://github.com/cmzmasek/forester/blob/843bf03568c99429994b370e3451b2716b624e8d/forester/java/src/org/forester/io/parsers/nhx/NHXtags.java)].

### Differences from BEAST-style

| Aspect         | BEAST `[&...]`                          | NHX `[&&NHX:...]`                    |
| -------------- | --------------------------------------- | ------------------------------------ |
| Opener         | `[&` (single `&`)                       | `[&&NHX` (double `&&`)               |
| Tag separator  | `,` (comma)                             | `:` (colon)                          |
| Value escaping | `{...}` for arrays, `"..."` for strings | `>` as internal separator            |
| Schema         | Ad-hoc, arbitrary keys                  | Fixed tag vocabulary                 |
| Node vs branch | Position-based (before/after `:`)       | Tag-based (`XN=` node, `XB=` branch) |

### Cross-format parsing

Forester's NHXParser recognizes both formats [[src](https://github.com/cmzmasek/forester/blob/843bf03568c99429994b370e3451b2716b624e8d/forester/java/src/org/forester/io/parsers/nhx/NHXParser.java)]: it checks for `&&NHX` prefix first, then falls back to BEAST-style parsing for `[&prob=...]` (MrBayes patterns), then stores unrecognized `[&...]` as raw comment strings.

## Dialect 3: Extended Newick (eNewick)

Origin: Morin and Moret 2006 (NetGen), formalized by <a id="cite-9"></a>[Cardona, Rossello, and Valiente 2008](https://doi.org/10.1186/1471-2105-9-532) [[9](#ref-9)].

### Grammar

Extends Newick to encode phylogenetic networks (DAGs) by duplicating hybrid nodes annotated with `#`:

```
Leaf    -> Name Hybrid
Hybrid  -> empty | '#' Type Integer
Type    -> empty | String              // H, LGT, R, or any string
```

Each hybrid node appears once per parent in the string. Parser merges nodes sharing the same `#` identifier into a single node with multiple parents.

### Reticulation types

| Type    | Event                 |
| ------- | --------------------- |
| `H`     | Hybridization         |
| `LGT`   | Lateral gene transfer |
| `R`     | Recombination         |
| (empty) | Unspecified           |

### Acceptor edge (`##`)

For LGT events, `##` (double hash) on exactly one copy marks the acceptor (main lineage) edge. Other copies are transfer edges. Convention used by Dendroscope [[repo](https://github.com/husonlab/dendroscope3)] [[src](https://github.com/husonlab/dendroscope3/blob/c2d35555003a9261120a1d3223fd00e3e180f46b/src/dendroscope/util/NewickInputDialog.java)] and SplitsTree [[repo](https://github.com/husonlab/splitstree6)] [[src](https://github.com/husonlab/splitstree6/blob/8236ada348055d9d6e7c4141c0c73511f38b075b/src/main/java/splitstree6/io/nexus/TreesNexusOutput.java)].

### Branch lengths on reticulate edges

Each copy carries its own branch length (the edge from that copy's parent to the hybrid node).

### Example

```
(A,B,((C,(Y)x#H1)c,(x#H1,D)d)e)f;
```

Node `x#H1` appears twice. Parser merges both into a single hybrid node with two parents (from lineage c and lineage d).

### NEXUS interaction

eNewick strings can appear inside NEXUS `Begin Trees;` blocks. SplitsTree adds `PROPERTIES reticulated=yes` as a custom extension to signal network content. Legacy parsers treat `x#H1` as an oddly-named leaf (backward compatible, but network semantics are lost).

### Backward compatibility

`#` is not in the standard forbidden character set for unquoted Newick labels in most implementations, so `x#H1` parses as a valid leaf name. This provides backward compatibility: the tree structure is readable, but the network is flattened to a tree with duplicate "leaves."

## Dialect 4: Rich Newick

Origin: PhyloNet [[repo](https://github.com/phylonet/PhyloNet)] (Rice University, Nakhleh Lab). Extends eNewick with rooting markers, tree probabilities, and inheritance probabilities.

### Grammar

From PhyloNet ANTLR grammar at [RichNewick_1_1.g](https://github.com/phylonet/PhyloNet/blob/2ee2b48cdf0d9fef61066d265c7690f2bca6fc2a/src/edu/rice/cs/bioinfo/library/language/richnewick/_1_1/reading/parsers/antlr/ast/) [[src](https://github.com/phylonet/PhyloNet/blob/2ee2b48cdf0d9fef61066d265c7690f2bca6fc2a/src/edu/rice/cs/bioinfo/library/language/richnewick/_1_1/reading/parsers/antlr/ast/)]:

```
ROOTAGE_QUALIFIER  = '[&' ('R'|'r'|'U'|'u') ']'
TREE_PROB          = '[&' ('W'|'w') ' '+ DECIMAL ']'
```

### Rooting markers

`[&R]` (rooted) and `[&U]` (unrooted) prefix before the tree string. Uses single `&` (same as BEAST comment prefix).

### Extra colon-delimited fields

Rich Newick extends per-node fields beyond `name:length`:

```
name:length:bootstrap:probability
```

For hybrid nodes with inheritance probabilities:

```
name#Hindex:branch_length:population_size:gamma
```

Where gamma (0 < gamma < 1) specifies the fraction of genome inherited from this parent edge.

Empty fields use `::` (colons present, values absent). This can break standard Newick parsers that attempt to parse extra colons as branch lengths.

### Relationship to other dialects

```
Standard Newick  ⊂  Extended Newick (eNewick)  ⊂  Rich Newick
                     (+#H markers)                 (+[&R/U] +[&W] +extra colon fields +gamma)
```

## Tool interoperability matrix

R = reads, W = writes, RW = both, blank = unsupported. Source evidence in footnotes.

| Tool        | Plain Newick | BEAST `[&k=v]` | NHX `[&&NHX:k=v]` | eNewick `#H` | Rich Newick |
| ----------- | ------------ | -------------- | ----------------- | ------------ | ----------- |
| BEAST1      | RW           | RW [a]         |                   |              |             |
| BEAST2      | RW           | RW [b]         |                   |              |             |
| MrBayes     | R            | W [c]          |                   |              |             |
| FigTree     | R            | RW [d]         |                   |              |             |
| IQ-TREE     | RW           | RW [e]         |                   |              |             |
| RAxML       | RW           |                |                   |              |             |
| FastTree    | RW           |                |                   |              |             |
| Dendroscope | RW           | R              |                   | RW [f]       |             |
| SplitsTree6 | RW           |                |                   | RW [f]       |             |
| PhyloNet    | RW           |                |                   | RW [g]       | RW [g]      |
| PAUP\*      | RW           |                |                   |              |             |
| Mesquite    | RW           |                |                   |              |             |
| ETE toolkit | RW           | R              | RW [h]            |              |             |
| BioPython   | RW           | RW [i]         |                   |              |             |
| DendroPy    | RW           | RW             | RW                |              |             |
| TreeTime v0 | RW           | W [j]          |                   |              |             |
| augur       | RW           |                |                   |              |             |

[a] `setCommentDelimiters('[', ']', '\0', '!', '&')` in [Importer.java](https://github.com/beast-dev/beast-mcmc/blob/248263332d365fb97b1df0f5c2d4e28debdaa804/src/dr/evolution/io/Importer.java). `parseMetaCommentPairs()` in [NexusImporter.java](https://github.com/beast-dev/beast-mcmc/blob/248263332d365fb97b1df0f5c2d4e28debdaa804/src/dr/evolution/io/NexusImporter.java). Supports `{v1,v2}` array values and nested arrays.
[b] ANTLR4 grammar at [NewickParser.g4](https://github.com/CompEvol/beast2/blob/9321c88df6e5e90da4db5c0fc4872576990c015c/src/beast/base/evolution/tree/treeparser/NewickParser.g4). Formally distinguishes `nodeMeta` from `lengthMeta`.
[c] Writes BEAST-style `[&prob=...,height_mean=...]` in consensus output via `PrintFigTreeNodeInfo` in [sumpt.c](https://github.com/NBISweden/MrBayes/blob/bb09fffbf9967cba62a4e8af2d487f443e9438d9/src/sumpt.c). Cannot read annotated output back in the general case.
[d] Reads/writes via JEBL library. Adds custom `begin figtree;` NEXUS block for display settings at [FigTreeNexusExporter.java](https://github.com/rambaut/figtree/blob/24c51adc9b4a61760b828d99511cf3449a7ba9d3/src/figtree/application/FigTreeNexusExporter.java).
[e] NCL layer reads `[&R]`/`[&U]` only in [nxstreesblock.cpp](https://github.com/Cibiv/IQ-TREE/blob/6776a95f15a2eccda2aa330497291dc246575995/ncl/nxstreesblock.cpp). Writes `[&gCF=...,sCF=...]` with double-quoted values for comma protection.
[f] Huson lab tools. Dendroscope: [NewickInputDialog.java](https://github.com/husonlab/dendroscope3/blob/c2d35555003a9261120a1d3223fd00e3e180f46b/src/dendroscope/util/NewickInputDialog.java). SplitsTree: [TreesNexusOutput.java](https://github.com/husonlab/splitstree6/blob/8236ada348055d9d6e7c4141c0c73511f38b075b/src/main/java/splitstree6/io/nexus/TreesNexusOutput.java) adds `PROPERTIES reticulated=yes`.
[g] [RnNewickPrinter.java](https://github.com/phylonet/PhyloNet/blob/2ee2b48cdf0d9fef61066d265c7690f2bca6fc2a/src/edu/rice/cs/bioinfo/library/language/richnewick/_1_1/reading/parsers/antlr/ast/) writes Rich Newick with `[&U]` prefix and `name#Hindex:bl:support:probability` fields.
[h] Reads NHX natively, BEAST-style with appropriate settings. Primary NHX consumer in the Python ecosystem.
[i] `Bio.Phylo` `.comment` attribute becomes `[&...]` in Nexus output. Does NOT parse `[&...]` back into structured data on read - comments stripped or stored as raw strings.
[j] Writes BEAST-style `[&mutations=...,date=...]` via BioPython `.comment` at [CLI_io.py](https://github.com/neherlab/treetime/blob/master/treetime/CLI_io.py). Cannot read annotations back (BioPython strips comments on Newick read). Local v0 source: `packages/legacy/treetime/treetime/CLI_io.py#L166-L205`.

### File extensions

| Extension                | Typical content                                     |
| ------------------------ | --------------------------------------------------- |
| `.nwk`, `.newick`, `.nh` | Plain or annotated Newick (ambiguous)               |
| `.nex`, `.nexus`, `.nxs` | NEXUS (may contain annotated Newick in TREES block) |
| `.trees`                 | NEXUS with multiple trees (BEAST posterior samples) |
| `.tree`, `.tre`          | Newick single tree                                  |
| `.treefile`, `.contree`  | IQ-TREE output (may contain `[&...]`)               |
| `.nhx`                   | NHX format (rare in practice)                       |
| `.enwk`                  | Extended Newick (very rare, PhyloNet)               |

No extension reliably distinguishes annotated from plain Newick. Auto-detection must be content-based, not extension-based.

## Common workflow chains

1. Bayesian phylogenetics: IQ-TREE/RAxML (`.nwk`) -> BEAST (`.trees` NEXUS) -> TreeAnnotator (`.tree` NEXUS) -> FigTree. All stages use BEAST-style `[&k=v]`
2. Nextstrain/augur: IQ-TREE/FastTree (`.nwk`) -> augur refine (calls TreeTime v0) -> augur export -> Auspice JSON. TreeTime v0 writes NEXUS with `[&mutations=...,date=...]` annotations
3. Gene tree reconciliation: RAxML (`.nwk`) -> Notung/RANGER (NHX) -> Forester/ETE. NHX with `[&&NHX:D=Y:S=human:B=90]`
4. Network phylogenetics: PhyloNet (Rich Newick) <-> SplitsTree/Dendroscope (eNewick). `#H`/`#LGT`/`#R` markers, `[&U]`/`[&R]` prefix

## Known compatibility problems

1. Comment stripping: BioPython `Phylo.read(file, 'newick')` discards `[...]` content. TreeTime v0 writes annotated Nexus but cannot read annotations back. The `bio` Rust crate fails on any `[...]` content
2. NHX vs BEAST comment clash: both use `[...]` but different internal syntax (`:` vs `,` separators, `&&NHX` vs `&` prefix). Forester explicitly checks for `&&NHX` first
3. Node vs branch ambiguity: no standard for which annotations belong to node vs edge. BEAST2 grammar is canonical (`nodeMeta` before `:`, `lengthMeta` after `:` before length). MrBayes puts branch attrs after length. Many tools do not distinguish
4. Comma-in-value escaping: BEAST uses `{v1,v2}`, IQ-TREE uses `"v1,v2"`. No standard mechanism
5. `[&R]`/`[&U]` recognition: part of Rich Newick, used by MrBayes, IQ-TREE, BEAST. Not widely recognized outside these tools
6. eNewick `#` in names: backward-compatible (legacy parsers treat as name character) but breaks if a node name legitimately contains `#`
7. Rich Newick extra colons: `name:length:bootstrap:probability` breaks standard parsers that try to parse additional colons as branch lengths

## Implications for v1 parser and writer

### Reader requirements

A superset parser that handles all dialects via content-based auto-detection:

- `[&&NHX` prefix -> NHX mode (colon-separated tags)
- `[&` prefix -> BEAST mode (comma-separated key=value)
- `#` in label -> eNewick hybrid markers
- `[&R]`/`[&U]`/`[&W]` before tree -> Rich Newick rooting/weight markers
- Plain `[...]` without `&` prefix -> standard comment (strip and ignore)

### Writer requirements

Configurable output dialect per `--nwk-style` flag:

- `plain`: no annotations, no comments. Maximum compatibility. Round-trips through every tool
- `annotated`: BEAST-style `[&key=value,...]`. Compatible with BEAST, FigTree, MrBayes consensus, IQ-TREE, DendroPy, BioPython
- `nhx`: NHX `[&&NHX:key=value:...]`. Compatible with Forester, ETE, DendroPy

Default `plain` because: the `bio` crate (current v1 reader) cannot parse annotations; annotated Newick from `optimize` output causes parse failure in `ancestral` (the bug that triggered this research); `.nexus` output provides no additional annotation capability.

The `--nwk-style` flag should apply to both `.nwk` and `.nexus` output because both embed the same Newick string. The `.nexus` container adds TAXA/TREES blocks but no annotation capability.

eNewick and Rich Newick are structural (DAG topology, extra colon fields), not just comment conventions. These belong on a separate `--output-tree-format` axis when network algorithms ship.

## Glossary

1. <a id="gloss-1"></a> **NEXUS.** Block-structured container file format for systematic data (<a id="cite-2b"></a>[Maddison, Swofford, and Maddison 1997](https://doi.org/10.1093/sysbio/46.4.590) [[2](#ref-2)]). Contains TAXA, CHARACTERS, TREES, and other blocks. Trees inside NEXUS are Newick strings. Comments `[...]` are formally discardable. [↩](#gloss-use-1)

2. <a id="gloss-2"></a> **NeXML.** XML-based replacement for NEXUS with proper annotation support via XML Schema (<a id="cite-3c"></a>[Vos et al. 2012](https://doi.org/10.1093/sysbio/sys025) [[3](#ref-3)]). Created by the NESCent working group after concluding NEXUS cannot be revised without breaking backward compatibility. [↩](#gloss-use-2)

3. <a id="gloss-3"></a> **Newick format.** Parenthetical notation for phylogenetic trees encoding topology, branch lengths, and node labels in a single string. Adopted 1986, named after Newick's restaurant in Dover, New Hampshire (<a id="cite-1c"></a>[Felsenstein 1986](https://phylipweb.github.io/phylip/newicktree.html) [[1](#ref-1)]). No formal specification document. [↩](#gloss-use-3)

## References

1. <a id="ref-1"></a> Felsenstein, Joseph. 1986. "The Newick Tree Format." PHYLIP documentation. https://phylipweb.github.io/phylip/newicktree.html [↩¹](#cite-1) [↩²](#cite-1b) [↩³](#cite-1c)

2. <a id="ref-2"></a> Maddison, David R., David L. Swofford, and Wayne P. Maddison. 1997. "Nexus: An Extensible File Format for Systematic Information." _Systematic Biology_ 46(4):590-621. https://doi.org/10.1093/sysbio/46.4.590 [↩¹](#cite-2a) [↩²](#cite-2b)

3. <a id="ref-3"></a> Vos, Rutger A., James P. Balhoff, Jason A. Caravas, et al. 2012. "NeXML: Rich, Extensible, and Verifiable Representation of Comparative Data and Metadata." _Systematic Biology_ 61(4):675-689. https://doi.org/10.1093/sysbio/sys025 [↩¹](#cite-3a) [↩²](#cite-3b) [↩³](#cite-3c)

4. <a id="ref-4"></a> Suchard, Marc A., Philippe Lemey, Guy Baele, Daniel L. Ayres, Alexei J. Drummond, and Andrew Rambaut. 2018. "Bayesian Phylogenetic and Phylodynamic Data Integration Using BEAST 1.10." _Virus Evolution_ 4(1):vey016. https://doi.org/10.1093/ve/vey016 [↩](#cite-4)

5. <a id="ref-5"></a> Bouckaert, Remco, Joseph Heled, Denise Kuhnert, Tim Vaughan, Chieh-Hsi Wu, Dong Xie, Marc A. Suchard, Andrew Rambaut, and Alexei J. Drummond. 2014. "BEAST 2: A Software Platform for Bayesian Evolutionary Analysis." _PLoS Computational Biology_ 10(4):e1003537. https://doi.org/10.1371/journal.pcbi.1003537 [↩](#cite-5)

6. <a id="ref-6"></a> Ronquist, Fredrik, and John P. Huelsenbeck. 2003. "MrBayes 3: Bayesian Phylogenetic Inference Under Mixed Models." _Bioinformatics_ 19(12):1572-1574. https://doi.org/10.1093/bioinformatics/btg180 [↩](#cite-6)

7. <a id="ref-7"></a> Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh. 2015. "IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies." _Molecular Biology and Evolution_ 32(1):268-274. https://doi.org/10.1093/molbev/msu300 [↩¹](#cite-7a) [↩²](#cite-7b)

8. <a id="ref-8"></a> Zmasek, Christian M., and Sean R. Eddy. 2001. "ATV: Display and Manipulation of Annotated Phylogenetic Trees." _Bioinformatics_ 17(4):383-384. https://doi.org/10.1093/bioinformatics/17.4.383 [↩](#cite-8)

9. <a id="ref-9"></a> Cardona, Gabriel, Francesc Rossello, and Gabriel Valiente. 2008. "Extended Newick: It Is Time for a Standard Representation of Phylogenetic Networks." _BMC Bioinformatics_ 9:532. https://doi.org/10.1186/1471-2105-9-532 [↩](#cite-9)
