# Multi-format tree I/O

v0 reads trees from Newick and Nexus (via `Bio.Phylo.read()` in [packages/legacy/treetime/treetime/treeanc.py#L328-L332](../../packages/legacy/treetime/treetime/treeanc.py#L328-L332)) and writes Newick, Nexus, and a basic Auspice JSON ([packages/legacy/treetime/treetime/CLI_io.py#L211-L226](../../packages/legacy/treetime/treetime/CLI_io.py#L211-L226)). v1 reads Newick and writes Newick in three comment styles, Nexus, Auspice v2 JSON, UShER MAT (protobuf and JSON), a graph JSON and Graphviz DOT. Every writer reads one struct of per-node facts. Parsers for Nexus, UShER MAT and PhyloXML exist as library crates, and no command reads these formats yet.

The motivation is interoperability with the phylogenetics tool landscape. TreeTime sits at the center of multiple workflows: it receives trees from phylogenetic inference tools (IQ-TREE, RAxML, FastTree) and produces timed trees consumed by downstream visualization and surveillance systems. Each of these systems has its own format:

- Pandemic surveillance (UShER, PANGOLIN, UCSC) operates on mutation-annotated trees in protobuf, storing millions of SARS-CoV-2 sequences. Writing MAT files lets UShER tools such as matUtils load TreeTime trees. Reading MAT files would let TreeTime date trees produced by UShER's parsimony placement.
- Genomic epidemiology (Nextstrain, Auspice) uses a JSON format embedding visualization metadata alongside the tree. TreeTime is already the engine behind `augur refine`; v1 writes Auspice JSON and augur node data JSON directly, so its results load into Auspice and into later augur steps without conversion.
- Comparative genomics (Archaeopteryx, Forester, ETE) uses PhyloXML for richly annotated trees carrying taxonomy, sequences, evolutionary events, geographic distributions, and protein domain architecture. Its structured annotation model carries data that Newick comment extensions cannot represent.
- Bayesian phylogenetics (BEAST, MrBayes, FigTree) uses Nexus and Newick with tool-specific comment conventions. These remain the baseline formats.

v0's format support was sufficient when TreeTime operated within the Nextstrain pipeline alone. v1 targets a broader set of workflows where trees arrive from and depart to different tools. [kb/proposals/io-format-coverage.md](../proposals/io-format-coverage.md) ranks the formats that v1 does not read or write yet.

## The phylogenetic format landscape

Phylogenetic tree formats have accumulated over four decades, each generation adding data dimensions that the previous could not represent.

**Newick** (Felsenstein, 1986 [[1]](#ref-1)) encodes topology and branch lengths in a single line of nested parentheses. Universal tool support, but no metadata mechanism, no schema, and ambiguous node/branch attribute handling. Ad-hoc extensions (NHX, BEAST `[&...]`, Rich Newick, eNewick) reuse the comment syntax `[...]` with mutually incompatible conventions.

**Nexus** (Maddison et al., 1997 [[2]](#ref-2)) wraps Newick trees in a block-structured container with taxa, character matrices, and distance blocks. Trees inside Nexus are still Newick strings, inheriting all Newick limitations.

**PhyloXML** (Han and Zmasek, 2009 [[3]](#ref-3)) and **NeXML** (Vos et al., 2012 [[4]](#ref-4)) address the metadata problem through XML schemas. PhyloXML focuses on tree annotation; NeXML modernizes the full Nexus scope with XSD validation. Neither has displaced Newick/Nexus due to ecosystem inertia.

**Auspice v2 JSON** (Hadfield et al., 2018 [[5]](#ref-5)) is a domain-specific format for real-time pathogen evolution visualization, embedding colorings, geographic coordinates, and genome annotations alongside the tree.

**Usher MAT protobuf** (Turakhia et al., 2021 [[6]](#ref-6)) represents mutation-annotated trees in a compact binary format designed for pandemic-scale datasets with millions of sequences.

The result is a fragmented landscape where most tools support Newick, many support Nexus, and specialized tools support one or two additional formats. Conversion between formats is common and lossy.

## Newick

Newick is a plain-text parenthetical format where tree topology is encoded by nesting and branch lengths follow a colon: `(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;` (Felsenstein, 1986 [[1]](#ref-1)). The grammar is minimal: `Tree -> Subtree ";"`, with subtrees being either leaves (name) or internal nodes (parenthesized children + name), each optionally carrying `:length`. Unrooted trees use an arbitrary root with 3 children; rooted binary trees have 2 children per internal node. No rooting marker exists in the base format (some extensions add `[&R]`/`[&U]`).

Underscores map to spaces in unquoted labels. Labels containing `()[],;:` or underscores must be single-quoted. The format has no mechanism for metadata beyond topology, labels, and branch lengths. This spawned a family of incompatible comment-based extensions:

- NHX (New Hampshire eXtended): `[&&NHX:key=value:key=value:...]` syntax for gene tree/species tree reconciliation metadata. Parsers unaware of NHX treat these as comments and ignore them.
- BEAST/MrBayes annotations: `[&key=value,key=value]` (single `&`). Attributes before the colon apply to the node, after the colon to the branch. Incompatible with NHX syntax.
- Extended Newick (eNewick): encodes phylogenetic networks (hybridization, lateral gene transfer) via `#` notation for reticulate nodes.
- Rich Newick: adds `[&U]`/`[&R]` prefix for unrooted/rooted designation, plus optional bootstrap and probability fields.

These extensions all reuse the comment mechanism `[...]` with different internal conventions. There is no way to distinguish node from branch attributes, no escaping for `]` or `,` inside annotation values, and no formal schema.

### v1 implementation

The `util-newick` crate reads Newick in two steps (`fn newick_from_reader()` in [packages/util-newick/src/parse.rs](../../packages/util-newick/src/parse.rs)). A flat `pest` grammar splits the input into tokens (`(`, `)`, `,`, `:` with a branch length, label, comment), and a loop builds the graph with an explicit stack of open parentheses. Neither step calls itself once per nesting level, so the depth of a tree is limited by memory, not by the thread stack: a ladder tree 100,000 levels deep reads on a 2 MiB stack. Nested `[` inside a comment are counted on the stack of the `pest` grammar (`PUSH`/`DROP`), which is also not recursive.

The reader accepts every input that v0 (Biopython) accepts in the cases below, and errors name the line and column:

- A missing final `;`, branch lengths such as `.5` and `+1`, a byte order mark, and comments between any two tokens, including before the first `(` and after the rooting comment `[&R]`. Comments after the final `;` are ignored; any other text after it is an error, as in v0
- A comment before `:` belongs to the node, and a comment after `:` belongs to the branch, also for unnamed nodes. `A:[&rate=1.5]` is a branch annotation without a branch length
- Inside comments, text in double quotes may contain `[`, `]` and `;`, so BEAST values such as `[&a="x]y"]` read correctly (v0 ends the comment at the first `]`). A single `"` in a free-text comment pairs with the next `"` in the file
- A numeric label on an internal node is read as a branch support value, not as a name, also when it is quoted (as in v0)
- A label is a plain name, also when it ends in `#` and digits, as in v0. eNewick hybrid nodes are read only when `NewickReadOptions::enewick` is set; the commands do not set it, so `EPI_ISL#402124` stays one sample name. With the option, a quoted label is never a hybrid node, and a quoted name takes its hybrid tag after the closing quote (`'x y'#H1`). Two occurrences of a hybrid node with different labels or with children in both, two edges between the same nodes, and cycles are errors

BEAST and NHX comments are kept as node and branch attributes (`fn classify_comment()` in [packages/util-newick/src/annotation.rs](../../packages/util-newick/src/annotation.rs)). Values in single or double quotes may contain commas. An unquoted number is kept as a `Number` when its text is the shortest form of its value and as `NumberText` with its original text otherwise (`0123`, `2020.50`, `1e5`), so writing it again gives the same text.

`fn tree_read()` in [packages/treetime-io/src/tree.rs](../../packages/treetime-io/src/tree.rs) reads the tree input of every command: a file that starts with `#NEXUS` is read as Nexus, any other file as Newick. `fn nwk_read()` in [packages/treetime-io/src/nwk.rs](../../packages/treetime-io/src/nwk.rs) builds the graph and the node names from the parse.

`fn nwk_write()` in the same file writes a checked rooted tree (`TreeView`), the node names, the branch lengths and an ordered list of `(key, value)` comments per node in one of three styles, selected with `--output-nwk-style`: `plain` (no comments, the default), `beast` (`[&key=value]`, strings in double quotes) and `nhx` (`[&&NHX:key=value]`). `fn nwk_node_comments()` in [packages/app-output/src/nwk_comments.rs](../../packages/app-output/src/nwk_comments.rs) builds the comments from the facts that the run has, in v0 order: branch mutations, then the date with two decimals, then the trait value. With more than one style, the files get the secondary extensions `.annotated` (beast) and `.nhx`.

Newick is the tree input format of every command (ancestral, clock, timetree, optimize, prune, mugration, homoplasy) and a tree output of every command.

## Nexus

Nexus (Maddison et al., 1997 [[2]](#ref-2)) is a block-structured container format with a `#NEXUS` header and `BEGIN...END` blocks. Standard blocks include `TAXA` (taxon count and labels), `CHARACTERS` (sequence alignments with format directives like `datatype=dna`, `missing=?`, `gap=-`), `TREES` (named Newick strings), `DISTANCES` (distance matrices), and `ASSUMPTIONS` (weight sets, character sets). Keywords are case-insensitive. Comments use `[...]`. Each predefined block type may appear once.

PAUP\*, MrBayes, Mesquite, MacClade, SplitsTree, BEAST, IQ-TREE all support Nexus. Trees inside Nexus are Newick strings, inheriting all Newick limitations. The format has no formal schema - validation depends on parser conventions and tools interpret edge cases differently.

### v1 implementation

`fn nex_write()` in [packages/treetime-io/src/nex.rs](../../packages/treetime-io/src/nex.rs) wraps the Newick text, in the same comment style, in a `#NEXUS` header with a `Taxa` block and a `Trees` block. `TaxLabels` lists the named leaves, and `NTax` counts the labels written. Every command writes `.nexus` output alongside `.nwk`.

Every command reads a Nexus tree file through `fn tree_read()` in [packages/treetime-io/src/tree.rs](../../packages/treetime-io/src/tree.rs). The file must contain exactly one tree, as `Bio.Phylo.read()` requires in v0; a file with several trees, for example a BEAST posterior sample, is an error that names the number of trees. v0 tries Newick first and reads Nexus only when Newick fails and the file extension is `.nexus` or `.nex`; v1 selects the reader by the `#NEXUS` header instead, so a Nexus file is read whatever its extension.

The Nexus parser (`fn nexus_from_string()` in [packages/util-newick/src/nexus.rs](../../packages/util-newick/src/nexus.rs)) is a `pest` grammar over blocks and commands that skips comments and quoted text when it looks for keywords and `;`. It reads `Translate` tables per `Trees` block, `Tree` and `UTree` commands with an optional `*`, comments between the tree name and `=` (BEAST 1 writes `tree STATE_0 [&lnP=...] = ...`), quoted tree names, and `End;` or `EndBlock;`. Other blocks and commands are skipped. A block without its end, or a command without its `;`, is an error.

## NeXML

NeXML (Vos et al., 2012 [[4]](#ref-4)) is an XML-based modernization of Nexus covering the full Nexus scope (taxa, character matrices, trees) with XSD schema validation and semantic web compatibility via RDFa-compatible ontology annotations. It is a sibling effort to PhyloXML: NeXML has broader scope (character matrices as first-class citizens), while PhyloXML is more tree-annotation-focused. Supported by DendroPy, Bio::Phylo (Perl), and the NeXML Java API.

v1 does not implement NeXML. It is mentioned here for context as the other XML-based phylogenetic format that readers will encounter when evaluating PhyloXML.

## PhyloXML

### What the format represents

PhyloXML is an XML language for phylogenetic trees and associated data (Han and Zmasek, 2009 [[3]](#ref-3)), defined by an XSD schema (version 1.20, namespace `http://www.phyloxml.org`). A PhyloXML document contains one or more `<phylogeny>` elements, each holding a recursive `<clade>` tree. Each clade can carry typed, structured annotations:

- Taxonomy: scientific name, common name, rank (domain through subspecies), authority, taxonomic identifier, synonyms
- Sequences: molecular sequences (DNA, RNA, protein), accessions, annotations, domain architecture with positions and confidence values
- Events: speciation, duplication, loss, lateral gene transfer, fusion - enabling gene tree/species tree reconciliation
- Distributions: geographic coordinates (latitude/longitude), polygons, geodetic datums
- Dates: temporal data with value, minimum, maximum, and units (e.g. "mya")
- Confidence: bootstrap support, posterior probability, or custom confidence types
- Binary characters: gained/lost/present/absent character tracking
- Properties: arbitrary typed key-value pairs using namespace-prefixed references (e.g. `NOAA:depth`), applicable to phylogenies, clades, nodes, branches, or annotations
- Display: branch color (RGB, cascading to descendants) and width
- Cross-references: clade-to-clade and sequence-to-sequence relations (orthology, paralogy, xenology)

The `<property>` element is the extensibility mechanism: it accepts any XSD datatype and any namespace-prefixed key, allowing custom annotations without schema modification. The `id_source`/`id_ref` mechanism enables internal cross-linking within a document.

### Tool support

PhyloXML is the native format for Archaeopteryx (Java tree viewer for phylogenomics, successor to ATV) and the Forester library. Biopython's `Bio.Phylo.PhyloXML` module provides full read/write. ETE toolkit, Dendroscope, and MEGA also support it.

### Tradeoffs

XML syntax is inherently verbose - a tree that takes one line in Newick requires many lines in PhyloXML. File sizes are larger, and parsing requires an XML parser rather than the simple recursive-descent parsers used for Newick. For workflows that need only topology and branch lengths, PhyloXML is overkill. Its value is in richly annotated trees where the typed annotation structure prevents the ad-hoc metadata encoding that plagues Newick comment extensions.

Unlike NeXML and Nexus, PhyloXML does not support character matrices or alignments as standalone blocks. Sequences attach to individual clades, not to a global alignment structure.

### v1 implementation

The `util-phyloxml` crate ([packages/util-phyloxml/src/types.rs](../../packages/util-phyloxml/src/types.rs)) implements the full PhyloXML type model: `Phyloxml`, `PhyloxmlPhylogeny`, `PhyloxmlClade`, `PhyloxmlTaxonomy`, `PhyloxmlSequence`, `PhyloxmlEvents`, `PhyloxmlDistribution`, `PhyloxmlDate`, `PhyloxmlProperty`, `PhyloxmlBinaryCharacters`, `PhyloxmlDomainArchitecture`, and supporting types. `fn phyloxml_read()` and `fn phyloxml_write()` in [packages/util-phyloxml/src/lib.rs](../../packages/util-phyloxml/src/lib.rs) read and write XML through `quick-xml` with serde.

No package depends on the crate, and no converter between PhyloXML and the graph exists. [kb/issues/N-io-phyloxml-crate-unused.md](../issues/N-io-phyloxml-crate-unused.md) asks whether to wire it into the commands or remove it.

## Usher MAT protobuf

### The mutation-annotated tree

UShER (Ultrafast Sample placement on Existing tRees) was developed at UC Santa Cruz by Turakhia, Thornlow, Hinrichs, De Maio, and colleagues for real-time phylogenetic placement during the SARS-CoV-2 pandemic (Turakhia et al., 2021 [[6]](#ref-6)). The companion matUtils toolkit provides utilities for manipulating MAT files (McBroome et al., 2021 [[7]](#ref-7)).

The core data structure is the **mutation-annotated tree (MAT)**: a phylogenetic tree where each branch carries the mutations inferred by maximum parsimony (Fitch's algorithm). Each mutation is recorded once on the ancestral branch where it occurred, rather than redundantly for every descendant. This "evolutionary compression" yields dramatic storage savings: a tree of 400,000 SARS-CoV-2 genomes requires 31 MB as MAT vs. 12 GB as MSA or 21 GB as VCF. Full sample genotypes can be reconstructed by traversing the root-to-tip mutation path for any leaf.

### Protobuf schema

The MAT is serialized using Google Protocol Buffers. The `parsimony.proto` schema has four components:

- `data` (top-level): a Newick string, per-node mutation lists in preorder traversal order, condensed nodes (groups of identical sequences collapsed into single representatives), and per-node metadata
- `mut`: individual mutation with `position`, `ref_nuc`, `par_nuc` (parent nucleotide), repeated `mut_nuc` (child nucleotides), and `chromosome`. Nucleotide encoding: 0=A, 1=C, 2=G, 3=T
- `condensed_node`: maps a node name to a list of identical leaf names, compacting polytomies of identical sequences
- `node_metadata`: per-node clade annotation strings

A second schema (`mutation_detailed.proto`) provides a more compact binary representation with packed mutation fields, plus structures for sample placement workflows.

### Ecosystem context

UShER operates within the pandemic genomics surveillance ecosystem:

- GISAID: primary global repository of SARS-CoV-2 sequences (over 14 million by 2023). UShER processes sequences deposited here.
- Nextstrain: real-time pathogen evolution visualization platform (Hadfield et al., 2018 [[5]](#ref-5)). The `augur refine` step calls TreeTime (Sagulenko et al., 2018 [[8]](#ref-8)) for molecular clock inference.
- PANGOLIN: SARS-CoV-2 lineage assignment tool. UShER serves as an alternative backend for lineage placement, offering faster and more accurate results than the machine learning model for closely related sequences.
- Taxonium: web-based viewer for large mutation-annotated trees, consuming the related `taxodium.proto` schema.
- UCSC Genome Browser: hosts the UShER web interface and maintains daily updated global SARS-CoV-2 phylogenetic trees.

The format's design reflects pandemic-scale constraints: millions of sequences arriving daily, trees updated continuously, and parsimony-based placement as the only tractable approach at this scale. For closely related sequences (SARS-CoV-2 genomes differ by only a handful of mutations), parsimony and likelihood methods converge in accuracy.

### Tradeoffs

The MAT format is compact and fast for its intended use case (mutation-annotated trees at pandemic scale). The binary protobuf encoding requires the protobuf runtime for reading and writing, and the format stores only mutations, not full sequences or probability distributions. Metadata beyond clade annotations requires the extended schemas. The format is specific to the UShER ecosystem and not a community standard like Newick or PhyloXML.

### v1 implementation

The `util-usher-mat` crate ([packages/util-usher-mat/src/lib.rs](../../packages/util-usher-mat/src/lib.rs)) compiles the `parsimony.proto`, `mutation_detailed.proto` and `taxodium.proto` schemas with `prost` and provides `fn usher_mat_pb_read()` and `fn usher_mat_pb_write()`.

`fn mat_tree()` in [packages/app-output/src/usher_mat.rs](../../packages/app-output/src/usher_mat.rs) builds the `UsherTree` from the checked tree: the Newick string and the per-node mutation lists in preorder. UShER has no gap state, so a deletion is written as missing data and an insertion is left out, with one warning per run ([kb/decisions/io-usher-mat-gaps-as-missing-data.md](io-usher-mat-gaps-as-missing-data.md)). The run writes the tree as protobuf (`.mat.pb`, through `fn usher_mat_pb_write_file()` in [packages/treetime-io/src/usher_mat.rs](../../packages/treetime-io/src/usher_mat.rs)) and as JSON of the same type model (`.mat.json`), which can be read without protobuf tooling.

No command reads a MAT: the parser exists, but no converter from a MAT to the graph and its partitions does ([kb/issues/N-io-usher-mat-input-not-wired.md](../issues/N-io-usher-mat-input-not-wired.md)).

## Auspice v2 JSON

### v0 output

v0 writes a basic Auspice JSON alongside Nexus output via `create_auspice_json()` ([packages/legacy/treetime/treetime/CLI_io.py#L277-L341](../../packages/legacy/treetime/treetime/CLI_io.py#L277-L341)). This output contains tree structure, mutations, dates, and divergence, and targets Nextstrain's Auspice viewer for pathogen evolution visualization. v0 cannot read Auspice JSON.

### The Auspice v2 schema

Nextstrain was co-developed by Trevor Bedford (Fred Hutchinson Cancer Research Center) and Richard Neher (University of Basel) for real-time tracking of pathogen evolution (Hadfield et al., 2018 [[5]](#ref-5)). The platform won the inaugural Open Science Prize in 2017. The data pipeline flows from raw sequences through Augur's (Huddleston et al., 2021 [[9]](#ref-9)) modular subcommands (`filter`, `align`, `tree`, `refine`, `ancestral`, `traits`, `export`) to the Auspice v2 JSON consumed by the viewer. TreeTime (Sagulenko et al., 2018 [[8]](#ref-8)) is the computational engine called by `augur refine` for molecular clock inference and ancestral reconstruction.

The Auspice v2 format (schema ID `https://nextstrain.org/schemas/dataset/v2`) combines metadata and tree in a single JSON file with four top-level fields:

- `version`: constant `"v2"`
- `meta`: dataset metadata including `updated` date, `panels` (tree, map, frequencies, entropy, measurements), `colorings` (typed color-by options with scales and legends), `geo_resolutions` (geographic trait definitions with latitude/longitude per deme), `genome_annotations` (nucleotide and CDS positions in 1-based GFF convention), `display_defaults` (layout, distance measure, default color-by), `filters`, `maintainers`, and `data_provenance`
- `tree`: recursive node structure where each node has `name`, `node_attrs` (divergence `div`, decimal date `num_date` with confidence intervals, arbitrary traits as `{"value": ..., "confidence": ...}` objects), `branch_attrs` (per-gene mutation lists, branch labels), and optional `children`
- `root_sequence`: optional mapping of gene/segment names to reference sequences

The v2 format replaced the earlier v1 format which split data across two files (`meta.json` + `tree.json`), used 0-based genome annotation positions, and lacked structured `node_attrs`/`branch_attrs`.

### Ecosystem role

Nextstrain/Auspice is the standard visualization platform in genomic epidemiology, powering nextstrain.org for continuous tracking of influenza, SARS-CoV-2, dengue, Ebola, Zika, measles, monkeypox, RSV, tuberculosis, and more. The platform is referenced by WHO, CDC, and national public health agencies.

### Tradeoffs

The format embeds rich visualization metadata directly in the data file, making it web-ready for the Auspice viewer. JSON overhead and cumulative per-node divergence storage produce large files for big trees. The format is Nextstrain-specific - BEAST, IQ-TREE, FigTree, and other phylogenetics tools do not produce or consume it. Dynamic trait properties (`node_attrs` pattern properties) make static typing difficult, requiring serde flatten or similar mechanisms.

### v1 implementation

The type model in [packages/treetime-io/src/auspice_types.rs](../../packages/treetime-io/src/auspice_types.rs) covers the Auspice v2 schema: `AuspiceTree`, `AuspiceTreeNode`, `AuspiceTreeNodeAttrs`, `AuspiceTreeBranchAttrs`, `AuspiceTreeMeta`, `AuspiceGenomeAnnotations`, colorings, and display defaults. The types derive serde in both directions, but no command reads Auspice JSON.

`fn auspice_tree()` in [packages/app-output/src/auspice.rs](../../packages/app-output/src/auspice.rs) builds the tree from the checked tree and the facts that the run has:

- `div`: summed from the root down along the divergence branch lengths, in the order of augur and v0, or taken from the values the command computed
- Dates: `num_date` with its confidence interval, plus the `Date` and `Excluded` (`bad_branch`) colorings, when the run inferred dates
- Trait: the value, the confidence of each state above 0.001 and the entropy, plus a coloring and a filter for the attribute, when the run reconstructed a discrete trait
- Mutations: nucleotide substitutions and amino-acid mutations per CDS in `branch_attrs`, the `Genotype` (`gt`) coloring when at least one branch shows a mutation ([kb/decisions/auspice-genotype-coloring-when-mutations-shown.md](auspice-genotype-coloring-when-mutations-shown.md)), the root sequences, and the genome annotations with the entropy panel

Every number in the node attributes must be finite; otherwise the error names the node and the field.

## Supporting formats

### Graph JSON

The graph type derives serde, and every command writes the graph (nodes, edges and graph-level data) as `{stem}.graph.json` through `json_write_file()` of `treetime-utils`, with no format-specific code. The file shows the graph as the run left it, also when it is not a tree. No command reads it.

### Graphviz DOT

`fn graphviz_write_file()` in [packages/treetime-io/src/graphviz.rs](../../packages/treetime-io/src/graphviz.rs) writes a directed graph with subgraphs for the roots, the internal nodes and the leaves, and with an invisible edge for every ordered pair of leaves (and of roots), which keeps them in one row. Each root appears once, and labels are escaped. Every command writes it as `{stem}.dot`.

## Output architecture

Each command builds one `AnnotatedGraph` ([packages/app-output/src/annotated_graph.rs](../../packages/app-output/src/annotated_graph.rs)) after its last change to the graph. The struct borrows the graph, the node names, the divergence branch lengths, the time branch lengths (timetree only) and the divergence values. It has three optional fact groups: `sequences` (root sequence, branch mutations, amino acids), `dates` and `traits`. A command sets a group only when the run produced those facts, and a writer decides what to write from the facts that are present, never from the command.

`fn write_graph_outputs()` in [packages/app-output/src/tree_output.rs](../../packages/app-output/src/tree_output.rs) writes the graph JSON and DOT files. `AnnotatedTreeView::new()` then checks that the graph is one rooted tree with `TreeView` ([packages/treetime-graph/src/tree_view.rs](../../packages/treetime-graph/src/tree_view.rs)): exactly one root, at most one parent per node, no cycle. Newick, Nexus, Auspice JSON, UShER MAT, augur node data JSON and the per-node tables take only the checked view, so they cannot run on a graph that is not a tree. When the check fails, the graph files are already written, and the error lists the tree outputs that were not written.

Output paths come from one output plan ([packages/app-output/src/output_plan.rs](../../packages/app-output/src/output_plan.rs)): `--output-all <dir>` writes the default set of the command, `--output-selection` restricts that set, and a per-file flag such as `--output-nwk` sets one path. The plan rejects two outputs with the same destination. An output that `--output-all` selects but the run has no data for is skipped with a debug message; a per-file flag for such an output fails.

The tree readers and writers open files through `read_file_with()` and `write_file_with()` of `treetime-utils`, so gz, bz2, xz and zst compression follows the file extension on input and output, and `-` means standard input or output.

### Current limitations

- The commands read trees from Newick only. UShER MAT, Auspice JSON, Nexus and PhyloXML input need conversion to Newick (and FASTA for sequences) first
- The Newick reader rejects eNewick hybrid nodes, although the graph can hold networks
- Graphviz DOT output grows with the square of the number of leaves, because of the invisible edges ([kb/issues/M-io-output-layer-scaling-naming-and-duplication.md](../issues/M-io-output-layer-scaling-naming-and-duplication.md))

## Practical impact

- TreeTime results load into Auspice and into later augur steps without conversion
- UShER tools such as matUtils can load TreeTime trees with their mutations
- BEAST-style tools (FigTree) read the dates, mutations and traits from the annotated Newick and Nexus files
- Trees from UShER, Auspice datasets and Nexus files need conversion to Newick before TreeTime can analyze them

## References

<a id="ref-1"></a>[1] Felsenstein J. "The Newick tree format." 1986. <http://evolution.genetics.washington.edu/phylip/newicktree.html>

<a id="ref-2"></a>[2] Maddison DR, Swofford DL, Maddison WP. "NEXUS: an extensible file format for systematic information." Systematic Biology 46(4):590-621, 1997. <https://doi.org/10.1093/sysbio/46.4.590>

<a id="ref-3"></a>[3] Han MV, Zmasek CM. "phyloXML: XML for evolutionary biology and comparative genomics." BMC Bioinformatics 10:356, 2009. <https://doi.org/10.1186/1471-2105-10-356>

<a id="ref-4"></a>[4] Vos RA, Balhoff JP, Caravas JA, et al. "NeXML: rich, extensible, and verifiable representation of comparative data and metadata." Systematic Biology 61(4):675-689, 2012. <https://doi.org/10.1093/sysbio/sys025>

<a id="ref-5"></a>[5] Hadfield J, Megill C, Bell SM, Huddleston J, Potter B, Callender C, Sagulenko P, Bedford T, Neher RA. "Nextstrain: real-time tracking of pathogen evolution." Bioinformatics 34(23):4121-4123, 2018. <https://doi.org/10.1093/bioinformatics/bty407>

<a id="ref-6"></a>[6] Turakhia Y, Thornlow B, Hinrichs AS, De Maio N, Gozashti L, Lanfear R, Haussler D, Corbett-Detig R. "Ultrafast Sample placement on Existing tRees (UShER) enables real-time phylogenetics for the SARS-CoV-2 pandemic." Nature Genetics 53(6):809-816, 2021. <https://doi.org/10.1038/s41588-021-00862-7>

<a id="ref-7"></a>[7] McBroome J, Thornlow B, Hinrichs AS, De Maio N, Goldman N, Haussler D, Corbett-Detig R, Turakhia Y. "A Daily-Updated Database and Tools for Comprehensive SARS-CoV-2 Mutation-Annotated Trees." Molecular Biology and Evolution 38(12):5819-5824, 2021. <https://doi.org/10.1093/molbev/msab264>

<a id="ref-8"></a>[8] Sagulenko P, Puller V, Neher RA. "TreeTime: Maximum-likelihood phylodynamic analysis." Virus Evolution 4(1):vex042, 2018. <https://doi.org/10.1093/ve/vex042>

<a id="ref-9"></a>[9] Huddleston J, Hadfield J, Sibley TR, et al. "Augur: a bioinformatics toolkit for phylogenetic analyses of human pathogens." Journal of Open Source Software 6(57):2906, 2021. <https://doi.org/10.21105/joss.02906>
