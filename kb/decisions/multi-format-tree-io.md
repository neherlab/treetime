# Multi-format tree I/O

v0 reads trees from Newick and Nexus (via `Bio.Phylo.read()` in [packages/legacy/treetime/treetime/treeanc.py#L328-L332](../../packages/legacy/treetime/treetime/treeanc.py#L328-L332)) and writes Newick, Nexus, and a basic Auspice JSON ([packages/legacy/treetime/treetime/CLI_io.py#L211-L226](../../packages/legacy/treetime/treetime/CLI_io.py#L211-L226)). v1 expands to eight formats with full read/write support and a common graph intermediate representation.

The motivation is interoperability with the phylogenetics tool landscape. TreeTime sits at the center of multiple workflows: it receives trees from phylogenetic inference tools (IQ-TREE, RAxML, FastTree) and produces timed trees consumed by downstream visualization and surveillance systems. Each of these systems has its own format:

- **Pandemic surveillance** (UShER, PANGOLIN, UCSC) operates on mutation-annotated trees in protobuf, storing millions of SARS-CoV-2 sequences. Reading MAT files lets TreeTime perform molecular clock inference on trees produced by UShER's parsimony placement.
- **Genomic epidemiology** (Nextstrain, Auspice) uses a JSON format embedding visualization metadata alongside the tree. TreeTime is already the engine behind `augur refine`; bidirectional Auspice JSON support allows reading back previously exported datasets for re-analysis.
- **Comparative genomics** (Archaeopteryx, Forester, ETE) uses PhyloXML for richly annotated trees carrying taxonomy, sequences, evolutionary events, geographic distributions, and protein domain architecture. PhyloXML support enables TreeTime output to carry structured annotations that Newick comment extensions cannot represent.
- **Bayesian phylogenetics** (BEAST, MrBayes, FigTree) uses Nexus and Newick with tool-specific comment conventions. These remain the baseline formats.

v0's format support was sufficient when TreeTime operated within the Nextstrain pipeline alone. v1 targets a broader set of workflows where trees arrive from and depart to different tools. The `convert` subcommand ([packages/treetime-cli/src/convert/convert.rs#L71-L97](../../packages/treetime-cli/src/convert/convert.rs#L71-L97)) is a byproduct: once the format adapters exist, exposing them as a standalone conversion tool costs nothing.

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

- **NHX** (New Hampshire eXtended): `[&&NHX:key=value:key=value:...]` syntax for gene tree/species tree reconciliation metadata. Parsers unaware of NHX treat these as comments and ignore them.
- **BEAST/MrBayes annotations**: `[&key=value,key=value]` (single `&`). Attributes before the colon apply to the node, after the colon to the branch. Incompatible with NHX syntax.
- **Extended Newick (eNewick)**: encodes phylogenetic networks (hybridization, lateral gene transfer) via `#` notation for reticulate nodes.
- **Rich Newick**: adds `[&U]`/`[&R]` prefix for unrooted/rooted designation, plus optional bootstrap and probability fields.

These extensions all reuse the comment mechanism `[...]` with different internal conventions. There is no way to distinguish node from branch attributes, no escaping for `]` or `,` inside annotation values, and no formal schema.

### v1 implementation

v1 reads Newick via the `bio` crate's parser and writes via a custom serializer ([packages/treetime-io/src/nwk.rs](../../packages/treetime-io/src/nwk.rs)). The `NodeFromNwk`/`NodeToNwk` and `EdgeFromNwk`/`EdgeToNwk` traits define how graph nodes and edges map to Newick labels, comments, and branch lengths. A `CommentProviders` mechanism allows injecting external annotations (e.g. ancestral state metadata) into Newick comments during output.

Newick is the primary tree input format for all v1 analysis commands (ancestral, clock, timetree, optimize, prune, mugration) and is produced as output by all of them. In the convert command, Newick serves as both input and output format with auto-detection from `.nwk`/`.newick` extensions.

## Nexus

Nexus (Maddison et al., 1997 [[2]](#ref-2)) is a block-structured container format with a `#NEXUS` header and `BEGIN...END` blocks. Standard blocks include `TAXA` (taxon count and labels), `CHARACTERS` (sequence alignments with format directives like `datatype=dna`, `missing=?`, `gap=-`), `TREES` (named Newick strings), `DISTANCES` (distance matrices), and `ASSUMPTIONS` (weight sets, character sets). Keywords are case-insensitive. Comments use `[...]`. Each predefined block type may appear once.

PAUP\*, MrBayes, Mesquite, MacClade, SplitsTree, BEAST, IQ-TREE all support Nexus. Trees inside Nexus are Newick strings, inheriting all Newick limitations. The format has no formal schema - validation depends on parser conventions and tools interpret edge cases differently.

### v1 implementation

v1 writes Nexus by wrapping the Newick writer with `#NEXUS`/`Begin Taxa`/`Begin Trees` envelope blocks ([packages/treetime-io/src/nex.rs](../../packages/treetime-io/src/nex.rs)). All analysis commands produce `.nexus` output alongside `.nwk`. Nexus reading is not implemented - the convert command panics with `unimplemented!()` for Nexus input. v0 supports Nexus reading as a fallback when Newick parsing fails and the file extension is `.nexus` or `.nex`.

## NeXML

NeXML (Vos et al., 2012 [[4]](#ref-4)) is an XML-based modernization of Nexus covering the full Nexus scope (taxa, character matrices, trees) with XSD schema validation and semantic web compatibility via RDFa-compatible ontology annotations. It is a sibling effort to PhyloXML: NeXML has broader scope (character matrices as first-class citizens), while PhyloXML is more tree-annotation-focused. Supported by DendroPy, Bio::Phylo (Perl), and the NeXML Java API.

v1 does not implement NeXML. It is mentioned here for context as the other XML-based phylogenetic format that readers will encounter when evaluating PhyloXML.

## PhyloXML

### What the format represents

PhyloXML is an XML language for phylogenetic trees and associated data (Han and Zmasek, 2009 [[3]](#ref-3)), defined by an XSD schema (version 1.20, namespace `http://www.phyloxml.org`). A PhyloXML document contains one or more `<phylogeny>` elements, each holding a recursive `<clade>` tree. Each clade can carry typed, structured annotations:

- **Taxonomy**: scientific name, common name, rank (domain through subspecies), authority, taxonomic identifier, synonyms
- **Sequences**: molecular sequences (DNA, RNA, protein), accessions, annotations, domain architecture with positions and confidence values
- **Events**: speciation, duplication, loss, lateral gene transfer, fusion - enabling gene tree/species tree reconciliation
- **Distributions**: geographic coordinates (latitude/longitude), polygons, geodetic datums
- **Dates**: temporal data with value, minimum, maximum, and units (e.g. "mya")
- **Confidence**: bootstrap support, posterior probability, or custom confidence types
- **Binary characters**: gained/lost/present/absent character tracking
- **Properties**: arbitrary typed key-value pairs using namespace-prefixed references (e.g. `NOAA:depth`), applicable to phylogenies, clades, nodes, branches, or annotations
- **Display**: branch color (RGB, cascading to descendants) and width
- **Cross-references**: clade-to-clade and sequence-to-sequence relations (orthology, paralogy, xenology)

The `<property>` element is the extensibility mechanism: it accepts any XSD datatype and any namespace-prefixed key, allowing custom annotations without schema modification. The `id_source`/`id_ref` mechanism enables internal cross-linking within a document.

### Tool support

PhyloXML is the native format for Archaeopteryx (Java tree viewer for phylogenomics, successor to ATV) and the Forester library. Biopython's `Bio.Phylo.PhyloXML` module provides full read/write. ETE toolkit, Dendroscope, and MEGA also support it.

### Tradeoffs

XML syntax is inherently verbose - a tree that takes one line in Newick requires many lines in PhyloXML. File sizes are larger, and parsing requires an XML parser rather than the simple recursive-descent parsers used for Newick. For workflows that need only topology and branch lengths, PhyloXML is overkill. Its value is in richly annotated trees where the typed annotation structure prevents the ad-hoc metadata encoding that plagues Newick comment extensions.

Unlike NeXML and Nexus, PhyloXML does not support character matrices or alignments as standalone blocks. Sequences attach to individual clades, not to a global alignment structure.

### v1 implementation

The `phyloxml` crate ([packages/phyloxml/src/types.rs](../../packages/phyloxml/src/types.rs)) implements the full PhyloXML type model: `Phyloxml`, `PhyloxmlPhylogeny`, `PhyloxmlClade`, `PhyloxmlTaxonomy`, `PhyloxmlSequence`, `PhyloxmlEvents`, `PhyloxmlDistribution`, `PhyloxmlDate`, `PhyloxmlProperty`, `PhyloxmlBinaryCharacters`, `PhyloxmlDomainArchitecture`, and supporting types. Reading uses `quick-xml` with serde deserialization. Writing uses the `xml-rs` crate for pretty-printed output.

The `treetime-io` crate ([packages/treetime-io/src/phyloxml.rs](../../packages/treetime-io/src/phyloxml.rs)) provides graph integration through `PhyloxmlToGraph` and `PhyloxmlFromGraph` traits, converting between the PhyloXML type model and the generic `Graph` structure.

A **PhyloXML-JSON** variant serializes the same type model as JSON instead of XML, providing a more compact and JavaScript-friendly representation. Both XML and JSON variants are available for input and output in the convert command.

## Usher MAT protobuf

### The mutation-annotated tree

UShER (Ultrafast Sample placement on Existing tRees) was developed at UC Santa Cruz by Turakhia, Thornlow, Hinrichs, De Maio, and colleagues for real-time phylogenetic placement during the SARS-CoV-2 pandemic (Turakhia et al., 2021 [[6]](#ref-6)). The companion matUtils toolkit provides utilities for manipulating MAT files (McBroome et al., 2021 [[7]](#ref-7)).

The core data structure is the **mutation-annotated tree (MAT)**: a phylogenetic tree where each branch carries the mutations inferred by maximum parsimony (Fitch's algorithm). Each mutation is recorded once on the ancestral branch where it occurred, rather than redundantly for every descendant. This "evolutionary compression" yields dramatic storage savings: a tree of 400,000 SARS-CoV-2 genomes requires 31 MB as MAT vs. 12 GB as MSA or 21 GB as VCF. Full sample genotypes can be reconstructed by traversing the root-to-tip mutation path for any leaf.

### Protobuf schema

The MAT is serialized using Google Protocol Buffers. The `parsimony.proto` schema has four components:

- **`data`** (top-level): a Newick string, per-node mutation lists in preorder traversal order, condensed nodes (groups of identical sequences collapsed into single representatives), and per-node metadata
- **`mut`**: individual mutation with `position`, `ref_nuc`, `par_nuc` (parent nucleotide), repeated `mut_nuc` (child nucleotides), and `chromosome`. Nucleotide encoding: 0=A, 1=C, 2=G, 3=T
- **`condensed_node`**: maps a node name to a list of identical leaf names, compacting polytomies of identical sequences
- **`node_metadata`**: per-node clade annotation strings

A second schema (`mutation_detailed.proto`) provides a more compact binary representation with packed mutation fields, plus structures for sample placement workflows.

### Ecosystem context

UShER operates within the pandemic genomics surveillance ecosystem:

- **GISAID**: primary global repository of SARS-CoV-2 sequences (over 14 million by 2023). UShER processes sequences deposited here.
- **Nextstrain**: real-time pathogen evolution visualization platform (Hadfield et al., 2018 [[5]](#ref-5)). The `augur refine` step calls TreeTime (Sagulenko et al., 2018 [[8]](#ref-8)) for molecular clock inference.
- **PANGOLIN**: SARS-CoV-2 lineage assignment tool. UShER serves as an alternative backend for lineage placement, offering faster and more accurate results than the machine learning model for closely related sequences.
- **Taxonium**: web-based viewer for large mutation-annotated trees, consuming the related `taxodium.proto` schema.
- **UCSC Genome Browser**: hosts the UShER web interface and maintains daily updated global SARS-CoV-2 phylogenetic trees.

The format's design reflects pandemic-scale constraints: millions of sequences arriving daily, trees updated continuously, and parsimony-based placement as the only tractable approach at this scale. For closely related sequences (SARS-CoV-2 genomes differ by only a handful of mutations), parsimony and likelihood methods converge in accuracy.

### Tradeoffs

The MAT format is compact and fast for its intended use case (mutation-annotated trees at pandemic scale). The binary protobuf encoding requires the protobuf runtime for reading and writing, and the format stores only mutations, not full sequences or probability distributions. Metadata beyond clade annotations requires the extended schemas. The format is specific to the UShER ecosystem and not a community standard like Newick or PhyloXML.

### v1 implementation

The `usher-mat-utils` crate ([packages/usher-mat-utils/src/lib.rs](../../packages/usher-mat-utils/src/lib.rs)) handles protobuf encoding and decoding via the `prost` crate. Three protobuf schemas are compiled: `parsimony.proto`, `mutation_detailed.proto`, and `taxodium.proto`.

The `treetime-io` crate ([packages/treetime-io/src/usher_mat.rs](../../packages/treetime-io/src/usher_mat.rs)) provides graph integration through `UsherRead` and `UsherWrite` traits. The `usher_to_graph()` function parses the embedded Newick string and attaches per-node mutations from the protobuf data. The `usher_from_graph()` function serializes the graph back to the `UsherTree` protobuf structure.

A **Usher MAT-JSON** variant serializes the same protobuf type model as JSON, useful for inspection and debugging of MAT files without protobuf tooling. Both protobuf and JSON variants are available for input and output.

## Auspice v2 JSON

### From output-only to bidirectional

v0 writes a basic Auspice JSON alongside Nexus output via `create_auspice_json()` ([packages/legacy/treetime/treetime/CLI_io.py#L277-L341](../../packages/legacy/treetime/treetime/CLI_io.py#L277-L341)). This output contains tree structure, mutations, dates, and divergence, but v0 cannot read Auspice JSON back. The output targets Nextstrain's Auspice viewer for pathogen evolution visualization.

v1 implements full bidirectional Auspice v2 JSON support: both reading and writing through a complete type model that covers the Auspice v2 schema.

### The Auspice v2 schema

Nextstrain was co-developed by Trevor Bedford (Fred Hutchinson Cancer Research Center) and Richard Neher (University of Basel) for real-time tracking of pathogen evolution (Hadfield et al., 2018 [[5]](#ref-5)). The platform won the inaugural Open Science Prize in 2017. The data pipeline flows from raw sequences through Augur's (Huddleston et al., 2021 [[9]](#ref-9)) modular subcommands (`filter`, `align`, `tree`, `refine`, `ancestral`, `traits`, `export`) to the Auspice v2 JSON consumed by the viewer. TreeTime (Sagulenko et al., 2018 [[8]](#ref-8)) is the computational engine called by `augur refine` for molecular clock inference and ancestral reconstruction.

The Auspice v2 format (schema ID `https://nextstrain.org/schemas/dataset/v2`) combines metadata and tree in a single JSON file with four top-level fields:

- **`version`**: constant `"v2"`
- **`meta`**: dataset metadata including `updated` date, `panels` (tree, map, frequencies, entropy, measurements), `colorings` (typed color-by options with scales and legends), `geo_resolutions` (geographic trait definitions with latitude/longitude per deme), `genome_annotations` (nucleotide and CDS positions in 1-based GFF convention), `display_defaults` (layout, distance measure, default color-by), `filters`, `maintainers`, and `data_provenance`
- **`tree`**: recursive node structure where each node has `name`, `node_attrs` (divergence `div`, decimal date `num_date` with confidence intervals, arbitrary traits as `{"value": ..., "confidence": ...}` objects), `branch_attrs` (per-gene mutation lists, branch labels), and optional `children`
- **`root_sequence`**: optional mapping of gene/segment names to reference sequences

The v2 format replaced the earlier v1 format which split data across two files (`meta.json` + `tree.json`), used 0-based genome annotation positions, and lacked structured `node_attrs`/`branch_attrs`.

### Ecosystem role

Nextstrain/Auspice is the standard visualization platform in genomic epidemiology, powering nextstrain.org for continuous tracking of influenza, SARS-CoV-2, dengue, Ebola, Zika, measles, monkeypox, RSV, tuberculosis, and more. The platform is referenced by WHO, CDC, and national public health agencies.

### Tradeoffs

The format embeds rich visualization metadata directly in the data file, making it web-ready for the Auspice viewer. JSON overhead and cumulative per-node divergence storage produce large files for big trees. The format is Nextstrain-specific - BEAST, IQ-TREE, FigTree, and other phylogenetics tools do not produce or consume it. Dynamic trait properties (`node_attrs` pattern properties) make static typing difficult, requiring serde flatten or similar mechanisms.

### v1 implementation

The type model in [packages/treetime-io/src/auspice_types.rs](../../packages/treetime-io/src/auspice_types.rs) covers the full Auspice v2 schema: `AuspiceTree`, `AuspiceTreeNode`, `AuspiceTreeNodeAttrs`, `AuspiceTreeBranchAttrs`, `AuspiceTreeMeta`, `AuspiceGenomeAnnotations`, colorings, geo resolutions, and display defaults. Tree traversal iterators (BFS, DFS pre/post) are provided on the node type.

The `auspice_to_graph()` and `auspice_from_graph()` functions ([packages/treetime-io/src/auspice.rs](../../packages/treetime-io/src/auspice.rs)) convert between the Auspice tree model and the generic `Graph` structure. The converter adapter ([packages/treetime-cli/src/convert/auspice.rs](../../packages/treetime-cli/src/convert/auspice.rs)) handles divergence accumulation (Auspice stores cumulative divergence, the graph stores per-edge branch lengths) and mutation formatting.

## Supporting formats

### PhyloGraph JSON

An internal graph serialization format that round-trips the full `ConverterGraph` structure (nodes, edges, graph-level data) through serde JSON. No custom format logic - the graph type derives `Serialize`/`Deserialize` and uses the generic `json_read_file()`/`json_write_file()` helpers from `treetime-io`. Useful for debugging and for lossless intermediate storage during multi-step format conversions. The `clock` command writes this format as `graph_output.json`.

### Graphviz DOT

The `graphviz_write_file()` function ([packages/treetime-io/src/graphviz.rs](../../packages/treetime-io/src/graphviz.rs)) generates directed graph output with subgraphs for roots, internal nodes, and leaves. Used by the `clock` command for tree visualization (`graph_output.dot`). Output only, not available in the convert command.

## I/O architecture

All format adapters convert through a common intermediate representation: `ConverterGraph`, a `Graph<ConverterNode, ConverterEdge, ConverterData>` ([packages/treetime-cli/src/convert/convert.rs#L25-L69](../../packages/treetime-cli/src/convert/convert.rs#L25-L69)) where nodes carry names, edges carry branch lengths and optional mutations, and graph-level data carries Auspice metadata and root sequence information. Format-specific traits (`AuspiceRead`/`AuspiceWrite`, `UsherRead`/`UsherWrite`, `PhyloxmlToGraph`/`PhyloxmlFromGraph`, `NodeFromNwk`/`NodeToNwk`) handle the conversion between each format's type model and the common graph. Adding a new format requires implementing only the format-specific traits, not N-to-N pairwise converters.

Transparent compression is supported on both input and output: gz, bz2, xz, and zst formats are auto-detected by file extension.

The `convert` subcommand exposes all eight format adapters as a standalone tool. The `TreeFormat` enum ([packages/treetime-cli/src/convert/args.rs#L8-L18](../../packages/treetime-cli/src/convert/args.rs#L8-L18)) defines `Auspice`, `MatJson`, `MatPb`, `Newick`, `Nexus`, `PhyloGraph`, `Phyloxml`, `PhyloxmlJson` with auto-detection from file extensions and explicit `--input-format`/`--output-format` overrides.

### Current limitations

- Analysis commands (ancestral, clock, timetree) still read Newick only and write Newick + Nexus - the new format adapters are not yet integrated into the analysis pipeline
- Nexus reading is not implemented (`convert_read_file` panics with `unimplemented!()` for `TreeFormat::Nexus`)
- Nexus extension auto-detection has a bug: the match arm `"nex | nexus"` is a literal string rather than separate match arms, so `.nex`/`.nexus` files are not auto-detected
- Graphviz DOT is not available as a convert command format

## Practical impact

- Trees from UShER pandemic surveillance can be imported for molecular clock inference after conversion to Newick
- Auspice JSON datasets can be read back for re-analysis, not just produced as one-way output
- PhyloXML's rich annotation structure (taxonomy, sequences, events, distributions) is available for TreeTime output
- Users can convert between formats without external tools: `treetime convert input.mat.pb -o output.phylo.xml`
- Future integration of the format adapters into analysis commands will allow direct input from any supported format

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
