# Convert command description is stale

The `convert` command's help text says "Read and write Usher MAT files" (`packages/treetime-cli/src/convert/args.rs` doc comment on `pub struct Args`), but the command is a general tree-format converter supporting 8 formats: Auspice, MatJson, MatPb, Newick, Nexus, PhyloGraph, Phyloxml, and PhyloxmlJson.

The description should reflect the full set of supported formats.
