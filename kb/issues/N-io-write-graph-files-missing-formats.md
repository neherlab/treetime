# shared graph writer missing PhyloXML, Auspice, and UShER MAT formats

The shared graph writer produces four output formats (Newick, NEXUS, JSON, Graphviz) for every tree-outputting command. The converter tool supports three additional formats that are not included: PhyloXML, Auspice JSON, and UShER MAT (JSON and protobuf).

## Constraints

Each missing format requires adapter traits that are not implemented on all graph payload types:

- PhyloXML: requires `PhyloxmlDataFromGraphData + PhyloxmlFromGraph<N, E, D>`. Currently only implemented for the converter's payload type.
- Auspice JSON: requires a command-specific `AuspiceWriter` trait adapter. Timetree has its own `write_auspice_json` with domain logic (confidence intervals, mutations).
- UShER MAT: requires a command-specific `UsherWriter` trait adapter. Only meaningful for converter and ancestral-like payloads.

## Options

- Implement `PhyloxmlFromGraph` for `NodeAncestral`/`EdgeAncestral`, `NodeTimetree`/`EdgeTimetree`, and `NodeClock`/`EdgeClock`, then add PhyloXML to the shared graph writer
- Auspice and UShER MAT may not be feasible as generic outputs due to command-specific adapter logic

## Locations

- Shared helper: `packages/treetime-io/src/graph.rs`
- PhyloXML traits: `packages/treetime-io/src/phyloxml.rs`
- Auspice traits: `packages/treetime-io/src/auspice.rs`
- UShER MAT traits: `packages/treetime-io/src/usher_mat.rs`
- Converter format enum: `packages/treetime-cli/src/convert/args.rs` (`TreeFormat`)

## Related

- [kb/proposals/output-format-selection.md](../proposals/output-format-selection.md) -- three-tier output format selection subsumes this issue
- [kb/tickets/io-output-selection-infrastructure.md](../tickets/io-output-selection-infrastructure.md) -- output selection infrastructure
- [kb/tickets/io-format-adapter-impls.md](../tickets/io-format-adapter-impls.md) -- format adapter implementations
- [N-timetree-auspice-json-incomplete.md](N-timetree-auspice-json-incomplete.md) - timetree auspice output gaps
- [N-ancestral-auspice-json-not-produced.md](N-ancestral-auspice-json-not-produced.md) - ancestral auspice output missing entirely
