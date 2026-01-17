# phyloxml

Rust library for reading and writing PhyloXML files - an XML format for phylogenetic trees and associated data.

## Overview

Parse and serialize phylogenetic data according to the PhyloXML 1.20 schema (http://www.phyloxml.org/1.20/phyloxml.xsd). The implementation covers common elements but is not exhaustive - unknown fields are preserved via serde's flatten mechanism.

## Features

- Read PhyloXML from any `std::io::Read` source
- Write PhyloXML with pretty-printed indentation to any `std::io::Write` sink
- Preserve unrecognized elements through `BTreeMap<String, serde_json::Value>` fields
- Support for core PhyloXML elements:
  - Phylogeny metadata (name, description, rooted status)
  - Clade hierarchy with branch lengths, confidence values, colors
  - Taxonomy (scientific names, common names, ranks)
  - Sequences with molecular data and domain architecture
  - Geographic distributions with coordinates
  - Dates, events, properties, and references

## Usage

```rust
use phyloxml::{phyloxml_read, phyloxml_write, Phyloxml};
use std::fs::File;

// Read
let file = File::open("tree.xml")?;
let data: Phyloxml = phyloxml_read(file)?;

for phylogeny in &data.phylogeny {
    println!("Tree: {:?}, rooted: {}", phylogeny.name, phylogeny.rooted);
}

// Write
let out = File::create("output.xml")?;
phyloxml_write(out, &data)?;
```

## API

### Functions

- `phyloxml_read(reader) -> Result<Phyloxml, DeError>` - Deserialize PhyloXML from a reader
- `phyloxml_write(writer, &Phyloxml) -> std::io::Result<()>` - Serialize to pretty-printed XML

### Key Types

- `Phyloxml` - Root container with vector of phylogenies
- `PhyloxmlPhylogeny` - Single phylogenetic tree with metadata
- `PhyloxmlClade` - Tree node with optional children, branch length, sequences, taxonomy
- `PhyloxmlSequence` - Molecular sequence data with accession and annotations
- `PhyloxmlTaxonomy` - Taxonomic classification
- `PhyloxmlDate` - Temporal information with optional ranges
- `PhyloxmlConfidence` - Support values with type annotation

## Dependencies

- `quick-xml` - XML parsing and serialization with serde support
- `serde` - Serialization framework
- `xml` - Pretty-printing output
