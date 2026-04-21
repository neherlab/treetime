# Stdin FASTA path reads one record and truncates multi-record alignments

The stdin FASTA input path in the ancestral command reads exactly one record from standard input and silently discards the rest. A multi-record alignment piped via stdin is truncated to its first sequence.

## Impact

Users piping FASTA from another tool (common in bioinformatics pipelines) get a single-sequence "alignment" without any error or warning. Ancestral reconstruction on a single sequence produces degenerate results.

## Fix

Read all records from stdin until EOF, matching the file-path behavior. Add a test that pipes a multi-record alignment through stdin and verifies all records are present in the parsed alignment.
