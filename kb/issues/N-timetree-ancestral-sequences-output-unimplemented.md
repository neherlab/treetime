# Timetree ancestral-sequence FASTA output is unimplemented

TreeTime v0 writes `ancestral_sequences.fasta` from reconstructed internal and tip sequences. The v1 timetree command can reconstruct ancestral sequences through its partitions but does not expose the equivalent FASTA output.

Exact output parity remains the approved default. Define the sequence-selection, ambiguity, and compression projection contract against v0 before creating an implementation ticket.
