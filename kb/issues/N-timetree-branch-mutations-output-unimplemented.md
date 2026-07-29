# Timetree branch-mutation table output is unimplemented

TreeTime v0 writes `branch_mutations.txt`, a tabular mapping of mutations to branches. V1 can emit mutation annotations in tree formats, but it does not write the separate table required by workflows that consume the v0 file contract.

Trace the v0 columns, ordering, coordinate convention, and ambiguity handling before creating an implementation ticket. Tree annotations do not establish tabular interoperability.
