# Timetree sequence-evolution model text output lacks a replacement decision

TreeTime v0 writes `sequence_evolution_model.txt`; v1 writes `gtr.json` when a model is available. No approved decision establishes that the JSON output preserves the text file's information and workflow interoperability or intentionally replaces its public contract.

Compare both schemas and downstream uses. Implement the text contract unless an approved parity decision establishes `gtr.json` as its replacement.
