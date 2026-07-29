# Timetree molecular-clock text output lacks a replacement decision

TreeTime v0 writes `molecular_clock.txt`; v1 writes a structured clock-model JSON document. No approved decision establishes that the JSON output preserves the text file's information and workflow interoperability or intentionally replaces its public contract.

Compare both schemas and downstream uses. Implement the text contract unless an approved parity decision establishes the JSON output as its replacement.
