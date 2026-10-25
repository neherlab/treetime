# Least-squares rerooting gives a different clock on the built SARS-CoV-2 tree

On `sc2/4500`, the tree from the prototype builder and the reference tree give the same clock when both are rooted on Wuhan-Hu-1. They give different clocks under least-squares rerooting ([kb/reports/feat-tree-infer-prototype.md](../reports/feat-tree-infer-prototype.md)):

| Tree | Rooting | Clock rate | Root date |
| --- | --- | --- | --- |
| Reference | outgroup | 1.04e-3 | 2020.02 |
| Built | outgroup | 9.95e-4 | 2019.80 |
| Reference | least squares | 9.1e-4 | 2019.81 |
| Built | least squares | 7.4e-4 | 2019.01 |

Collapsing the empty branches of the reference tree leaves its least-squares estimate unchanged, so the polytomies of the built tree are not the cause. 1364 of 5100 tips have root-to-tip distances that differ by more than 5 substitutions between the two trees.

## Impact

A user who builds the tree natively and runs `timetree` with the default rerooting gets a clock rate about 19% lower and a root about 0.8 years earlier than with the ML tree. It is unknown which estimate is closer to the truth.

## Fix

1. Find which tips move the least-squares root, and whether substitution-reversion pairs on the backbone lengthen their root-to-tip paths.
2. Repeat on simulations with known dates and rates, comparing the built tree, the true topology, and an IQ-TREE tree.
3. If the built tree is biased, consider date-aware tie-breaking among equally parsimonious placements ([kb/proposals/native-parsimony-tree-builder.md](../proposals/native-parsimony-tree-builder.md), open research questions).
