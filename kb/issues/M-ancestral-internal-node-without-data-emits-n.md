# Internal nodes emit `N` where no descendant has data; v0 emits a base

At an alignment column where every descendant leaf of an internal node is `N`, v1 writes `N` for that internal node in both dense and sparse reconstruction. v0 writes the most likely base, which with no data is the most likely state under the model's equilibrium frequencies. With the ambiguous tip states kept (augur `--keep-ambiguous`), the two engines therefore report the changes to `N` on different branches, and their mutation counts differ by up to a factor of 7.6 on Nextstrain builds.

## Reconstructed states

Tree `((L1,L2)X,L3)root` with columns 3-4 of every leaf set to `N`:

| Engine            | X    | root |
| ----------------- | ---- | ---- |
| v0 (`--gtr JC69`) | `CC` | `CC` |
| v1 dense          | `NN` | `NN` |
| v1 sparse         | `NN` | `NN` |

- v0 uses the 5-state `nuc` alphabet in which `N` is the all-ones profile ([packages/legacy/treetime/treetime/seq_utils.py#L20-L35](../../packages/legacy/treetime/treetime/seq_utils.py#L20-L35)) and assigns internal sequences by `prof2seq` argmax ([packages/legacy/treetime/treetime/treeanc.py#L909-L919](../../packages/legacy/treetime/treetime/treeanc.py#L909-L919)), so it never emits `N` at an internal node
- v1 marks such columns as unknown from the children (`fn compute_node_ranges` in [packages/treetime/src/seq/indel.rs](../../packages/treetime/src/seq/indel.rs)) and fills them with `N`

[kb/decisions/ancestral-marginal-tip-reconstruction-and-imputation.md](../decisions/ancestral-marginal-tip-reconstruction-and-imputation.md) describes this behavior ("Internal unknown states"), but no decision compares it with v0 or approves the difference.

## Effect on the mutation lists

With `--report-ambiguous`, which keeps the changes to `N` in the outputs, v1 reports one change to `N` on the branch above the unknown internal node, where v0 reports one change to `N` on the branch of every tip below it. In a four-tip tree whose tips `A` and `B` both lack the first two bases:

```bash
printf '(((A:0.1,B:0.1):0.1,C:0.1):0.05,D:0.3);\n' > tree.nwk
printf '>A\nNNGAAATTTCCCGGG\n>B\n---AAATTTCCCGGG\n>C\nATGAAATTTCCCG--\n>D\nATGAAATTTCCAGNN\n' > aln.fasta
treetime ancestral --tree tree.nwk --alignment aln.fasta --method-anc marginal --report-ambiguous \
  --gap-fill only-terminal --output-node-data v1.json
```

v1 writes `A1N` and `T2N` on the parent of `A` and `B`, and `G3N` on `B` (whose terminal gaps `--gap-fill only-terminal` fills with `N`). `augur ancestral --keep-ambiguous` with v0 writes `A1N` and `T2N` on `A`, `A1N`, `T2N`, and `G3N` on `B`, and nothing on their parent. Dense and sparse reconstructions of v1 give the same lists.

The Nextstrain H5N1 cattle-outbreak builds call `augur ancestral --keep-ambiguous`, and their sequences share long runs of missing data within clades. The mutations without gap changes, summed over the nodes of both trees:

| Build | Mutations, v0 | Mutations, v1 | Shared | Jaccard index |
| --- | ---: | ---: | ---: | ---: |
| cattle outbreak, genome | 2 599 536 | 343 834 | 234 715 | 0.087 |
| cattle outbreak, PB2 | 344 234 | 88 108 | 75 423 | 0.21 |
| cattle outbreak, HA | 297 333 | 82 679 | 73 470 | 0.24 |
| D1.1, genome | 496 933 | 472 665 | 462 561 | 0.91 |

The D1.1 build uses the same setting and agrees, so the size of the effect follows how much missing data the tips of a clade share.

## Decision required

> [!IMPORTANT]
> **Decision required.** How an internal node is reconstructed at a column where no descendant has data. Options:
>
> - **Keep `N`** and record the difference from v0 in `kb/decisions/`: no state is inferred without data, and a change to `N` appears once per clade. Mutation lists with ambiguous tip states then differ from augur
> - **Emit the most likely base as v0 does**: the model's equilibrium state, and the changes to `N` appear on every tip, as augur reports them
> - **Keep `N` in the reconstruction, and write the changes to `N` per tip in the augur outputs**: the outputs match augur, while the sequences keep `N`; the written mutations then no longer rebuild the internal sequences from the root
>
> Evidence: the tables above. The rule that the root plus the written mutations rebuilds every node ([kb/decisions/ancestral-marginal-tip-reconstruction-and-imputation.md](../decisions/ancestral-marginal-tip-reconstruction-and-imputation.md), "One filter") holds for the first two options.
