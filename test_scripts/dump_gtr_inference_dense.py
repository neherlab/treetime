#!/usr/bin/env python3
"""Generate golden master fixture for GTR inference from dense representation.

Runs marginal ancestral reconstruction, extracts sequences via argmax from profiles,
counts mutations, and infers GTR model. Outputs JSON fixture for Rust tests.

Usage:
    ./dev/docker/python test_scripts/dump_gtr_inference_dense.py
"""
from __future__ import annotations

import json
import logging
import lzma
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from Bio import AlignIO

from ancestral_dense import ancestral_dense, init_sequences_dense
from infer_gtr import infer_gtr
from lib import graph_from_nwk_str
from payload import EdgePayload, NodePayload
from treetime import GTR

if TYPE_CHECKING:
    from lib import Graph

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)


def prof2seq(profile: np.ndarray, alphabet: str = "ACGT") -> str:
    """Extract sequence from profile matrix via argmax."""
    return "".join(alphabet[i] for i in np.argmax(profile, axis=1))


def count_composition(seq: str, alphabet: str = "ACGT") -> np.ndarray:
    """Count character composition of sequence."""
    composition = np.zeros(len(alphabet))
    for ch in seq:
        if ch in alphabet:
            composition[alphabet.index(ch)] += 1.0
    return composition


def get_mutation_counts_dense(graph: Graph) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract mutation counts from dense profiles (argmax sequences)."""
    root = graph.get_one_root()
    root_profile = root.payload().dense_sequences[0].profile.dis
    root_seq = prof2seq(root_profile)
    root_state = count_composition(root_seq)

    n = 4  # ACGT
    nij = np.zeros((n, n))
    ti = np.zeros(n)
    char_to_index = {nuc: i for i, nuc in enumerate("ACGT")}

    for edge in graph.get_edges():
        branch_length = edge.payload().branch_length or 0.0
        source = graph.get_node(edge.source())
        target = graph.get_node(edge.target())

        parent_profile = source.payload().dense_sequences[0].profile.dis
        child_profile = target.payload().dense_sequences[0].profile.dis

        parent_seq = prof2seq(parent_profile)
        child_seq = prof2seq(child_profile)

        child_composition = count_composition(child_seq)
        ti += child_composition * branch_length

        for pos, (p_ch, c_ch) in enumerate(zip(parent_seq, child_seq, strict=True)):
            if p_ch != c_ch and p_ch in char_to_index and c_ch in char_to_index:
                i = char_to_index[p_ch]  # parent (ref)
                j = char_to_index[c_ch]  # child (qry)
                nij[j, i] += 1.0
                ti[j] -= 0.5 * branch_length
                ti[i] += 0.5 * branch_length

    return nij, ti, root_state


def main() -> None:
    data_dir = Path("data/flu/h3n2/20")
    tree_path = data_dir / "tree.nwk"
    aln_path = data_dir / "aln.fasta.xz"
    output_path = Path("packages/treetime/src/gtr/__tests__/fixtures/gtr_inference_dense.json")

    log.info("Loading tree from %s", tree_path)
    nwk_str = tree_path.read_text()
    graph = graph_from_nwk_str(
        nwk_string=nwk_str,
        node_payload_factory=NodePayload,
        edge_payload_factory=EdgePayload,
    )

    log.info("Loading alignment from %s", aln_path)
    with lzma.open(aln_path, "rt") as fh:
        aln = {seq.id: str(seq.seq).upper() for seq in AlignIO.read(fh, "fasta")}

    log.info("Running marginal ancestral reconstruction with JC69 (nuc_nogap)")
    gtr = GTR.standard(model="JC69", alphabet="nuc_nogap")
    init_sequences_dense(graph, [aln], [gtr])
    log_lh = ancestral_dense(graph)
    log.info("LogLH after reconstruction: %f", log_lh)

    log.info("Extracting mutation counts from dense profiles")
    nij, ti, root_state = get_mutation_counts_dense(graph)

    log.info("Inferring GTR model")
    result = infer_gtr(nij, ti, root_state, pc=1.0)

    fixture = {
        "description": "Golden master for GTR inference from dense representation",
        "dataset": "flu/h3n2/20",
        "initial_gtr": "JC69",
        "mutation_counts": {
            "nij": nij.tolist(),
            "Ti": ti.tolist(),
            "root_state": root_state.tolist(),
        },
        "inferred_gtr": {
            "W": result["W"].tolist(),
            "pi": result["pi"].tolist(),
            "mu": float(result["mu"]),
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    log.info("Writing fixture to %s", output_path)
    with output_path.open("w") as fh:
        json.dump(fixture, fh, indent=2)

    log.info("nij:\n%s", nij)
    log.info("Ti: %s", ti)
    log.info("root_state: %s", root_state)
    log.info("W:\n%s", result["W"])
    log.info("pi: %s", result["pi"])
    log.info("mu: %s", result["mu"])


if __name__ == "__main__":
    main()
