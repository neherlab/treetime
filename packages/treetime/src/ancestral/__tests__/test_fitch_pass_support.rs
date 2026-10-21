#[cfg(test)]
pub(super) mod tests {
  pub(crate) mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::partition::fitch::partition::PartitionFitch;
    use crate::seq::alignment::NodeSeqInput;
    use crate::test_utils::{get_common_length, leaf_seq_inputs};
    use eyre::Report;
    use indoc::indoc;
    use maplit::btreemap;
    use rayon::ThreadPoolBuilder;
    use std::collections::BTreeMap;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::node_names::NodeNames;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::AlignmentRecord;

    pub(crate) struct FitchInput {
      pub graph: Graph,
      pub names: NodeNames,
      pub partition: PartitionFitch,
      pub node_inputs: BTreeMap<GraphNodeKey, NodeSeqInput>,
    }

    pub(crate) fn fitch_input() -> Result<FitchInput, Report> {
      let alphabet = Alphabet::default();
      let aln: Vec<AlignmentRecord> = fasta_read(
        indoc! {r#"
          >A
          ACATCGCCNNA--GAC
          >B
          GCATCCCTGTA-NG--
          >C
          CCGGCGATGTRTTG--
          >D
          TCGGCCGTGTRTTG--
        "#}
        .as_bytes(),
        &alphabet,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
      let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;".as_slice())?;
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let partition = PartitionFitch {
        alphabet,
        length: get_common_length(&aln)?,
        nodes: btreemap! {},
        edges: btreemap! {},
      };
      let node_inputs = leaf_seq_inputs(&graph, &names, aln);
      Ok(FitchInput {
        graph,
        names,
        partition,
        node_inputs,
      })
    }

    pub(crate) fn run_in_pool<T: Send>(
      threads: usize,
      run: impl FnOnce() -> Result<T, Report> + Send,
    ) -> Result<T, Report> {
      ThreadPoolBuilder::new().num_threads(threads).build()?.install(run)
    }
  }
}
