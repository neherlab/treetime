use eyre::Report;
use std::collections::BTreeMap;
use std::path::PathBuf;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::FastaWriter;
use treetime_primitives::Seq;
use treetime_utils::make_internal_error;

pub(crate) struct NucFastaWriters {
  partitions: BTreeMap<usize, NucFasta>,
}

impl NucFastaWriters {
  pub(crate) fn create(
    paths: BTreeMap<usize, PathBuf>,
    mut descs: BTreeMap<usize, BTreeMap<GraphNodeKey, Option<String>>>,
  ) -> Result<Self, Report> {
    let partitions = paths
      .into_iter()
      .map(|(partition, path)| {
        let fasta = NucFasta {
          writer: FastaWriter::create(&path)?,
          descs: descs.remove(&partition).unwrap_or_default(),
          path,
        };
        Ok((partition, fasta))
      })
      .collect::<Result<_, Report>>()?;
    Ok(Self { partitions })
  }

  pub(crate) fn write(&mut self, partition: usize, key: GraphNodeKey, name: &str, seq: &Seq) -> Result<(), Report> {
    let Some(fasta) = self.partitions.get_mut(&partition) else {
      return make_internal_error!("No reconstructed nucleotide FASTA was planned for partition {partition}");
    };
    let desc = fasta.descs.get(&key).and_then(Option::as_deref);
    fasta.writer.write(name, desc, seq)
  }

  pub(crate) fn finish(self) -> Result<Vec<PathBuf>, Report> {
    self
      .partitions
      .into_values()
      .map(|fasta| {
        fasta.writer.finish()?;
        Ok(fasta.path)
      })
      .collect()
  }
}

struct NucFasta {
  path: PathBuf,
  writer: FastaWriter,
  descs: BTreeMap<GraphNodeKey, Option<String>>,
}
