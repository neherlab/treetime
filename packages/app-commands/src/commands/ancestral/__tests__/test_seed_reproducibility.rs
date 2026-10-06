#[cfg(test)]
mod tests {
  use eyre::{OptionExt, Report, WrapErr};
  use pretty_assertions::assert_eq;
  use std::fs;
  use tempfile::tempdir;

  const TREE: &str = "((A:0.4,B:0.4)AB:0.4,(C:0.4,D:0.4)CD:0.4)root;";
  const ALIGNMENT: &str = ">A\nACGTACGTAC\n>B\nCAGTTCGAAC\n>C\nGCTAACTTGA\n>D\nTGCATGCAGT\n";
  const SEED_PREFIX: &str = "Sampling from the profile is stochastic; seed ";

  #[test]
  fn test_seed_reproducibility_logged_seed_reproduces_unseeded_sampling() -> Result<(), Report> {
    let dir = tempdir().wrap_err("When creating a temporary directory")?;
    fs::write(dir.path().join("tree.nwk"), TREE).wrap_err("When writing the tree fixture")?;
    fs::write(dir.path().join("aln.fasta"), ALIGNMENT).wrap_err("When writing the alignment fixture")?;

    let (unseeded_fasta, unseeded_messages) = helpers::run(dir.path(), "unseeded", None)?;
    let seed = unseeded_messages
      .iter()
      .find_map(|message| message.strip_prefix(SEED_PREFIX))
      .and_then(|rest| rest.split_whitespace().next())
      .ok_or_eyre("an unseeded sampling run must log its seed")?
      .parse::<u64>()
      .wrap_err("When parsing the logged seed")?;
    let (seeded_fasta, seeded_messages) = helpers::run(dir.path(), "seeded", Some(seed))?;

    let seed_lines = |messages: &[String]| {
      messages
        .iter()
        .filter(|message| message.starts_with(SEED_PREFIX))
        .cloned()
        .collect::<Vec<_>>()
    };
    assert_eq!(unseeded_fasta, seeded_fasta);
    assert_eq!(seed_lines(&unseeded_messages), seed_lines(&seeded_messages));
    Ok(())
  }

  mod helpers {
    use crate::commands::ancestral::args::{SampleModeCli, TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
    use crate::commands::ancestral::run::run_ancestral_reconstruction;
    use crate::commands::shared::alignment::AlignmentArgs;
    use crate::commands::shared::model::{GtrModelNameCli, ModelArgs};
    use crate::commands::shared::seed::SeedArgs;
    use crate::job::{JobEvent, JobProgress};
    use eyre::{Report, WrapErr};
    use parking_lot::Mutex;
    use std::fs;
    use std::path::Path;
    use treetime::cancel::NoopCancel;
    use treetime::progress::{LogEvent, NoopProgress};

    pub(super) fn run(dir: &Path, name: &str, seed: Option<u64>) -> Result<(String, Vec<String>), Report> {
      let fasta_path = dir.join(format!("{name}.fasta"));
      let args = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![dir.join("aln.fasta")],
        },
        tree: Some(dir.join("tree.nwk")),
        model_args: ModelArgs {
          model: GtrModelNameCli::JC69,
          ..ModelArgs::default()
        },
        include_leaves: true,
        output_reconstructed_nuc_fasta: Some(fasta_path.clone()),
        sample_from_profile: SampleModeCli::All,
        seed_args: SeedArgs { seed },
        ..TreetimeAncestralArgsRaw::default()
      })?;
      let messages = Mutex::new(vec![]);
      let log = JobProgress::new(|event: JobEvent| {
        if let JobEvent::Log {
          data: LogEvent { message, .. },
        } = event
        {
          messages.lock().push(message);
        }
      });
      run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &log)?;
      let fasta = fs::read_to_string(&fasta_path).wrap_err("When reading the reconstructed sequences")?;
      Ok((fasta, messages.into_inner()))
    }
  }
}
