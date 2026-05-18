//! Tests for Fitch GTR inference.

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::gtr::infer_gtr::common::{InferGtrOptions, infer_gtr_impl};
  use crate::ancestral::gtr_inference::get_mutation_counts_fitch;
  use crate::pretty_assert_ulps_eq;
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use ndarray::array;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
  }

  #[test]
  fn test_get_mutation_counts_fitch() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;

    let counts_actual = get_mutation_counts_fitch(&graph, &fitch)?;
    pretty_assert_ulps_eq!(
      counts_actual.nij,
      array![[0., 0., 0., 0.], [2., 0., 0., 1.], [3., 2., 0., 0.], [0., 1., 1., 0.]],
      epsilon = 1e-9
    );
    pretty_assert_ulps_eq!(counts_actual.root_state, array![4.0, 3.0, 3.0, 4.0], epsilon = 1e-9);
    pretty_assert_ulps_eq!(
      counts_actual.Ti,
      array![
        1.9800000227987766,
        2.9450000282377005,
        2.515000017359853,
        2.640000019222498
      ],
      epsilon = 1e-9
    );

    Ok(())
  }

  #[test]
  fn test_infer_gtr_fitch() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;

    let counts = get_mutation_counts_fitch(&graph, &fitch)?;
    let actual = infer_gtr_impl(
      &counts,
      &InferGtrOptions {
        pc: 0.1,
        ..InferGtrOptions::default()
      },
    )?;

    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      array![
        [0.0,                  2.175112392280934,  2.956016575101041,  0.18620301184748264],
        [2.175112392280934,    0.0,                1.405280909422987,  1.4146569619049196 ],
        [2.956016575101041,    1.405280909422987,  0.0,                0.744903148823837  ],
        [0.18620301184748264,  1.4146569619049196, 0.744903148823837,  0.0                ]
      ],
      &actual.W,
      epsilon = 1e-7
    );

    pretty_assert_ulps_eq!(
      array![
        0.14878846342301447,
        0.2405153630584697,
        0.31239203299951585,
        0.29830414051899984
      ],
      &actual.pi,
      epsilon = 1e-7
    );
    pretty_assert_ulps_eq!(0.9471364348810554, actual.mu, epsilon = 1e-7);

    Ok(())
  }
}
