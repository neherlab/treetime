#[cfg(test)]
mod tests {
  use crate::tree_ir::mutation::{NUC_GENE, TreeIrSub};
  use pretty_assertions::assert_eq;
  use treetime_primitives::AsciiChar;

  fn nuc(c: char) -> AsciiChar {
    AsciiChar::try_from_char(c).unwrap()
  }

  #[test]
  fn test_tree_ir_mutation_nuc_string_roundtrip() {
    let sub = TreeIrSub {
      gene: NUC_GENE.to_owned(),
      position: 123,
      parent: nuc('A'),
      child: nuc('T'),
    };
    assert_eq!("A123T", sub.to_auspice_string());
    assert_eq!(sub, TreeIrSub::from_auspice_string(NUC_GENE, "A123T").unwrap());
  }

  #[test]
  fn test_tree_ir_mutation_aa_string_roundtrip() {
    let sub = TreeIrSub {
      gene: "M".to_owned(),
      position: 1,
      parent: nuc('M'),
      child: nuc('I'),
    };
    assert_eq!("M1I", sub.to_auspice_string());
    assert_eq!(sub, TreeIrSub::from_auspice_string("M", "M1I").unwrap());
  }

  #[test]
  fn test_tree_ir_mutation_too_short_is_error() {
    let _error = TreeIrSub::from_auspice_string(NUC_GENE, "A").unwrap_err();
  }

  #[test]
  fn test_tree_ir_mutation_non_numeric_position_is_error() {
    let _error = TreeIrSub::from_auspice_string(NUC_GENE, "AxT").unwrap_err();
  }
}
