#![allow(clippy::len_without_is_empty)]

use crate::io::fasta::read_many_fasta;
use eyre::Report;
use itertools::Itertools;
use ndarray::{s, Array1, Array2, ArrayBase, ArrayView1, Axis, Data, Ix1};
use polars::export::arrow::array::Array;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::Path;

#[derive(Clone, Debug)]
pub struct Sequence {
  index: usize,
  seq_name: String,
  seq: Array1<char>,
}

#[derive(Clone, Debug)]
pub struct SequenceCompressed {
  index: usize,
  seq_name: String,
  compressed_alignment: Array1<char>,
  compressed_to_full_sequence_map: Vec<Array1<i32>>,
}

#[derive(Clone, Debug, Default)]
pub struct AlignmentPattern {
  length: usize,
  patterns: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct SequenceData {
  aln: Array2<char>,
  aln_transpose: Array2<char>,
  aln_compressed: Array2<char>,
  records: Vec<Sequence>,
  records_compressed: Vec<SequenceCompressed>,
  seq_names_to_indices: HashMap<String, usize>,
  len_full: usize,
  len_compressed: usize,
  multiplicity: Array1<f64>,
  compression_map: Array1<usize>,
  decompression_map: HashMap<usize, Array1<usize>>,
}

impl SequenceData {
  pub fn new<P: AsRef<Path>>(fasta_paths: &[P], ambiguous_letter: char) -> Result<Self, Report> {
    let fasta_records = read_many_fasta(fasta_paths)?;

    // TODO: does this need to be mode complex?
    let mut len_full = fasta_records[0].seq.len();

    let seq_names_to_indices: HashMap<String, usize> = fasta_records
      .iter()
      .map(|record| (record.seq_name.clone(), record.index))
      .collect();

    let records = fasta_records
      .iter()
      .map(|record| Sequence {
        index: record.index,
        seq_name: record.seq_name.clone(),
        seq: Array1::<char>::from_iter(record.seq.chars()),
      })
      .collect_vec();

    // Combine all sequences into a matrix of shape (num_seq x seq_len)
    let aln_vec: Vec<ArrayView1<char>> = records.iter().map(|record| record.seq.view()).collect();
    let aln: Array2<char> = ndarray::stack(Axis(0), &aln_vec)?;

    // Transpose to shape shape (seq_len x num_seq)
    let aln_transpose: Array2<char> = aln.t().to_owned();

    // TODO: `additional_constant_sites`: where does it come from?
    let additional_constant_sites: Option<usize> = None;

    #[allow(clippy::self_assignment)]
    if let Some(additional_constant_sites) = additional_constant_sites {
      add_additional_constant_sites();
      // TODO: full length should change here
      len_full = len_full;
    }

    let (aln_compressed, alignment_patterns) = compress(&aln_transpose, ambiguous_letter)?;

    let multiplicity = compute_multiplicity(&alignment_patterns);

    let compression_map = compute_compression_map(&alignment_patterns, len_full);
    let decompression_map = compute_decompression_map(&alignment_patterns);

    Ok(Self {
      aln,
      aln_transpose,
      aln_compressed,
      records,
      records_compressed: vec![],
      seq_names_to_indices,
      len_full,
      len_compressed: multiplicity.len(),
      multiplicity,
      compression_map,
      decompression_map,
    })
  }

  #[inline]
  pub const fn len_full(&self) -> usize {
    self.len_full
  }

  #[inline]
  pub const fn len_compressed(&self) -> usize {
    self.len_compressed
  }

  /// Inverse of the uncompressed sequence length. Length scale for short branches.
  #[inline]
  pub fn one_mutation(&self) -> f64 {
    1.0 / (self.len_full() as f64)
  }

  #[inline]
  pub fn multiplicity(&self) -> ArrayView1<f64> {
    self.multiplicity.view()
  }

  /// Retrieve compressed sequence by name
  pub fn get_compressed(&self, seq_name: &str) -> Option<ArrayView1<char>> {
    self
      .seq_names_to_indices
      .get(seq_name)
      .map(|index| self.aln_compressed.slice(s![*index, ..]))
  }

  /// Retrieve full sequence by name
  pub fn get_full(&self, seq_name: &str) -> Option<Array1<char>> {
    self
      .get_compressed(seq_name)
      .map(|compressed| self.decompress(&compressed))
  }

  pub fn decompress<S>(&self, compressed_seq: &ArrayBase<S, Ix1>) -> Array1<char>
  where
    S: Data<Elem = char>,
  {
    // TODO: account for `additional_constant_sites`
    // if include_additional_constant_sites:
    //     L = self.full_length
    // else:
    //     L = self.full_length - self.additional_constant_sites
    self.compression_map.iter().map(|i| compressed_seq[*i]).collect()
  }
}

fn compress(
  aln_transpose: &Array2<char>,
  ambiguous_letter: char,
) -> Result<(Array2<char>, HashMap<String, AlignmentPattern>), Report> {
  let mut alignment_patterns = HashMap::<String, AlignmentPattern>::new();
  let mut compressed_aln_transpose: Vec<Array1<char>> = vec![];
  let variable_positions = Array1::<usize>::from_iter(0..aln_transpose.shape()[0]);
  for pi in variable_positions {
    let mut pattern = aln_transpose.slice(s!(pi, ..)).to_owned();

    // if the column contains only one state and ambiguous nucleotides, replace
    // those with the state in other strains right away
    let unique_letters = replace_ambiguous_in_place(&mut pattern, ambiguous_letter);

    let mut str_pattern = pattern.iter().join("");
    if unique_letters.len() > 1 {
      let _ = write!(str_pattern, "_{pi}");
    }

    let alignment_pattern = alignment_patterns.entry(str_pattern).or_insert_with(|| {
      let length = compressed_aln_transpose.len();
      compressed_aln_transpose.push(pattern);
      AlignmentPattern {
        length,
        patterns: vec![],
      }
    });
    alignment_pattern.patterns.push(pi);
  }

  let aln_compressed_transposed = ndarray::stack(
    Axis(0),
    &compressed_aln_transpose.iter().map(ArrayBase::view).collect_vec(),
  )?;

  let aln_compressed = aln_compressed_transposed.t().to_owned();

  Ok((aln_compressed, alignment_patterns))
}

fn replace_ambiguous_in_place(pattern: &mut Array1<char>, ambiguous_letter: char) -> Vec<char> {
  let mut unique_letters = pattern.iter().unique().copied().collect_vec();
  if unique_letters.len() == 2 && unique_letters.contains(&ambiguous_letter) {
    let other_letter = unique_letters
      .into_iter()
      .filter(|c| c == &ambiguous_letter)
      .collect_vec()[0];
    pattern.mapv_inplace(|c| if c == ambiguous_letter { other_letter } else { c });
    unique_letters = vec![other_letter];
  }
  unique_letters
}

/// Count how many times each column is repeated in the real alignment
fn compute_multiplicity(alignment_patterns: &HashMap<String, AlignmentPattern>) -> Array1<f64> {
  let mut multiplicity = Array1::<usize>::zeros(alignment_patterns.len());
  for AlignmentPattern { length, patterns } in alignment_patterns.values() {
    multiplicity[*length] = patterns.len();
  }
  multiplicity.mapv(|x| x as f64)
}

// Bind positions in full length sequence to that of the compressed sequence
fn compute_compression_map(
  alignment_patterns: &HashMap<String, AlignmentPattern>,
  full_length: usize,
) -> Array1<usize> {
  let mut compression_map = Array1::<usize>::zeros(full_length);
  for AlignmentPattern { length, patterns } in alignment_patterns.values() {
    for i in patterns {
      compression_map[*i] = *length;
    }
  }
  compression_map
}

// Bind position in compressed sequence to the array of positions in full length sequence
fn compute_decompression_map(alignment_patterns: &HashMap<String, AlignmentPattern>) -> HashMap<usize, Array1<usize>> {
  let mut decompression_map = HashMap::<usize, Array1<usize>>::new();
  for AlignmentPattern { length, patterns } in alignment_patterns.values() {
    decompression_map.insert(*length, Array1::<usize>::from_iter(patterns.clone()));
  }
  decompression_map
}

fn add_additional_constant_sites() {
  unimplemented!("add_additional_constant_sites: not yet implemented");
}
