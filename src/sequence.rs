//! Module for sequences

use log::error;
use crate::snps::get_expected_changes_uc;
use anyhow::{bail, Result};
use phf::{phf_map, Map};
use std::collections::HashMap;

/// Universal Genetic Code
pub static UNIVERSAL_CODE: Map<&'static str, &'static str> = phf_map! {
    "TTT" => "F", "TCT" => "S", "TAT" => "Y", "TGT" => "C",
    "TTC" => "F", "TCC" => "S", "TAC" => "Y", "TGC" => "C",
    "TTA" => "L", "TCA" => "S", "TAA" => "*", "TGA" => "*",
    "TTG" => "L", "TCG" => "S", "TAG" => "*", "TGG" => "W",
    "CTT" => "L", "CCT" => "P", "CAT" => "H", "CGT" => "R",
    "CTC" => "L", "CCC" => "P", "CAC" => "H", "CGC" => "R",
    "CTA" => "L", "CCA" => "P", "CAA" => "Q", "CGA" => "R",
    "CTG" => "L", "CCG" => "P", "CAG" => "Q", "CGG" => "R",
    "ATT" => "I", "ACT" => "T", "AAT" => "N", "AGT" => "S",
    "ATC" => "I", "ACC" => "T", "AAC" => "N", "AGC" => "S",
    "ATA" => "I", "ACA" => "T", "AAA" => "K", "AGA" => "R",
    "ATG" => "M", "ACG" => "T", "AAG" => "K", "AGG" => "R",
    "GTT" => "V", "GCT" => "A", "GAT" => "D", "GGT" => "G",
    "GTC" => "V", "GCC" => "A", "GAC" => "D", "GGC" => "G",
    "GTA" => "V", "GCA" => "A", "GAA" => "E", "GGA" => "G",
    "GTG" => "V", "GCG" => "A", "GAG" => "E", "GGG" => "G",
};

/// Alternative, faster
pub fn universal_code(codon: &[u8]) -> Result<u8> {
    let aa = match codon {
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" => b'R',
        b"GGT" | b"GGA" | b"GGC" | b"GGG" => b'G',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        b"TAT" | b"TAC" => b'Y',
        b"TTT" | b"TTC" => b'F',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"TGT" | b"TGC" => b'C',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"AGT" | b"AGC" => b'S',
        b"AGA" | b"AGG" => b'R',
        b"ATG" => b'M',
        b"TGG" => b'W',
        _ => bail!("Cannot find codon {:?}", &codon),
    };
    Ok(aa)
}

/// Binary constant for Adenine
pub const ADENINE: u8 = b'A';
/// Binary constant for Thymine
pub const THYMINE: u8 = b'T';
/// Binary constant for Cytosine
pub const CYTOSINE: u8 = b'C';
/// Binary constant for Guanine
pub const GUANINE: u8 = b'G';

/// All nucleotides
static NUCLEOTIDES: [u8; 4] = [ADENINE, CYTOSINE, THYMINE, GUANINE];

/// map of complements
pub const REV_COMP: Map<u8, u8> = phf_map! {
    b'T' => b'A',
    b'A' => b'T',
    b'C' => b'G',
    b'G' => b'C',
};

/// Given an iterator over a sequence, returns it's
/// complement, based on `REV_COMP`.
pub fn complement<'a, S>(seq: S) -> Vec<u8>
where
    S: Iterator<Item = &'a u8>,
{
    seq.map(|c| REV_COMP[c]).collect()
}

pub type Sequence = Vec<u8>;

/// A Sequence record
///
/// TODO: Expand and use for Qualities?
#[derive(Debug, Default)]
pub struct SequenceRecord {
    // The sequence header, until the first whitespace, `>` removed
    pub id: String,
    // The sequence
    pub seq: Sequence,
    // This keeps the rest of the header
    pub attributes: String,
}

impl SequenceRecord {
    /// Returns if the sequence is empty
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
    /// Checks if the values are `NUCLEOTIDES`
    pub fn is_correct(&self) -> bool {
        self.seq
            .iter()
            .filter(|c| !matches!(**c, ADENINE | THYMINE | GUANINE | CYTOSINE))
            .count()
            == 0
    }
    /// Returns a Vector with the reverse complement of the
    /// sequence.
    pub fn reverse_comp(&self) -> Vec<u8> {
        reverse_comp(&self.seq)
    }
}

/// Returns a Vector with the reverse complement of the
/// sequence.
pub fn reverse_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => {
                error!("Char {} is not valid", c.to_string());
                b'X'
            },
        })
        .collect()
}

/// A structure to hold Aminoacid Changes counts
#[derive(Debug, Default)]
pub struct AaChanges {
    pub syn: u32,
    pub nonsyn: u32,
}

/// Returns the Expected changes matrix
///
/// It simulates that a change occur in the codon and finds how many changes
/// results in a synonymous and non-synonymous aminoacid change.
///
/// TODO: try and use a `&'static str` instead of `String`, compare speed
pub fn get_exp_matrix(trans_table: &Map<&'static str, &'static str>) -> HashMap<String, AaChanges> {
    let mut syn_matrix: HashMap<String, AaChanges> = trans_table
        .keys()
        .map(|codon| (codon.to_string(), AaChanges::default()))
        .collect();

    for codon in trans_table.keys() {
        for phase in 0..=2 {
            for nucleotide in &NUCLEOTIDES {
                // to change one of the characters, it's necessary to convert
                // to a Vec. using .bytes() gets an iterator
                let mut codon2 = Vec::from_iter(codon.bytes());
                codon2[phase] = *nucleotide;
                // now we need to reconstruct a string to perform the comparison
                // and ::from_utf8_lossy is a safe one to use in the case or
                // aa/nuc sequences
                let codon2 = String::from_utf8_lossy(&codon2);
                // no change at all, skip - needs to DeRef the &str
                if *codon == codon2 {
                    continue;
                }
                if trans_table[codon] == trans_table[&codon2] {
                    // no change in aminoacid, so a synonymous
                    syn_matrix.get_mut(*codon).unwrap().syn += 1;
                } else {
                    syn_matrix.get_mut(*codon).unwrap().nonsyn += 1;
                }
            }
        }
    }
    syn_matrix
}

/// Given a sequence and the result of `get_exp_matrix` returns the number
/// of expected synonymous and non-synonymous changes for the sequence.
pub fn get_seq_exp_changes(seq: &[u8]) -> (u32, u32) {
    let mut syn: u32 = 0;
    let mut nonsyn: u32 = 0;

    for chunk in seq.chunks(3) {
        // reached a point where no complete codon is present
        if chunk.len() < 3 {
            break;
        }

        let (s, n) = get_expected_changes_uc(chunk);

        syn += s;
        nonsyn += n;
    }
    (syn, nonsyn)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_get_seq_exp_changes() {
        let seq = b"ACTGAT";
        let syn_matrix = get_exp_matrix(&UNIVERSAL_CODE);
        let (syn, nonsyn) = get_seq_exp_changes(seq);
        assert_eq!(syn, syn_matrix["ACT"].syn + syn_matrix["GAT"].syn);
        assert_eq!(nonsyn, syn_matrix["ACT"].nonsyn + syn_matrix["GAT"].nonsyn);
    }

    #[test]
    fn test_get_seq_exp_changes_uneven_codons() {
        let seq = b"ACTGATF";
        let syn_matrix = get_exp_matrix(&UNIVERSAL_CODE);
        let (syn, nonsyn) = get_seq_exp_changes(seq);
        assert_eq!(syn, syn_matrix["ACT"].syn + syn_matrix["GAT"].syn);
        assert_eq!(nonsyn, syn_matrix["ACT"].nonsyn + syn_matrix["GAT"].nonsyn);
    }

    #[test]
    fn test_complement() {
        let original = b"ACTGATATATGCGCGCATCTC";
        let reverse_c = b"GAGATGCGCGCATATATCAGT";
        let c = complement(original.iter().rev());
        assert_eq!(reverse_c.to_vec(), c)
    }

    #[test]
    fn test_complement_slice() {
        let original = b"ACTGATATATGCGCGCATCTC";
        let reverse_c = b"GAGATGCGCGCATATATC";
        let c = complement(original[3..].iter().rev());
        assert_eq!(reverse_c.to_vec(), c)
    }

    #[test]
    fn test_record_is_correct() {
        let record = SequenceRecord {
            id: "test".to_string(),
            seq: "ACTG".bytes().collect(),
            ..Default::default()
        };
        assert!(record.is_correct());
    }

    #[test]
    fn test_record_is_correct_fail() {
        let record = SequenceRecord {
            id: "test".to_string(),
            seq: "ACTGX".bytes().collect(),
            ..Default::default()
        };
        assert!(!record.is_correct());
    }

    #[test]
    fn test_universal_code() {
        for (key, value) in UNIVERSAL_CODE.entries() {
            let result = universal_code(key.as_bytes()).unwrap();
            assert_eq!(&std::str::from_utf8(&[result]).unwrap(), value)
        }
    }
}
