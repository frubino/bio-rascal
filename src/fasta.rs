//! Module used to read and write Fasta files.

use super::sequence::SequenceRecord;
use anyhow::{Context, Result};
use log::info;
use std::fs::File;
use std::io::{BufReader, Read};
use std::{io::BufRead, path::Path};

/// Parses a Fasta header, splits by whitespace and keeps the first
/// element ad ID - This is similar in behaviour as most tools. The
/// rest of the header is also returned, for possible attributes to
/// be parsed by the caller.
fn parse_fasta_seq_header(line: &str) -> (String, String) {
    let header = line[1..].split_once(' ');
    match header {
        None => (line[1..].to_string(), "".to_string()),
        Some((name, attrs)) => (name.to_string(), attrs.to_string()),
    }
}

/// Parses a Fasta file, given a Path to a file.
pub fn parse_fasta_file<P>(file_name: P) -> Result<Vec<SequenceRecord>>
where
    P: AsRef<Path>,
{
    let mut records: Vec<SequenceRecord> = vec![];

    let file_handle = File::open(file_name)
        .map(BufReader::new)
        .context("Cannot open the Fasta file")?;
    let mut name: String = String::new();
    let mut seq: Vec<u8> = vec![];

    for line in file_handle.lines() {
        let line = line?;

        if line.starts_with('>') {
            let new_header = line;

            if !seq.is_empty() && !name.is_empty() {
                let (seq_id, attrs) = parse_fasta_seq_header(&name);
                records.push(SequenceRecord {
                    id: seq_id,
                    seq,
                    attributes: attrs,
                });
                seq = vec![];
            }
            name = new_header;
        } else {
            seq.extend(line.trim().to_uppercase().as_bytes());
        }
    }

    if !seq.is_empty() && !name.is_empty() {
        let (seq_id, attrs) = parse_fasta_seq_header(&name);
        records.push(SequenceRecord {
            id: seq_id,
            seq,
            attributes: attrs,
        });
    }

    Ok(records)
}

/// Fasta file, used to iterate over a file
pub struct FastaReader {
    bufreader: BufReader<Box<dyn Read>>,
    curr_seq: Vec<u8>,
    curr_name: String,
}

impl FastaReader {
    /// A new FastaRead from a file name/path
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<Self> {
        //info!("Opening file {}", file_name.as_ref().display());
        let file_handle = super::io::open_file_base(file_name).context("Cannot open fasta file")?;
        let bufreader = BufReader::new(file_handle);
        Ok(FastaReader {
            bufreader,
            curr_name: "".into(),
            curr_seq: vec![],
        })
    }

    /// TODO: Document
    pub fn from_reader(reader: Box<dyn Read>) -> Self {
        FastaReader {
            bufreader: BufReader::new(reader),
            curr_name: "".into(),
            curr_seq: vec![],
        }
    }

    /// Reads the next Fasta sequence. internally called from `next`
    pub fn read_next_sequence(&mut self) -> Option<SequenceRecord> {
        let mut line = String::new();

        while let Ok(read_len) = self.bufreader.read_line(&mut line) {
            line = line.trim().into();
            //println!("{}, {}: {}, {:?}", line, read_len, self.curr_name, self.curr_seq);
            if read_len == 0 {
                //println!("{} - {:?}", self.curr_name, self.curr_seq);
                if !self.curr_seq.is_empty() && !self.curr_name.is_empty() {
                    let (seq_id, attrs) = parse_fasta_seq_header(&self.curr_name);
                    let seq = self.curr_seq.clone();
                    self.curr_seq.clear();
                    self.curr_name = line.clone();
                    return Some(SequenceRecord {
                        id: seq_id,
                        seq,
                        attributes: attrs,
                    });
                } else {
                    return None;
                }
            }

            if line.starts_with('>') {
                if !self.curr_seq.is_empty() && !self.curr_name.is_empty() {
                    let (seq_id, attrs) = parse_fasta_seq_header(&self.curr_name);
                    let seq = self.curr_seq.clone();
                    self.curr_seq.clear();
                    self.curr_name = line.clone();
                    return Some(SequenceRecord {
                        id: seq_id,
                        seq,
                        attributes: attrs,
                    });
                }
                self.curr_name = line.clone();
            } else {
                self.curr_seq.extend(line.trim().to_uppercase().as_bytes());
            }
            line.clear();
        }

        None
    }
}

impl Iterator for FastaReader {
    type Item = SequenceRecord;
    fn next(&mut self) -> Option<Self::Item> {
        self.read_next_sequence()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::env::var;
    use std::path::PathBuf;

    #[test]
    fn test_header_parse() {
        let header = ">name attr1=value attr2=value";
        let (name, attrs) = parse_fasta_seq_header(header);
        assert_eq!(name, "name");
        assert_eq!(attrs, "attr1=value attr2=value")
    }
    #[test]
    fn test_header_parse1() {
        let header = ">name";
        let (name, attrs) = parse_fasta_seq_header(header);
        assert_eq!(name, "name");
        assert_eq!(attrs, "")
    }

    #[test]
    fn test_fasta_file_attrs() {
        let base_dir = PathBuf::from(var("CARGO_MANIFEST_DIR").unwrap());
        let fasta_file = base_dir.join("src/test-data/test-fasta1.fa");
        let records = parse_fasta_file(fasta_file).unwrap();
        assert_eq!(records[0].id, "name1");
        assert_eq!(records[1].id, "name2");
        assert_eq!(records[2].id, "name3");
        assert_eq!(records[0].seq.len(), 18);
        assert_eq!(records[1].seq.len(), 33);
        assert_eq!(records[2].seq.len(), 12);
        assert_eq!(records[0].seq, b"ACTGATGATGAGACGATA");
        assert_eq!(records[1].seq, b"ACTAGATCAGACATAGCATCAGCATGCATGCAT");
        assert_eq!(records[2].seq, b"ACACATGCATCA");
        assert_eq!(records[0].attributes, "attr1=value");
        assert_eq!(records[1].attributes, "attr2=value");
        assert_eq!(records[2].attributes, "attr3=value3");
    }

    #[test]
    fn test_fastareader() {
        let base_dir = PathBuf::from(var("CARGO_MANIFEST_DIR").unwrap());
        let fasta_file = base_dir.join("src/test-data/test-fasta1.fa");

        let records: Vec<SequenceRecord> = FastaReader::new(&fasta_file).unwrap().collect();
        assert_eq!(records[0].id, "name1");
        assert_eq!(records[1].id, "name2");
        assert_eq!(records[2].id, "name3");
        assert_eq!(records[0].seq.len(), 18);
        assert_eq!(records[1].seq.len(), 33);
        assert_eq!(records[2].seq.len(), 12);
        assert_eq!(records[0].seq, b"ACTGATGATGAGACGATA");
        assert_eq!(records[1].seq, b"ACTAGATCAGACATAGCATCAGCATGCATGCAT");
        assert_eq!(records[2].seq, b"ACACATGCATCA");
        assert_eq!(records[0].attributes, "attr1=value");
        assert_eq!(records[1].attributes, "attr2=value");
        assert_eq!(records[2].attributes, "attr3=value3");
    }
}
