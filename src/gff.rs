//! Module for using annotations and GFF files

use crate::sequence::{get_seq_exp_changes, reverse_comp, universal_code, Sequence};
use crate::taxon::ROOT_TAXON;
use anyhow::{anyhow, bail, Context, Result};
use log::info;
use std::collections::HashMap;
use std::fmt::{self, Display};
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use urlencoding::{decode, encode};
use uuid::Uuid;

/// Strand, can be positive `+` or negative `-`
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Strand {
    Positive,
    Negative,
}

impl Default for Strand {
    /// Assumes Positive as Default
    fn default() -> Strand {
        Strand::Positive
    }
}

impl Strand {
    /// Returns the value
    /// 
    /// Defaults to `Strand::Positive`
    pub fn from_value(field: &str) -> Strand {
        match field {
            "-" => Strand::Negative,
            _ => Strand::Positive,
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Strand::Negative => write!(f, "-"),
            Strand::Positive => write!(f, "+"),
        }
    }
}

/// Unimplemented for now, needs to understand which feature types to use and
/// and if it makes sense to use/keep (Mainly memory, using `Other(String)`)
/// is not a good match, alternative is to use interned strings.
pub enum FeatureType {
    CDS,
    Gene,
    Other(String),
}

/// Phase can be 0, 1, 2 for CDS features or '.' for others, this last
/// one is made equivalent to 0
#[derive(Debug, Clone, Copy)]
pub enum Phase {
    Phase0,
    Phase1,
    Phase2,
    NoPhase,
}

impl Phase {
    // Using a u32 to avoid a cast when adding to Annotation.start
    pub fn to_value(&self) -> u32 {
        match *self {
            Phase::Phase0 | Phase::NoPhase => 0,
            Phase::Phase1 => 1,
            Self::Phase2 => 2,
        }
    }

    /// Parses a Phase field and returns the correct Value
    pub fn from_value(field: &str) -> Result<Phase> {
        match field {
            "0" => Ok(Phase::Phase0),
            "1" => Ok(Phase::Phase1),
            "2" => Ok(Phase::Phase2),
            "." => Ok(Phase::NoPhase),
            _ => Err(anyhow!("Unknown phase {}", field)),
        }
    }
}

impl fmt::Display for Phase {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Phase::NoPhase => write!(f, "."),
            Phase::Phase0 => write!(f, "0"),
            Phase::Phase1 => write!(f, "1"),
            Phase::Phase2 => write!(f, "2"),
        }
    }
}

impl Default for Phase {
    /// NoPhase is the Default
    fn default() -> Self {
        Phase::NoPhase
    }
}

/// Defines an annotation, with fields from a GFF line, plus some information
/// that is usually put in other contexts
#[derive(Debug, Default, Clone)]
pub struct Annotation {
    /// The sequence the annotation is described for
    pub seq_id: String,
    /// Source of the annotation, usually the software used
    pub source: String,
    /// Can be of several types, it's a string. CDS, gene, exon are all used
    pub feature_type: String,
    /// 1-based start of the annotation
    pub start: u32,
    /// End of the annotation, inclusive
    pub end: u32,
    /// Usually the score of the software.
    pub score: f64,
    /// The strand on which the annotation is defined
    pub strand: Strand,
    /// `Phase`
    pub phase: Phase,
    /// Taxon ID, part of the attributes
    pub taxon_id: u32,
    /// Unique ID for annotation, if it is not found, one will be created. Part of the attributes
    pub uid: Uuid,
    /// The rest of attributes on the last column of the file
    pub attributes: HashMap<String, String>,
}

impl Annotation {
    /// Returns the value of an attribute, useful for cases
    /// where `taxon_id` or `uid` are requested.
    pub fn get_attr<S>(&self, key: S) -> Option<String>
    where
        S: AsRef<str>,
    {
        let key = key.as_ref();
        match key {
            "taxon_id" => match self.taxon_id {
                ROOT_TAXON.. => Some(self.taxon_id.to_string()),
                _ => None,
            },
            "uid" => Some(self.uid.to_string()),
            _ => self.attributes.get(key).cloned(),
        }
    }
    /// Returns the length of the annotation
    pub fn length(&self) -> u32 {
        self.end - self.start + 1
    }
    /// Returns the region as used in `samtools` and others
    pub fn get_region(&self) -> String {
        format!("{}:{}-{}", self.seq_id, self.start, self.end)
    }

    pub fn get_seq(&self, seq: &Sequence) -> Sequence {
        let start = (self.start - 1) as usize;
        let end = self.end as usize;
        if self.strand == Strand::Positive {
            seq[start..end].iter().copied().collect()
        } else {
            reverse_comp(&seq[start..end])
        }
    }

    pub fn get_exp_syn(&self, seq: &Sequence) -> (u32, u32) {
        get_seq_exp_changes(&self.get_seq(&seq))
    }

    pub fn new<L: AsRef<str>>(line: L) -> Result<Self> {
        parse_gff_line(line, true)
    }

    /// Given a position on the chromosome, returns the position relative
    /// to the Annotation start
    /// 
    /// A value of 0 is not accepted and will return the same error as
    /// a position outside a valid position (but outside the annotation)
    pub fn get_relative_pos(&self, pos: u32) -> Result<u32> {
        if pos > self.end || pos < 1 {
            bail!("Position outside the annotation")
        } else {
            Ok(pos - self.start + 1)
        }
    }

    /// Returns `true` if the position is within the `Annotation`
    pub fn contains(&self, pos: u32) -> bool {
        pos >= self.start && pos <= self.end
    }
    
    /// This is useful to get the codon at a specific position
    /// TODO: check against: 7a04c052-e562-4891-9a98-425abcd36b70 sample H_3
    pub fn get_codon_at<'a>(&self, pos: u32, seq: &'a [u8]) -> Result<&'a [u8]> {
        let rel_pos = self.get_relative_pos(pos).context("Cannot get relative position")?;
        
        let codon_index = (rel_pos + self.phase.to_value() - 1) / 3;
        
        let start: usize = (self.start - 1 + self.phase.to_value() + (codon_index * 3)) as usize;
        let end = start + 3;
        if end > seq.len() {
            bail!("Position is out of boundaries")
        }

        Ok(&seq[start..end])
    }

    /// Returns if a nucleotide change is synonymous or non-synonymous
    ///
    /// Requires the chromosome, change position and change. An error is
    /// returned if the position is out of the `Annotation` or the codons
    /// are not present in the `universal_code`.
    pub fn is_syn(&self, seq: &[u8], pos: u32, change: &u8) -> Result<bool> {
        // First find the relative position
        let rel_pos = self.get_relative_pos(pos).context("Cannot get relative position")?;
        // position in the codon using the relative position and the phase/frame
        let codon_change: usize = ((rel_pos + self.phase.to_value() - 1) % 3) as usize;
        // Builds a vector, as we need to clone it and substitute the value,
        // moreover the reverse/complement size needs to be know because of
        // the invert iterator
        let mut codon = self.get_codon_at(pos, seq).context("Cannot get codon")?.to_vec();
        let mut var_codon = codon.clone();
        var_codon[codon_change] = *change;

        // Reverse complement needed
        if self.strand == Strand::Negative {
            codon = reverse_comp(&codon);
            var_codon = reverse_comp(&var_codon);
        }

        let codon_aa = universal_code(&codon)?;
        let var_codon_aa = universal_code(&var_codon)?;

        Ok(codon_aa == var_codon_aa)
    }
}

impl Display for Annotation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut attributes = String::with_capacity(150);

        if self.taxon_id >= ROOT_TAXON {
            attributes.push_str(";taxon_id=\"");
            attributes.push_str(&self.taxon_id.to_string());
            attributes.push('"');
        }

        for (key, value) in self.attributes.iter() {
            attributes.push(';');
            attributes.push_str(key);
            attributes.push_str("=\"");
            attributes.push_str(&encode(value));
            attributes.push('"');
        }
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tuid=\"{}\"{}",
            self.seq_id,
            self.source,
            self.feature_type,
            self.start,
            self.end,
            self.score,
            self.strand,
            self.phase,
            self.uid,
            attributes
        )
    }
}

/// Structure to implement an `Iterator` over a GFF file
pub struct GffReader {
    bufreader: BufReader<Box<dyn Read>>,
}

impl GffReader {
    /// Creates a GffReader from a path
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<Self> {
        info!("Opening file {}", file_name.as_ref().display());
        let bufreader = super::io::open_file_base(file_name)?;
        Ok(GffReader {
            bufreader: BufReader::new(bufreader),
        })
    }

    /// Creates a GffReader from a BufReader
    pub fn from_reader(reader: Box<dyn Read>) -> Self {
        GffReader {
            bufreader: BufReader::new(reader),
        }
    }
}

impl Iterator for GffReader {
    type Item = Annotation;

    /// Iterates over the internal GFF file buffer untile the end of the file
    /// or a Fasta sequence is found
    fn next(&mut self) -> Option<Self::Item> {
        let mut buffer = String::new();

        if let Ok(read_size) = self.bufreader.read_line(&mut buffer) {
            // end of file, stop reading
            if read_size == 0 {
                None
            } else if buffer.starts_with('#') {
                // Goes to next iteration
                self.next()
            // starts of sequence, stop reading
            } else if buffer.starts_with('>') {
                None
            } else {
                match Annotation::new(buffer) {
                    Err(_) => None,
                    Ok(annotation) => Some(annotation),
                }
            }
        // Error, stop reading
        } else {
            None
        }
    }
}

/// Parses a GFF line, corresponding to an annotation. Attributes can be skipped
/// by passing `false` to `include_attributes`.
pub fn parse_gff_line<L: AsRef<str>>(line: L, include_attributes: bool) -> Result<Annotation> {
    // Trims the line, splits and collect the fields separated by `\t`
    let fields: Vec<&str> = line.as_ref().trim().split('\t').collect();

    let mut annotation: Annotation = Annotation {
        seq_id: fields[0].to_string(),
        source: fields[1].to_string(),
        feature_type: fields[2].to_string(),
        start: fields[3].parse::<u32>().context("Parsing Start")?,
        end: fields[4].parse::<u32>().context("Parsing End")?,
        score: fields[5].parse::<f64>().unwrap_or(0.),
        strand: Strand::from_value(fields[6]),
        phase: Phase::from_value(fields[7]).context("Parsing Phase")?,
        ..Default::default()
    };

    for attribute in fields[8].split(';') {
        let v: Vec<&str> = attribute
            .split('=')
            .map(|el| el.trim_matches('"'))
            .collect();
        let key: &str = v[0];
        let value = decode(v[1]).unwrap();
        match key {
            "taxon_id" => {
                let taxon_id = value.parse::<i32>().context("Parsing Taxon ID")?;
                if taxon_id < 0 {
                    // Cases where the taxon_id is negative
                    annotation.taxon_id = (u32::MAX as i32 + taxon_id) as u32;
                } else {
                    annotation.taxon_id = taxon_id as u32;
                }
            }
            "uid" => annotation.uid = Uuid::try_parse(value.as_ref())?,
            _ => {
                if include_attributes {
                    annotation
                        .attributes
                        .insert(key.to_string(), value.to_string());
                }
            }
        }
    }
    // if no `uid` is defined for the annotation, creates one
    if annotation.uid.is_nil() {
        annotation.uid = Uuid::new_v4();
    }

    Ok(annotation)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::env::var;
    use std::fs::File;
    use std::path::PathBuf;

    #[test]
    fn test_parse_gff_line() {
        let test_line = include_str!("test-data/test-line.gff");
        // needs to use trim, because .lines() removes them, but here is included
        let annotation = parse_gff_line(test_line.trim().to_string(), true).unwrap();
        assert_eq!(annotation.taxon_id, 1263104);
        assert_eq!(
            annotation.uid.to_string(),
            "465473df-3447-4f17-b84a-d575935929b7"
        );
        assert_eq!(
            annotation
                .get_attr("identity")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            53.85
        );
    }

    #[test]
    fn test_parse_gff_line_taxon_id_negative() {
        let test_line = include_str!("test-data/test-line.gff").replace("1263104", "-1");
        // needs to use trim, because .lines() removes them, but here is included
        let annotation = parse_gff_line(test_line.trim().to_string(), true).unwrap();
        assert_eq!(annotation.taxon_id, u32::MAX - 1);
        assert_eq!(
            annotation.uid.to_string(),
            "465473df-3447-4f17-b84a-d575935929b7"
        );
        assert_eq!(
            annotation
                .get_attr("identity")
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            53.85
        );
    }

    #[test]
    fn test_get_attr() {
        let mut annotation: Annotation = Annotation {
            ..Default::default()
        };
        annotation
            .attributes
            .insert("test".to_string(), "test".to_string());
        assert_eq!(annotation.get_attr("test").unwrap(), "test");
    }

    #[test]
    fn test_gffreader() {
        let base_dir = PathBuf::from(var("CARGO_MANIFEST_DIR").unwrap());
        let file_name = base_dir.join("src/test-data/test-line2.gff");
        let reader = GffReader::new(&file_name).unwrap();

        let annotations: Vec<Annotation> = reader.collect();
        assert_eq!(
            annotations[0].uid.to_string(),
            "465473df-3447-4f17-b84a-d575935929b7"
        );
        assert_eq!(
            annotations[1].uid.to_string(),
            "465473df-3447-4f17-b84a-d575935929b1"
        );
        assert_eq!(annotations.len(), 2);
    }

    #[test]
    fn test_gffreader_from_reader() {
        let base_dir = PathBuf::from(var("CARGO_MANIFEST_DIR").unwrap());
        let file_name = base_dir.join("src/test-data/test-line2.gff");
        let reader =
            GffReader::from_reader(Box::new(BufReader::new(File::open(&file_name).unwrap())));

        let annotations: Vec<Annotation> = reader.collect();
        assert_eq!(
            annotations[0].uid.to_string(),
            "465473df-3447-4f17-b84a-d575935929b7"
        );
        assert_eq!(
            annotations[1].uid.to_string(),
            "465473df-3447-4f17-b84a-d575935929b1"
        );
        assert_eq!(annotations.len(), 2);
    }
    #[test]
    fn test_annotation_get_relative_pos() {
        let mut a = Annotation::default();
        a.start = 2;
        a.end = 10;
        assert_eq!(a.get_relative_pos(3).unwrap(), 2);
        assert!(a.get_relative_pos(12).is_err());
        // not effected by phase
        a.phase = Phase::Phase1;
        assert_eq!(a.get_relative_pos(3).unwrap(), 2);
    }

    #[test]
    fn test_annotation_is_syn() {
        let seq = b"ATTAATG";
        let mut annotation = Annotation::default();
        annotation.start = 1;
        annotation.end = 7;
        annotation.phase = Phase::Phase0;
        assert!(!annotation.is_syn(seq, 1, &b'T').unwrap());
        assert!(annotation.is_syn(seq, 1, &b'A').unwrap());
        assert!(!annotation.is_syn(seq, 2, &b'G').unwrap());
        assert!(!annotation.is_syn(seq, 4, &b'G').unwrap());
        annotation.phase = Phase::Phase1;
        assert!(!annotation.is_syn(seq, 1, &b'G').unwrap());
    }

    #[test]
    fn test_annotation_contains() {
        let annotation = Annotation {
            start: 2,
            end: 10,
            ..Default::default()
        };
        assert!(annotation.contains(2));
        assert!(annotation.contains(10));
        assert!(!annotation.contains(11));
        assert!(!annotation.contains(1));
    }
    
    #[test]
    fn test_get_codon_at() {
        let seq = b"TTATTGTTATTG";
        let mut a = Annotation::default();
        a.start = 1;
        a.end = 12;
        assert_eq!(a.get_codon_at(1, seq).unwrap(), b"TTA");
        assert_eq!(a.get_codon_at(2, seq).unwrap(), b"TTA");
        assert_eq!(a.get_codon_at(3, seq).unwrap(), b"TTA");
        assert_eq!(a.get_codon_at(4, seq).unwrap(), b"TTG");
        assert_eq!(a.get_codon_at(5, seq).unwrap(), b"TTG");
        assert_eq!(a.get_codon_at(6, seq).unwrap(), b"TTG");
        assert_eq!(a.get_codon_at(7, seq).unwrap(), b"TTA");
        assert_eq!(a.get_codon_at(8, seq).unwrap(), b"TTA");
        assert_eq!(a.get_codon_at(9, seq).unwrap(), b"TTA");
        assert_eq!(a.get_codon_at(10, seq).unwrap(), b"TTG");
        assert_eq!(a.get_codon_at(11, seq).unwrap(), b"TTG");
        assert_eq!(a.get_codon_at(12, seq).unwrap(), b"TTG");
        // check boundaries
        assert!(a.get_codon_at(0, seq).is_err());
        assert!(a.get_codon_at(13, seq).is_err());
        // checks that it respects the Phase
        a.phase = Phase::Phase1;
        assert_eq!(a.get_codon_at(1, seq).unwrap(), b"TAT");
        // checks out of boundaries
        a.phase = Phase::Phase0;
        let seq = b"TTATTGTTATT";
        assert_eq!(a.get_codon_at(9, seq).unwrap(), b"TTA");
        assert!(a.get_codon_at(10, seq).is_err());
    }
}
