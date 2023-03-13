//! Module for SNPs and VCF files
//!
use super::sequence::Sequence;
use super::taxon::ROOT_TAXON;
use anyhow::{bail, Context, Result};
use log::error;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use uuid::Uuid;

pub fn get_expected_changes_uc(codon: &[u8]) -> (u32, u32) {
    match codon {
        b"ACA" | b"ACC" | b"ACG" | b"ACT" | b"CCA" | b"CCC" | b"CCG" | b"CCT" | b"CGC" | b"CGT"
        | b"CTC" | b"CTT" | b"GCA" | b"GCC" | b"GCG" | b"GCT" | b"GGA" | b"GGC" | b"GGG"
        | b"GGT" | b"GTA" | b"GTC" | b"GTG" | b"GTT" | b"TCA" | b"TCC" | b"TCG" | b"TCT" => (3, 6),
        b"AAA" | b"AAC" | b"AAG" | b"AAT" | b"AGC" | b"AGT" | b"CAA" | b"CAC" | b"CAG" | b"CAT"
        | b"GAA" | b"GAC" | b"GAG" | b"GAT" | b"TAC" | b"TAG" | b"TAT" | b"TGA" | b"TGC"
        | b"TGT" | b"TTC" | b"TTT" => (1, 8),
        b"AGA" | b"AGG" | b"ATA" | b"ATC" | b"ATT" | b"TAA" | b"TTA" | b"TTG" => (2, 7),
        b"CGA" | b"CGG" | b"CTA" | b"CTG" => (4, 5), //ok
        b"ATG" | b"TGG" => (0, 9),
        _ => {
            error!("Cannot find codon {}", String::from_utf8_lossy(codon));
            (0, 0)
        },
    }
}

pub trait CalculatePnPs {
    fn get_pn(&self) -> f64;
    fn get_ps(&self) -> f64;
    fn get_pnps(&self) -> f64 {
        let ps = self.get_ps();
        let pn = self.get_pn();
        pn / ps
    }
}

#[derive(Debug, Default, Serialize, Deserialize, Clone, Copy)]
pub struct PnPs {
    pub uid: Uuid,
    pub coverage: u32,
    pub syn: u32,
    pub nonsyn: u32,
    pub exp_syn: u32,
    pub exp_nonsyn: u32,
}

impl CalculatePnPs for PnPs {
    fn get_ps(&self) -> f64 {
        if self.syn == 0 || self.exp_syn == 0 {
            f64::NAN
        } else {
            self.syn as f64 / self.exp_syn as f64
        }
    }
    fn get_pn(&self) -> f64 {
        if self.nonsyn == 0 || self.exp_nonsyn == 0 {
            f64::NAN
        } else {
            self.nonsyn as f64 / self.exp_nonsyn as f64
        }
    }
}

/// Structure to group [PnPs] values using keys
#[derive(Debug, Default)]
pub struct GroupPnPs<'a> {
    /// Gene ID
    pub gene_id: String,
    /// Taxon ID, refering to a taxonomy
    pub taxon_id: u32,
    /// Alternative to [GroupPnPs::taxon_id]
    pub taxon_lineage: String,
    /// Vector of [PnPs]
    pub pnps: Vec<&'a PnPs>,
}

impl<'a> GroupPnPs<'a> {
    pub fn new() -> Self {
        GroupPnPs {
            taxon_id: ROOT_TAXON,
            ..Default::default()
        }
    }
}

impl<'a> CalculatePnPs for GroupPnPs<'a> {
    fn get_ps(&self) -> f64 {
        let syn = self.pnps.iter().fold(0u32, |acc, p| acc + p.syn);
        let exp_syn = self.pnps.iter().fold(0u32, |acc, p| acc + p.exp_syn);
        if syn == 0 || exp_syn == 0 {
            f64::NAN
        } else {
            syn as f64 / exp_syn as f64
        }
    }

    fn get_pn(&self) -> f64 {
        let nonsyn = self.pnps.iter().fold(0u32, |acc, p| acc + p.nonsyn);
        let exp_nonsyn = self.pnps.iter().fold(0u32, |acc, p| acc + p.exp_nonsyn);
        if nonsyn == 0 || exp_nonsyn == 0 {
            f64::NAN
        } else {
            nonsyn as f64 / exp_nonsyn as f64
        }
    }
}

/// Stores information about the `INFO` section of a [VcfRecord]
#[derive(Debug, Default)]
pub struct VcfRecordInfo {
    /// Indicates that the variant is an INDEL.
    pub indel: bool,
    /// Maximum number of raw reads supporting an indel
    pub idv: u32,
    /// Maximum fraction of raw reads supporting an indel
    pub imf: f64,
    /// Raw read depth
    pub dp: u32,
    /// Allele count in genotypes, for each ALT allele, in the same order as listed
    pub ac: Vec<u32>,
    /// total number of alleles in called genotypes
    pub an: u32,
    pub mq: u32,
    pub others: HashMap<String, String>,
}

impl VcfRecordInfo {
    pub fn from_string<F: AsRef<str>>(line: F) -> Result<Self> {
        let mut vri = VcfRecordInfo::default();

        for field in line.as_ref().split(';') {
            match field {
                // since default for bool is `false`
                // this seems to be a specific value, out of standard for bcftools
                // since it doesn't have a value, we keep it this way
                "INDEL" => vri.indel = true,
                _ => {
                    let (tag, value) = field.split_once('=').unwrap_or(("", ""));
                    match tag {
                        "" => continue,
                        "IDV" => vri.idv = value.parse().context("Failed to parse `IDV`")?,
                        "IMF" => vri.imf = value.parse().context("Failed to parse `IMF`")?,
                        "DP" => vri.dp = value.parse().context("Failed to parse `DP`")?,
                        "AC" => {
                            vri.ac = value.split(',').map(|e| e.parse().ok()).flatten().collect()
                        }
                        "AN" => vri.an = value.parse().context("Failed to parse `AN`")?,
                        "MQ" => vri.mq = value.parse().context("Failed to parse `MQ`")?,
                        _ => _ = vri.others.insert(tag.to_string(), value.to_string()),
                    }
                }
            }
        }

        Ok(vri)
    }
}

/// Last column in a [VcfRecord] and is a per-sample field
#[derive(Debug, Default)]
pub struct GenotypeFields {
    pub sample_name: String,
    pub gt: u8,
    pub others: HashMap<String, String>,
}

impl GenotypeFields {
    /// Parses the fields of a sample
    pub fn from_string<F: AsRef<str>>(info: F, sample: F, sample_name: F) -> Result<Self> {
        let mut g = GenotypeFields::default();
        g.sample_name = sample_name.as_ref().to_string();
        for (i, v) in info.as_ref().split(':').zip(sample.as_ref().split(':')) {
            match i {
                "GT" => {
                    g.gt = match v {
                        // cannot resolve or its the REF nucleotide
                        "." | "0" => 0,
                        _ => v.parse().context("Failed to parse GT")?,
                    }
                }
                _ => _ = g.others.insert(i.to_string(), v.to_string()),
            }
        }
        Ok(g)
    }
}

/// Implements a VCF file records, use [VcfReader] to get them
#[derive(Debug, Default)]
pub struct VcfRecord {
    /// Chromosome
    pub chrom: String,
    /// 1-based position on [VcfRecord::chrom]
    pub pos: u32,
    /// Unused
    pub id: String,
    /// Nucleotide on the Reference 
    pub ref_c: Sequence,
    /// Alternative (SNPs) to the Reference 
    pub alt_c: Vec<Sequence>,
    /// Quality of the SNP
    pub qual: f64,
    /// If the filter is passed
    pub filt: String,
    pub info: VcfRecordInfo,
    pub samples: Vec<GenotypeFields>,
}

impl VcfRecord {
    pub fn parse_vcf_line<L: AsRef<str>>(line: L, sample_names: &Vec<String>) -> Result<Self> {
        let fields: Vec<&str> = line.as_ref().trim().split('\t').collect();
        let chrom = fields[0].to_string();
        let pos: u32 = fields[1].parse().context("Failed to parse `POS`")?;
        let id = fields[2].to_string();
        let ref_c: Sequence = fields[3].as_bytes().to_vec();
        let alt_c: Vec<Sequence> = fields[4]
            .split(',')
            .map(|f| f.as_bytes().to_vec())
            .collect();
        let qual: f64 = fields[5].parse().context("Failed to parse `QUAL`")?;
        let filt = fields[6].to_string();
        let info: VcfRecordInfo = VcfRecordInfo::from_string(fields[7])?;
        let mut samples: Vec<GenotypeFields> = vec![];
        for (field, sample_name) in fields[9..].iter().zip(sample_names.iter()) {
            samples.push(GenotypeFields::from_string(fields[8], field, sample_name)?);
        }

        Ok(VcfRecord {
            chrom,
            pos,
            id,
            ref_c,
            alt_c,
            qual,
            filt,
            info,
            samples,
        })
    }
    pub fn get_sample_snps(&self) -> Vec<(String, u8)> {
        self.samples
            .iter()
            // zips the names (columns) and fields
            //.zip(self.samples.iter())
            // only keeps thos that have a SNP defined (>0)
            .filter(|g| g.gt > 0)
            // maps the index to the ALT character, but we assume is a SNP, INDEL WILL cause problems
            // TODO: checks or find way to control this situation
            // also, the index needs to be substracted, maybe work around another way?
            .map(|g| {
                (
                    g.sample_name.clone(),
                    self.alt_c[(g.gt - 1) as usize][0].clone(),
                )
            })
            .collect()
    }
}

pub struct VcfReader {
    bufreader: BufReader<Box<dyn Read>>,
    pub file_format: String,
    pub sample_names: Vec<String>,
}

impl VcfReader {
    pub fn new<P: AsRef<Path>>(file_name: P) -> Result<Self> {
        let reader = super::io::open_file_base(file_name)?;
        VcfReader::from_reader(reader)
    }
    pub fn from_reader(reader: Box<dyn Read>) -> Result<Self> {
        let bufreader = BufReader::new(reader);
        let mut s = VcfReader {
            bufreader,
            file_format: String::new(),
            sample_names: vec![],
        };
        s.read_header()?;
        Ok(s)
    }

    fn read_header(&mut self) -> Result<()> {
        // reads the first line, which is the file format
        let mut buffer = String::new();
        self.bufreader.read_line(&mut buffer)?;
        if buffer.starts_with("##fileformat=") {
            self.file_format = buffer.trim()[13..].to_string();
        } else {
            bail!("Unexpected file format line {}", buffer)
        }

        loop {
            buffer.clear();
            self.bufreader.read_line(&mut buffer)?;
            if buffer.starts_with("##") {
                continue;
            } else if buffer.starts_with("#CHROM") {
                self.sample_names = buffer
                    .trim()
                    .split('\t')
                    .skip(9)
                    .map(|f| f.to_string())
                    .collect();
                break;
            } else {
                bail!("Unexpected line {}", buffer)
            }
        }
        Ok(())
    }
}

impl Iterator for VcfReader {
    type Item = VcfRecord;

    fn next(&mut self) -> Option<Self::Item> {
        let mut buffer = String::with_capacity(130);
        match self.bufreader.read_line(&mut buffer) {
            Err(_) | Ok(0) => None,
            Ok(_) => match VcfRecord::parse_vcf_line(&buffer, &self.sample_names) {
                Err(err) => {
                    error!("Error {}: {}", err, &buffer);
                    None
                }
                Ok(value) => Some(value),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::sequence::{get_exp_matrix, UNIVERSAL_CODE};

    use super::*;

    #[test]
    fn test_pnps_none() {
        let p = PnPs::default();
        assert!(p.get_ps().is_nan());
        assert!(p.get_pn().is_nan());
        assert!(p.get_pnps().is_nan());
    }

    #[test]
    fn test_pnps() {
        let p = PnPs {
            syn: 2,
            exp_syn: 2,
            nonsyn: 4,
            exp_nonsyn: 4,
            ..Default::default()
        };
        assert_eq!(p.get_pn(), 1.);
        assert_eq!(p.get_ps(), 1.);
        assert_eq!(p.get_pnps(), 1.);
    }

    #[test]
    fn test_group_pnps() {
        let pnps1 = PnPs {
            syn: 2,
            exp_syn: 4,
            nonsyn: 2,
            exp_nonsyn: 1,
            ..Default::default()
        };
        let pnps2 = PnPs {
            syn: 2,
            exp_syn: 0,
            nonsyn: 0,
            exp_nonsyn: 1,
            ..Default::default()
        };
        let gp = GroupPnPs {
            pnps: vec![&pnps1, &pnps2],
            ..Default::default()
        };
        assert_eq!(gp.get_ps(), 1.);
        assert_eq!(gp.get_pn(), 1.);
        assert_eq!(gp.get_pnps(), 1.);
    }

    #[test]
    fn test_group_pnps_none() {
        let pnps1 = PnPs {
            syn: 2,
            exp_syn: 0,
            nonsyn: 2,
            exp_nonsyn: 1,
            ..Default::default()
        };
        let pnps2 = PnPs {
            syn: 2,
            exp_syn: 0,
            nonsyn: 0,
            exp_nonsyn: 1,
            ..Default::default()
        };
        let gp = GroupPnPs {
            pnps: vec![&pnps1, &pnps2],
            ..Default::default()
        };
        assert!(gp.get_ps().is_nan());
        assert_eq!(gp.get_pn(), 1.);
        assert!(gp.get_pnps().is_nan());
    }

    #[test]
    fn test_group_pnps_new() {
        let gp = GroupPnPs::new();
        assert_eq!(gp.taxon_id, ROOT_TAXON);
    }

    #[test]
    fn test_vcf_record_parse_vcf_line() {
        let test_line = include_str!("test-data/test-snps-line.vcf");
        let sample_names = vec!["sample_name".to_string(); 18];
        let record = VcfRecord::parse_vcf_line(test_line, &sample_names).unwrap();
        assert_eq!(
            record.chrom,
            "Frog-Gut_megahit_1000-data_no-controls-refine-bin.10.fa|k141_10668202|rEe2x"
        );
        assert_eq!(record.pos, 3867);
        assert_eq!(record.qual, 5.4431);
        assert_eq!(record.info.dp, 18);
        assert_eq!(record.info.ac, [1]);
        assert_eq!(record.info.an, 3);
        assert_eq!(record.info.mq, 20);
        assert_eq!(record.samples[1].gt, 0);
        assert_eq!(record.samples[2].gt, 0);
        assert_eq!(record.samples[17].gt, 1);
        assert_eq!(record.samples[17].sample_name, "sample_name");
    }

    #[test]
    fn test_universal_code_expected() {
        let syn_matrix = get_exp_matrix(&UNIVERSAL_CODE);

        for (codon, aa_changes) in syn_matrix.iter() {
            let ac = (aa_changes.syn, aa_changes.nonsyn);
            assert_eq!(get_expected_changes_uc(codon.as_bytes()), ac);
        }
    }
}
