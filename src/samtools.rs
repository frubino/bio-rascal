//! Module to deal with files produced by `samtools`, generally used to supplent or avoid the use
//! of other libraries and use the output of `samtools` when possible.
use super::io::open_file;
use anyhow::Result;
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

/// Structure used to old the indices and coverage that are the output of `samtools depth`.
/// This is useful to get the coverage at a specific interval, without having to use a more
/// complex crate to get the read counts and infer the coverage from that.
#[derive(Debug, Default)]
pub struct Depth {
    /// Contains the 1-based indices, second column of the file
    indices: Vec<u32>,
    /// Contains the 1-based coverage, third column of the file
    ///
    /// Using u16 saves about 1/3 of memory in a test
    coverage: Vec<u32>,
}

impl Depth {
    /// Calculates the coverage at the range specified. The indices are 1-based as in `samtools`
    /// and the average coverage is calculated by also adding `0`s for positions not covered.
    /// The return value is the mean as an interger.
    pub fn coverage_at(&self, start: &u32, end: &u32) -> u32 {
        let count = end - start + 1;
        self.indices
            .iter()
            .zip(self.coverage.iter())
            .filter(|x| x.0 >= start && x.0 <= end)
            .map(|x| x.1)
            .sum::<u32>()
            / count as u32
    }
    /// Adds a position and its coverage to the structure
    pub fn add_value(&mut self, pos: u32, cov: u32) {
        self.indices.push(pos);
        self.coverage.push(cov);
    }
    
    pub fn get_coverage_list(&self, start: &u32, end: &u32) -> Vec<(u32, u32)> {
        self.indices
            .iter()
            .zip(self.coverage.iter())
            .filter(|x| x.0 >= start && x.0 <= end)
            .map(|(p, c)| (*p, *c)).collect()
    }
}

/// Alias for an HashMap that contains the sequence Id (first column in `samtools depth` output)
/// and corresponding `Depth`
pub type DepthMap = HashMap<String, Depth>;

/// Reads a gzipped depth file, created using `samtools depth [FILE]`. Accepts the Path of the file
/// and return a `DepthMap` or an Error.
///
/// In case the value of coverage cannot be parsed, it is assumed to go over
/// the boundary for `u16`, so `u16::MAX` is used.
pub fn read_depth_file<P: AsRef<Path>>(file_name: P) -> Result<DepthMap> {
    //info!("Reading depth file: {}", file_name.as_ref().display());
    let file_handle = open_file(file_name.as_ref())?;

    let mut depth_map = DepthMap::default();

    for line in file_handle.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        let seq_id = fields[0].to_string();
        let pos: u32 = fields[1].parse()?;
        // Assumes that the integer is out of bound and put the maximum value
        let cov: u16 = fields[2].parse().unwrap_or(u16::MAX);

        let record = depth_map.entry(seq_id).or_default();
        record.add_value(pos, cov);
    }

    Ok(depth_map)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_depth_add_value() {
        let mut d = Depth::default();
        d.add_value(1, 3);
        assert_eq!(d.indices[0], 1);
        assert_eq!(d.coverage[0], 3);
    }

    #[test]
    fn test_depth_coverage1() {
        let mut d = Depth::default();
        let pos: Vec<u32> = Vec::from_iter(1..=5);
        let cov: Vec<u16> = vec![4; 5];
        for (p, c) in pos.iter().zip(cov.iter()) {
            d.add_value(*p, *c);
        }

        assert_eq!(d.coverage_at(&1, &5), 4);
    }

    #[test]
    fn test_depth_coverage2() {
        let mut d = Depth::default();
        let pos: Vec<u32> = Vec::from_iter(1..=5);
        let cov: Vec<u16> = vec![4; 5];
        for (p, c) in pos.iter().zip(cov.iter()) {
            d.add_value(*p, *c);
        }

        assert_eq!(d.coverage_at(&1, &6), 3);
    }
}
