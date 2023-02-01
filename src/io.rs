//! Helper I/O functions
//!

use anyhow::Result;
use flate2::read::GzDecoder;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Checks if the file extension is `gz` and if true, wraps the
/// file into a `GzDecoder` before passing it to `BufReader`.
pub fn open_file<P: AsRef<Path>>(file_name: P) -> Result<Box<dyn BufRead>> {
    let file_name = file_name.as_ref();
    let gz_ext = OsStr::new("gz");
    let file = File::open(file_name)?;
    if file_name
        .extension()
        .filter(|ext| ext.eq(&gz_ext))
        .is_some()
    {
        Ok(Box::new(BufReader::new(GzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Checks if the file extension is `gz` and if true, wraps the
/// file into a `GzDecoder` before returning a boxed `Read`.
pub fn open_file_base<P: AsRef<Path>>(file_name: P) -> Result<Box<dyn Read>> {
    let file_name = file_name.as_ref();
    let gz_ext = OsStr::new("gz");
    let file = File::open(file_name)?;
    if file_name
        .extension()
        .filter(|ext| ext.eq(&gz_ext))
        .is_some()
    {
        Ok(Box::new(GzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}
