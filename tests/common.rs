use bio_rascal::taxon::Taxonomy;
use anyhow::Result;
use std::env::var;
use std::path::PathBuf;

pub fn setup_taxonomy() -> Result<Taxonomy> {
    let base_dir = PathBuf::from(var("CARGO_MANIFEST_DIR")?);
    let file_name = base_dir.join( "tests/test-data/taxonomy.json.gz");
    eprintln!("Loading from {:?}", file_name);
    Taxonomy::read_from_file(file_name)
}
