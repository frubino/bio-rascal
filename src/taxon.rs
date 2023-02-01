//! Module for Taxonomy
//!
use anyhow::{anyhow, Result};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use log::{error, info};
use serde::{Deserialize, Serialize};
use serde_json::{from_reader, to_writer};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{prelude::*, BufReader, BufWriter};
use std::path::Path;
use time::{Date, OffsetDateTime};

/// Taxonomic Rank, in some taxonomies superkingdom is equivalent to kingdom
#[derive(Debug, PartialEq, Eq, Serialize, Deserialize, PartialOrd, Ord)]
pub enum Rank {
    /// Specific to NCBI
    SuperKingdom,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    Strain,
    /// Specific to NCBI
    NoRank,
    /// Used when no other Rank fits
    Other,
}

impl Default for Rank {
    /// The default Rank is NoRank
    fn default() -> Self {
        Rank::NoRank
    }
}

impl Rank {
    /// Returns the one letter prefix associated with the Rank
    ///
    /// `NoRank` and `Other` become `x`
    pub fn to_prefix(&self) -> &str {
        match self {
            Rank::Class => "c",
            Rank::Family => "f",
            Rank::Genus => "g",
            Rank::SuperKingdom | Rank::Kingdom => "k",
            Rank::Order => "o",
            Rank::Phylum => "p",
            Rank::Species => "s",
            Rank::Strain => "t",
            _ => "x",
        }
    }
    /// Uses the character used in lineages to get the associated Rank
    pub fn from_prefix(rank_char: &str) -> Rank {
        match rank_char {
            "c" => Rank::Class,
            "f" => Rank::Family,
            "g" => Rank::Genus,
            "k" => Rank::Kingdom,
            "o" => Rank::Order,
            "p" => Rank::Phylum,
            "s" => Rank::Species,
            "t" => Rank::Strain,
            _ => Self::NoRank,
        }
    }
    /// Gets the Rank from the name
    ///
    /// These are from NCBI Taxonomy, but subfamilies and similar
    /// become `Other`
    pub fn from_name(rank: &str) -> Rank {
        match rank {
            "superkingdom" | "SuperKingdom" => Rank::SuperKingdom,
            "kingdom" | "Kingdom" => Rank::Kingdom,
            "phylum" | "Phylum" => Rank::Phylum,
            "class" | "Class" => Rank::Class,
            "order" | "Order" => Rank::Order,
            "family" | "Family" => Rank::Family,
            "genus" | "Genus" => Rank::Genus,
            "species" | "Species" => Rank::Species,
            "strain" | "Strain" => Rank::Strain,
            "no rank" | "NoRank" => Rank::NoRank,
            _ => Rank::Other,
        }
    }
    /// Returns the name associated with the Taxonomy
    pub fn to_name(&self) -> &str {
        match self {
            Rank::Class => "class",
            Rank::Family => "family",
            Rank::Genus => "genus",
            Rank::Kingdom => "kingdom",
            Rank::NoRank => "no rank",
            Rank::Order => "order",
            Rank::Phylum => "phylum",
            Rank::Other => "other",
            Rank::Species => "species",
            Rank::Strain => "strain",
            Rank::SuperKingdom => "superkingdom",
        }
    }
}

/// Constant for the scientific name value when parsing NCBI
/// Taxonomy
const SCIENTIFICNAME: &str = "scientific name";
/// Column value for the synonyms in the names (NCBI)
const SYNONYM: &str = "synonym";
/// Root taxon, derived by NCBI
pub const ROOT_TAXON: u32 = 1;

/// Definition of a Taxon, derived from NCBI Taxonomy and uses an unsigned integer for the identifier.
/// The parent ID of a Taxon is kept as well as the scientific name `name` and eventual synonymns for
/// the Taxon.
#[derive(Debug, Default, Serialize, Deserialize)]
pub struct Taxon {
    /// ID of the Taxon
    pub id: u32,
    /// ID of the parent Taxon
    pub parent_id: u32,
    /// Rank of the Taxon
    pub rank: Rank,
    // Scientific Name
    pub name: String,
    // Synonyms
    pub synonyms: Vec<String>,
}

impl Taxon {
    /// Generates a Taxon from a line of the `nodes.dmp` file.
    /// This won't fill the names, since they are in another file
    pub fn parse_ncbi_node(line: &str) -> Result<Self> {
        let mut taxon = Taxon {
            ..Default::default()
        };

        for (index, field) in line.split("\t|\t").map(|el| el.trim()).take(3).enumerate() {
            match index {
                0 => taxon.id = field.parse()?,
                1 => taxon.parent_id = field.parse()?,
                2 => taxon.rank = Rank::from_name(field),
                _ => break,
            }
        }
        Ok(taxon)
    }

    /// Returns the string that is used to get the lineage,
    pub fn to_lineage(&self) -> String {
        format!("{}__{}", self.rank.to_prefix(), self.name.replace(' ', "_"))
    }
}

/// Taxonomy structure, keeps all Taxa and provides methods to access it
#[derive(Debug, Serialize, Deserialize)]
pub struct Taxonomy {
    pub creation_date: Date,
    // Map of TaxonId and the relative `Taxon` structure
    nodes: HashMap<u32, Taxon>,
    // Map of TaxonId that where merged in another Taxon
    merged_nodes: HashMap<u32, u32>,
}

impl Default for Taxonomy {
    /// The default, sets the data structure and uses today as `creation_date`.
    fn default() -> Self {
        Taxonomy {
            creation_date: OffsetDateTime::now_utc().date(),
            nodes: HashMap::new(),
            merged_nodes: HashMap::new(),
        }
    }
}

impl Taxonomy {
    pub fn get_ranked_taxon(&self, taxon_id: u32, rank: Rank) -> Option<&Taxon> {
        // if the taxon_id is not found, return None
        let mut taxon = match self.get_taxon(&taxon_id) {
            None => return None,
            Some(taxon) => {
                if taxon.rank <= rank {
                    return Some(taxon)
                }
                taxon
            },
        };
        loop {
            taxon = match self.get_taxon(&taxon.parent_id){
                None => break,
                Some(taxon) => taxon,
            };
            if taxon.rank <= rank {
                return Some(taxon)
            } else if taxon.id == ROOT_TAXON {
                break;
            }
        }

        None
    }

    /// Inserts a `Taxon` into the taxonomy
    pub fn inset_taxon(&mut self, taxon: Taxon) {
        self.nodes.insert(taxon.id, taxon);
    }
    /// Parses the NCBI merged taxa from `merged.dmp`, it needs to be run
    /// after the taxonomy has been initialised using `parse_ncbi_nodes.
    fn parse_ncbi_merged<P: AsRef<Path>>(&mut self, merged_file: P) -> Result<()> {
        let reader = BufReader::new(File::open(merged_file)?);
        for line in reader.lines() {
            let (old_id, new_id) = parse_ncbi_merged_line(&line.unwrap())?;
            self.merged_nodes.insert(old_id, new_id);
        }
        Ok(())
    }
    /// Parses NCBI Taxonomy, file `nodes.dmp` and returns a Taxonomy instance.
    fn parse_ncbi_nodes<P: AsRef<Path>>(nodes_file: P) -> Result<Self> {
        let reader = BufReader::new(File::open(nodes_file).unwrap());
        let mut taxonomy = Taxonomy {
            ..Default::default()
        };
        for line in reader.lines() {
            let taxon = Taxon::parse_ncbi_node(&line?)?;
            taxonomy.nodes.insert(taxon.id, taxon);
        }
        Ok(taxonomy)
    }
    /// Parses the `names.dmp` file for NCBI Taxonomy. Requires that the Taxonomy
    /// be initialised already using `parse_ncbi_nodes.
    fn parse_ncbi_names<P: AsRef<Path>>(&mut self, names_file: P) -> Result<()> {
        //todo!("Da continuare")
        if self.nodes.is_empty() {
            return Err(anyhow!("Empty Taxonomy"));
        }
        let reader = BufReader::new(File::open(names_file).unwrap());

        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split('|').map(|el| el.trim()).collect();
            let taxon_id: u32 = fields[0].parse()?;
            let mut node = match self.nodes.get_mut(&taxon_id) {
                Some(node) => node,
                None => return Err(anyhow!("Taxon Id {} Not Found", taxon_id)),
            };
            let name: String = fields[1].to_string();
            if fields[3] == SCIENTIFICNAME {
                node.name = name
            } else if fields[3] == SYNONYM {
                node.synonyms.push(name)
            }
        }
        Ok(())
    }
    /// Parses the files used in NCBI Taxonomy and found in `taxdump.tar.gz`. The files are
    /// `nodes.dmp`, `merged.dmp` and `names.dmp`. The are other files, but those are contain
    /// the informatio that is parsed and used by this data structure. A Result, with the
    /// Taxonomy instance is returned.
    pub fn from_ncbi<P: AsRef<Path>>(nodes_file: P, names_file: P, merged_file: P) -> Result<Self> {
        info!(
            "Reading NCBI taxonomy from files in {}",
            nodes_file.as_ref().parent().unwrap().display()
        );
        let mut taxonomy = Taxonomy::parse_ncbi_nodes(nodes_file)?;
        taxonomy.parse_ncbi_names(names_file)?;
        taxonomy.parse_ncbi_merged(merged_file)?;
        Ok(taxonomy)
    }

    /// Returns if a Taxon ID is inside the Taxonomy
    pub fn has_taxon_id(&self, taxon_id: &u32) -> bool {
        self.nodes.contains_key(taxon_id) || self.merged_nodes.contains_key(taxon_id)
    }

    /// Returns an Option to a reference of a Taxon, using the ID. If the taxon is one of the
    /// merged ones, the Taxon returned is the one that is now present in the Taxonomy, so the
    /// ID will not be the same as one requested.
    pub fn get_taxon(&self, taxon_id: &u32) -> Option<&Taxon> {
        let taxon = self.nodes.get(taxon_id);
        // if it's not found, check if it's part of the merged nodes
        if taxon.is_none() {
            match self.merged_nodes.get(taxon_id) {
                // if None, let the caller decide
                None => return None,
                // otherwise, pass over the correct one
                Some(new_id) => return self.nodes.get(new_id),
            };
        }
        taxon
    }

    /// Same as `get_taxon`, but a mutable reference is returned
    ///
    /// No checks are done on the correspondence between Taxon.id
    /// and the internal HashMap index, take care changing attributes.
    ///
    /// This is used when parsing NCBI taxonomy, since names are stored
    /// in another file.
    pub fn get_taxon_mut(&mut self, taxon_id: &u32) -> Option<&mut Taxon> {
        let taxon_id = match self.get_taxon(taxon_id) {
            Some(taxon) => taxon.id,
            None => return None,
        };
        self.nodes.get_mut(&taxon_id)
    }

    /// Given a list of possible ancestors and a Taxon ID, returns if the latter
    /// is a descendent of one of the ancestors IDs passed.
    pub fn is_ancestor(&self, taxon_id: &u32, anc_ids: &[u32]) -> bool {
        let taxon = match self.get_taxon(taxon_id) {
            None => {
                error!("Cannot find id {}", taxon_id);
                return false;
            }
            Some(taxon) => taxon,
        };

        let anc_ids: HashSet<&u32> = HashSet::from_iter(anc_ids);
        // checks that is not present in anc_ids
        if anc_ids.contains(&taxon.id) {
            return true;
        }

        let mut parent_id = &taxon.parent_id;

        while let Some(taxon) = self.nodes.get(parent_id) {
            if taxon.id == ROOT_TAXON || taxon.parent_id == ROOT_TAXON {
                return false;
            }
            if anc_ids.contains(&taxon.id) || anc_ids.contains(&taxon.parent_id) {
                return true;
            }
            parent_id = &taxon.parent_id;
        }

        false
    }

    /// Given a mutable reference to a HashSet of IDs, the method will extend it with
    /// all possible children.
    pub fn find_children(&self, anc_ids: &mut HashSet<u32>) {
        let new_items: HashSet<u32> = self
            .iter_taxa()
            .filter(|taxon| !anc_ids.contains(&taxon.id))
            .filter(|taxon| anc_ids.contains(&taxon.parent_id))
            .map(|taxon| taxon.id)
            .collect();
        if !new_items.is_empty() {
            anc_ids.extend(new_items);
            self.find_children(anc_ids);
        }
    }

    /// Returns an iterator to the Taxon references
    pub fn iter_taxa(&self) -> impl Iterator<Item = &Taxon> {
        self.nodes.values()
    }

    /// Returns an iterator to the merged nodes key->value
    pub fn iter_merged(&self) -> impl Iterator<Item = (&u32, &u32)> {
        self.merged_nodes.iter()
    }

    /// Given a Taxon ID, return all its ancestors as a Vector of Taxon references. If
    /// `only_ranked` is true, ranks that are either `NoRank` or `Other` are skipped
    pub fn get_taxon_lineage(&self, taxon_id: &u32, only_ranked: bool) -> Option<Vec<&Taxon>> {
        // first check that the taxon_id is present and not the root
        let taxon = match self.get_taxon(taxon_id) {
            Some(taxon) => taxon,
            None => return None,
        };
        let mut taxon_ids: Vec<&Taxon> = vec![taxon];
        let mut parent_id = &taxon.parent_id;

        while let Some(taxon) = self.get_taxon(parent_id) {
            // if only ranked check that that is not one of those to skip
            if !(only_ranked && (taxon.rank == Rank::Other || taxon.rank == Rank::NoRank)) {
                taxon_ids.push(taxon);
            }
            // change the parent Id
            parent_id = &taxon.parent_id;
            if *parent_id == ROOT_TAXON {
                break;
            }
        }
        taxon_ids.reverse();
        Some(taxon_ids)
    }

    /// Uses `get_taxon_lineage` to return a String representing the lineage of the
    /// requested Taxon ID. Only ranked Taxa are used.
    pub fn get_taxon_lineage_string(&self, taxon_id: &u32) -> Result<String> {
        let lineage = match self.get_taxon_lineage(taxon_id, true) {
            Some(lineage) => lineage,
            None => return Err(anyhow!("Taxon not found")),
        };
        let lineage: Vec<String> = lineage
            .into_iter()
            .map(|taxon| taxon.to_lineage())
            .collect();
        Ok(lineage.join("|"))
    }

    /// Finds Taxa by their name.
    pub fn find_by_name(&self, name: &str) -> Option<Vec<&Taxon>> {
        let list: Vec<&Taxon> = self
            .nodes
            .values()
            .filter(|taxon| taxon.name == name)
            .collect();
        if list.is_empty() {
            return None;
        }
        Some(list)
    }

    /// Finds all taxa with a specific name and Rank
    pub fn find_by_name_rank(&self, name: &str, rank: &Rank) -> Option<&Taxon> {
        let list: Vec<&Taxon> = self
            .nodes
            .values()
            .filter(|taxon| taxon.name == name && taxon.rank == *rank)
            .collect();
        match list.len() {
            0 => None,
            1 => Some(list[0]),
            _ => panic!("Name {} and Rank {:?} found in multiple Taxa", name, rank),
        }
    }

    /// Finds a Taxon by name, but matches the synonyms.
    pub fn find_by_synonym(&self, name: &String) -> Option<Vec<&Taxon>> {
        let list: Vec<&Taxon> = self
            .nodes
            .values()
            .filter(|&taxon| taxon.synonyms.contains(name))
            .collect();
        if list.is_empty() {
            return None;
        }
        Some(list)
    }

    /// Writes the Taxonomy into a gzipped JSON file using `serde`
    pub fn write_to_file<P: AsRef<Path>>(&self, file_name: P) -> Result<()> {
        info!("Writing taxonomy to file {}", file_name.as_ref().display());
        let json_handle = GzEncoder::new(File::create(file_name)?, Compression::default());
        to_writer(BufWriter::new(json_handle), &self)?;
        Ok(())
    }

    /// Reads the Taxonomy file from a gzipped JSON file using `serde`
    pub fn read_from_file<P: AsRef<Path>>(file_name: P) -> Result<Self> {
        info!(
            "Reading taxonomy from file {}",
            file_name.as_ref().display()
        );
        let json_handle = GzDecoder::new(File::open(file_name)?);
        let taxonomy: Taxonomy = from_reader(BufReader::new(json_handle))?;
        Ok(taxonomy)
    }

    /// Read the taxonomy exported from MGKit.
    /// The output of the command used is `taxon-utils export -t tab taxonomy.pickle taxonomy.tsv`
    pub fn read_from_mgkit<P: AsRef<Path>>(file_name: P) -> Result<Self> {
        let file_handle = super::io::open_file(file_name)?;
        let mut taxonomy = Taxonomy::default();

        for line in file_handle.lines() {
            let line = line?;
            if line.starts_with('#') {
                dbg!("{}", line);
                continue;
            } else if line.starts_with('M') {
                let (old_id, new_id) = parse_mgkit_taxon_merged(&line)?;
                taxonomy.merged_nodes.insert(old_id, new_id);
            } else {
                let mut taxon = Taxon::default();
                for (index, field) in line.splitn(5, '\t').skip(1).enumerate() {
                    match index {
                        0 => taxon.id = field.parse()?,
                        1 => taxon.parent_id = field.parse()?,
                        2 => taxon.rank = Rank::from_name(field),
                        3 => taxon.name = field.to_string(),
                        _ => break,
                    }
                }
                taxonomy.nodes.insert(taxon.id, taxon);
            }
        }
        Ok(taxonomy)
    }
}

/// Parses the lines that are for merged nodes, returning the old and the new ID.
fn parse_mgkit_taxon_merged(line: &str) -> Result<(u32, u32)> {
    let mut it = line.split('\t').skip(1).take(2);
    let (old_id, new_id) = (it.next().unwrap(), it.next().unwrap());
    Ok((old_id.parse()?, new_id.parse()?))
}

/// Parses a line from NCBI Taxonomy `merged.dmp` file, which contains IDs
/// that were merged into another Taxon. Returns a Result that unwraps into
/// a tuple, the first element is the old taxon ID and the second element the
/// ID it was merged into.s
fn parse_ncbi_merged_line(line: &str) -> Result<(u32, u32)> {
    let mut old_id = 0u32;
    let mut new_id = 0u32;
    for (index, field) in line.split('|').map(|el| el.trim()).enumerate() {
        match index {
            0 => old_id = field.parse()?,
            1 => new_id = field.parse()?,
            _ => break,
        }
    }
    Ok((old_id, new_id))
}

/// Parses PhyloPhlan Taxonomy, in fact can be used to read any Taxonomy that uses `RANK__NAME`
/// as elements and `|` (pipe as separator)
/// Requires a `Path` to the file and the Column number where the lineage is stored.
/// It skips empty lines, lines starting with a `#` or that raise an Error
///
/// TODO: Makes it more general and test with more Taxonomies
pub fn parse_phylophlan_taxonomy<P: AsRef<Path>>(
    taxonomy_file: P,
    column_n: usize,
) -> Result<Taxonomy> {
    // this buffer is used to keep the rank__name available for fast lookup
    // taxonomy is generated by then consuming this hashmap
    let mut buffer: HashMap<String, Taxon> = HashMap::new();
    let reader = super::io::open_file(taxonomy_file)?;

    // this must be global to keep track
    let mut curr_taxon_id = ROOT_TAXON;

    // starts reading the file and skips the comments or line that generate errors
    for line in reader.lines() {
        let line = match line {
            Ok(line) => line,
            Err(_) => String::from("#"),
        };
        if line.starts_with('#') {
            continue;
        }
        // in each line a field is there for the names, while the next has the Ids for the
        // NCBI taxonomy for entries that have it. Since it's not important we won't use it
        // then brings back the iterator to one element and removes whitespace
        let field: &str = line.split('\t').nth(column_n).unwrap().trim();

        // the first element is always kingdom
        // if the Rank is Kingdom, the parent_id must be 1 (root)
        let mut parent_id = ROOT_TAXON;

        for taxon_field in field.split('|') {
            let taxon_field: String = taxon_field.into();
            let taxon = buffer.get(&taxon_field);
            match taxon {
                // Already present, skips but changes the parent_id so next iteration will have
                // the correct one if needed
                Some(taxon) => {
                    parent_id = taxon.id;
                    continue;
                }
                // Will add a new Taxon on this iteration, so increments for the taxon_id value
                None => curr_taxon_id += 1,
            }
            if let Some((rank_char, name)) = taxon_field.split_once("__") {
                // splits the name into its components and gets a new taxon
                let rank = Rank::from_prefix(rank_char);
                let taxon = Taxon {
                    id: curr_taxon_id,
                    parent_id,
                    rank,
                    name: name.replace('_', "").to_string(),
                    ..Default::default()
                };
                // parent Id is changed, in case this is not the last in the lineage
                parent_id = curr_taxon_id;
                // inserts into the buffer
                buffer.insert(taxon_field, taxon);
            } else {
                // in case some error arise, skips this line
                continue;
            }
        }
    }
    // Creates the taxonomy and populate it
    let mut taxonomy = Taxonomy::default();
    // by consuming the buffer values
    for value in buffer.into_values() {
        taxonomy.nodes.insert(value.id, value);
    }
    // and adds the root Taxon
    taxonomy.nodes.insert(
        1,
        Taxon {
            id: 1,
            parent_id: 1,
            rank: Rank::NoRank,
            name: "root".to_string(),
            ..Default::default()
        },
    );
    // finally returns the taxonomy and moves ownership to the calling
    Ok(taxonomy)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_empty_taxonomy() {
        let mut taxonomy = Taxonomy {
            ..Default::default()
        };
        let result = taxonomy.parse_ncbi_names("test");
        assert!(result.is_err());
    }

    #[test]
    fn test_rank() {
        assert_eq!(Rank::default(), Rank::NoRank);
        assert_eq!(Rank::from_name("kingdom"), Rank::Kingdom);
        assert_eq!(Rank::Kingdom.to_name(), "kingdom");
        assert_eq!(Rank::Kingdom.to_prefix(), "k");
        assert_eq!(Rank::Other.to_prefix(), "x");
        assert_eq!(Rank::NoRank.to_prefix(), "x");
        assert_eq!(Rank::from_prefix("g"), Rank::Genus);
        assert_eq!(Rank::from_prefix("x"), Rank::NoRank);
    }

    #[test]
    fn test_taxon_parse_ncbi_node() {
        let line = "2173\t|\t2172\t|\tno rank\t|\tmore";
        let taxon = Taxon::parse_ncbi_node(line).unwrap();
        assert_eq!(taxon.id, 2173);
        assert_eq!(taxon.parent_id, 2172);
        assert_eq!(taxon.rank, Rank::NoRank);
    }

    #[test]
    fn test_taxon_parse_ncbi_node_fail() {
        let line = "2173\t|\tno rank\t|\tmore";
        let taxon = Taxon::parse_ncbi_node(line);
        assert!(taxon.is_err());
    }

    #[test]
    fn test_taxon_to_lineage() {
        let taxon = Taxon {
            name: "Genus species".to_string(),
            rank: Rank::Species,
            ..Default::default()
        };
        assert_eq!(taxon.to_lineage(), "s__Genus_species");
    }

    #[test]
    fn test_taxonomy_has_taxon_id() {
        let taxon = Taxon {
            id: 2,
            ..Default::default()
        };
        let mut taxonomy = Taxonomy::default();
        taxonomy.inset_taxon(taxon);
        assert!(taxonomy.has_taxon_id(&2));
    }

    #[test]
    fn test_taxonomy_get_taxon() {
        let taxon = Taxon {
            id: 2,
            ..Default::default()
        };
        let mut taxonomy = Taxonomy::default();
        taxonomy.inset_taxon(taxon);
        assert!(taxonomy.get_taxon(&2).is_some());
        assert!(taxonomy.get_taxon(&3).is_none());
    }

    #[test]
    fn test_taxonomy_get_taxon_mut() {
        let taxon = Taxon {
            id: 2173,
            parent_id: 2172,
            ..Default::default()
        };
        let mut taxonomy = Taxonomy::default();
        taxonomy.inset_taxon(taxon);
        let mut t = taxonomy.get_taxon_mut(&2173).unwrap();
        assert_eq!(t.parent_id, 2172);
        t.parent_id = 2;
        let t = taxonomy.get_taxon(&2173).unwrap();
        assert_eq!(t.parent_id, 2);
    }
    #[test]
    fn test_rank_ord() {
        assert!(Rank::SuperKingdom < Rank::Kingdom);
        assert!(Rank::Kingdom < Rank::Phylum);
        assert!(Rank::Phylum < Rank::Class);
        assert!(Rank::Class < Rank::Order);
        assert!(Rank::Order < Rank::Family);
        assert!(Rank::Family < Rank::Genus);
        assert!(Rank::Genus < Rank::Species);
        assert!(Rank::Species < Rank::Strain);
        assert!(Rank::Strain < Rank::NoRank);
        assert!(Rank::NoRank < Rank::Other);
    }
}
