mod common;

#[test]
fn test_ranked_taxon() {
    let taxonomy = common::setup_taxonomy().unwrap();

    let ranked = taxonomy.get_ranked_taxon(2173, bio_rascal::taxon::Rank::Genus);
    // 2172
    assert_eq!(ranked.unwrap().id, 2172);
    let ranked = taxonomy.get_ranked_taxon(2173, bio_rascal::taxon::Rank::SuperKingdom);
    // 2157
    assert_eq!(ranked.unwrap().id, 2157);

    let ranked = taxonomy.get_ranked_taxon(9606, bio_rascal::taxon::Rank::Phylum);
    // 7711
    assert_eq!(ranked.unwrap().id, 7711);

    let ranked = taxonomy.get_ranked_taxon(960621412, bio_rascal::taxon::Rank::Phylum);
    // None
    assert!(ranked.is_none());
    
    let ranked = taxonomy.get_ranked_taxon(2172, bio_rascal::taxon::Rank::Species);
    // 2172
    assert_eq!(ranked.unwrap().id, 2172);

}