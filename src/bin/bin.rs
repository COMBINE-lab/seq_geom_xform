use seq_geom_parser::FragmentGeomDesc;
use seq_geom_xform::FragmentGeomDescExt;

use needletail::{parse_fastx_file, FastxReader, Sequence};

use regex::bytes::Regex;
use std::env;

fn main() {
    let gd = std::env::args().nth(1).unwrap();
    let filename1 = std::env::args().nth(2).unwrap();
    let filename2 = std::env::args().nth(3).unwrap();

    let geo = FragmentGeomDesc::try_from(gd.as_str()).unwrap();

    let mut parsed_records = seq_geom_xform::SeqPair::new();

    if let Ok(mut geo_re) = geo.as_regex() {
        println!(
            "geometry as regex = Read1 : {:?}, Read2 : {:?}",
            geo_re.r1_re, geo_re.r2_re
        );

        let mut reader = parse_fastx_file(&filename1).expect("valid path/file");
        let mut reader2 = parse_fastx_file(&filename2).expect("valid path/file");

        while let (Some(record), Some(record2)) = (reader.next(), reader2.next()) {
            let seqrec = record.expect("invalid record");
            let seqrec2 = record2.expect("invalid record");
            
            geo_re.parse_into(seqrec.sequence(), seqrec2.sequence(), &mut parsed_records);
            println!("parsed (r1, r2) = ({}, {})", parsed_records.s1, parsed_records.s2);
        }
    }

    //let r = regex::Regex::new(r"^([ACGTNacgtn]{16})([ACGTNacgtn]{12})[ACGTNacgtn]*").unwrap();
    //let s = "TGGTGGCCAGCGCCCCCTGCTGGCGCCGGGGCACTGCAGGGCCCTCTTGCTTACTGTATAGTGGTGGCACGCCGCCTGCTGGCAGCTAGGG";
    //if let Some(cap) = r.captures(s) {
    //println!("captures = {:?}", cap);
    //}
}
