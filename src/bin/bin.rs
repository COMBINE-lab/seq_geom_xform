use seq_geom_xform::FragmentGeomDescExt;
use seq_geom_parser::FragmentGeomDesc;

use needletail::{parse_fastx_file, Sequence, FastxReader};

use regex::bytes::Regex;
use std::env;

fn main() { 
    
    let gd = std::env::args().nth(1).unwrap();
    let filename1 = std::env::args().nth(2).unwrap();
    let filename2 = std::env::args().nth(3).unwrap();

    let geo = FragmentGeomDesc::try_from(gd.as_str()).unwrap();

    if let Ok(geo_re) = geo.as_regex() {
        println!("geometry as regex = Read1 : {:?}, Read2 : {:?}", geo_re.r1_re, geo_re.r2_re);

        let mut clocs = geo_re.r1_re.capture_locations();
        let mut clocs2 = geo_re.r2_re.capture_locations();

        let mut reader = parse_fastx_file(&filename1).expect("valid path/file");
        let mut reader2 = parse_fastx_file(&filename2).expect("valid path/file");

        while let (Some(record), Some(record2)) = (reader.next(), reader2.next()) {
            let seqrec = record.expect("invalid record");
            let seqrec2 = record2.expect("invalid record");
            let m1 = geo_re.r1_re.captures_read(&mut clocs, seqrec.sequence());
            let m2 = geo_re.r2_re.captures_read(&mut clocs2, seqrec2.sequence());
            
            let s1 = std::str::from_utf8(seqrec.sequence()).unwrap();
            let s2 = std::str::from_utf8(seqrec2.sequence()).unwrap();

            if clocs.len() > 1 {
                for cl in 1..clocs.len() {
                    if let Some(g) = clocs.get(cl) {
                        println!("read 1, capture group {} : {:?}, {}", cl, g, s1.get(g.0..g.1).unwrap());
                    }
                }
            } else if clocs.len() == 1 {
                println!("m1 = {:#?}", &clocs);
            }

            if clocs2.len() > 1 {
                for cl in 1..clocs2.len() {
                    if let Some(g) = clocs2.get(cl) {
                        println!("read 1, capture group {} : {:?}, {}", cl, g, s2.get(g.0..g.1).unwrap());
                    }
                }
 
            } else if clocs2.len() == 1 {
                println!("m2 = {:#?}", &clocs2);
            }

        }
    }

    //let r = regex::Regex::new(r"^([ACGTNacgtn]{16})([ACGTNacgtn]{12})[ACGTNacgtn]*").unwrap();
    //let s = "TGGTGGCCAGCGCCCCCTGCTGGCGCCGGGGCACTGCAGGGCCCTCTTGCTTACTGTATAGTGGTGGCACGCCGCCTGCTGGCAGCTAGGG"; 
    //if let Some(cap) = r.captures(s) {
        //println!("captures = {:?}", cap);
    //}
}
