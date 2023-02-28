use seq_geom_xform::FragmentGeomDescExt;
use seq_geom_parser::FragmentGeomDesc;
use regex;
use std::env;

fn main() { 
    
    let gd = std::env::args().nth(1).unwrap();
    let geo = FragmentGeomDesc::try_from(gd.as_str()).unwrap();
    if let Ok(geo_re) = geo.as_regex() {
        println!("geometry as regex = Read1 : {:?}, Read2 : {:?}", geo_re.r1_re, geo_re.r2_re);
    }

    //let r = regex::Regex::new(r"^([ACGTNacgtn]{16})([ACGTNacgtn]{12})[ACGTNacgtn]*").unwrap();
    //let s = "TGGTGGCCAGCGCCCCCTGCTGGCGCCGGGGCACTGCAGGGCCCTCTTGCTTACTGTATAGTGGTGGCACGCCGCCTGCTGGCAGCTAGGG"; 
    //if let Some(cap) = r.captures(s) {
        //println!("captures = {:?}", cap);
    //}
}
