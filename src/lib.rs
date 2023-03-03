use anyhow::{Context, Result};
use regex::bytes::{CaptureLocations, Regex};
use seq_geom_parser::{FragmentGeomDesc, GeomLen, GeomPiece, NucStr};

/*
#[derive(Debug, Display, Clone, Copy)]
enum TypedGenomLen{

}
*/

#[derive(Debug)]
pub struct FragmentRegexDesc {
    pub r1_cginfo: Vec<GeomPiece>,
    pub r2_cginfo: Vec<GeomPiece>,
    pub r1_re: Regex,
    pub r2_re: Regex,
    r1_clocs: CaptureLocations,
    r2_clocs: CaptureLocations,
}

pub struct SeqPair {
    pub s1: String,
    pub s2: String,
}

impl SeqPair {
    pub fn new() -> Self {
        SeqPair { s1: String::new(), s2: String::new() }
    }

    fn clear(&mut self) {
        self.s1.clear();
        self.s2.clear();
    }
}

impl FragmentRegexDesc {
    pub fn parse_into(&mut self, r1: &[u8], r2: &[u8], sp: &mut SeqPair) {
        sp.clear();
        let m1 = self.r1_re.captures_read(&mut self.r1_clocs, r1);
        let m2 = self.r2_re.captures_read(&mut self.r2_clocs, r2);

        let s1 = unsafe { std::str::from_utf8_unchecked(r1) };
        let s2 = unsafe { std::str::from_utf8_unchecked(r2) };

        // walk over the capture groups 
        
        // if we expected to capture the whole thing
        if self.r1_cginfo.len() > 1 {
           sp.s1 = String::from(s1);
        } else {
            for cl in 1..self.r1_clocs.len() {
                if let Some(g) = self.r1_clocs.get(cl) {
                    sp.s1 += s1.get(g.0..g.1).unwrap();
                }
            }
        }

        // if we expected to capture the whole thing
        if self.r2_cginfo.len() > 1 {
           sp.s2 = String::from(s2);
        } else {
            for cl in 1..self.r2_clocs.len() {
                if let Some(g) = self.r2_clocs.get(cl) {
                    sp.s2 += s2.get(g.0..g.1).unwrap();
                }
            }
        }
    }
}

/// Extension methods for FragmentGeomDesc
pub trait FragmentGeomDescExt {
    fn as_regex(&self) -> Result<FragmentRegexDesc, anyhow::Error>;
}

fn geom_piece_as_regex_string(gp: &GeomPiece) -> (String, Option<GeomPiece>) {
    let mut rep = String::from("");
    let mut geo = None;
    match gp {
        // single lengths
        GeomPiece::Discard(GeomLen::Bounded(x)) => {
            rep += &format!(r#"[ACGTN]{{{}}}"#, x);
            // don't need to capture
        }
        GeomPiece::Barcode(GeomLen::Bounded(x)) => {
            rep += &format!(r#"([ACGTN]{{{}}})"#, x);
            geo = Some(gp.clone());
        }
        GeomPiece::Umi(GeomLen::Bounded(x)) => {
            rep += &format!(r#"([ACGTN]{{{}}})"#, x);
            geo = Some(gp.clone());
        }
        GeomPiece::ReadSeq(GeomLen::Bounded(x)) => {
            rep += &format!(r#"([ACGTN]{{{}}})"#, x);
            geo = Some(gp.clone());
        }
        // length ranges
        GeomPiece::Discard(GeomLen::BoundedRange(l, h)) => {
            rep += &format!(r#"[ACGTN]{{{},{}}}"#, l, h);
            // don't need to capture
        }
        GeomPiece::Barcode(GeomLen::BoundedRange(l, h)) => {
            rep += &format!(r#"([ACGTN]{{{},{}}})"#, l, h);
            geo = Some(gp.clone());
        }
        GeomPiece::Umi(GeomLen::BoundedRange(l, h)) => {
            rep += &format!(r#"([ACGTN]{{{},{}}})"#, l, h);
            geo = Some(gp.clone());
        }
        GeomPiece::ReadSeq(GeomLen::BoundedRange(l, h)) => {
            rep += &format!(r#"([ACGTN]{{{},{}}})"#, l, h);
            geo = Some(gp.clone());
        }
        // fixed sequence
        GeomPiece::Fixed(NucStr::Seq(s)) => {
            // no caputre group because no need to capture this
            // right now
            rep += &format!(r#"{}"#, s);
        }
        // unbounded pieces
        GeomPiece::Discard(GeomLen::Unbounded) => {
            rep += &format!(r#"[ACGTN]*"#);
        }
        GeomPiece::Barcode(GeomLen::Unbounded) => {
            rep += &format!(r#"([ACGTN]*)"#);
            geo = Some(gp.clone());
        }
        GeomPiece::Umi(GeomLen::Unbounded) => {
            rep += &format!(r#"([ACGTN]*)"#);
            geo = Some(gp.clone());
        }
        GeomPiece::ReadSeq(GeomLen::Unbounded) => {
            rep += &format!(r#"([ACGTN]*)"#);
            geo = Some(gp.clone());
        }
    }
    (rep, geo)
}

impl FragmentGeomDescExt for FragmentGeomDesc {
    fn as_regex(&self) -> Result<FragmentRegexDesc, anyhow::Error> {
        let mut r1_re_str = String::from("");
        let mut r1_cginfo = Vec::<GeomPiece>::new();
        for geo_piece in &self.read1_desc {
            let (str_piece, geo_len) = geom_piece_as_regex_string(geo_piece);
            r1_re_str += &str_piece;
            if let Some(elem) = geo_len {
                r1_cginfo.push(elem);
            }
        }

        let mut r2_re_str = String::from("");
        let mut r2_cginfo = Vec::<GeomPiece>::new();
        for geo_piece in &self.read2_desc {
            let (str_piece, geo_len) = geom_piece_as_regex_string(geo_piece);
            r2_re_str += &str_piece;
            if let Some(elem) = geo_len {
                r2_cginfo.push(elem);
            }
        }

        let r1_re = Regex::new(&r1_re_str)
            .with_context(|| format!("Could not compile {} into regex description", r1_re_str))?;
        let r2_re = Regex::new(&r2_re_str)
            .with_context(|| format!("Could not compile {} into regex description", r2_re_str))?;
        println!("{:?}", r1_cginfo);
        println!("{:?}", r2_cginfo);

        let cloc1 = r1_re.capture_locations();
        let cloc2 = r2_re.capture_locations();

        Ok(FragmentRegexDesc {
            r1_cginfo,
            r2_cginfo,
            r1_re,
            r2_re,
            r1_clocs: cloc1,
            r2_clocs: cloc2
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
