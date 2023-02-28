use seq_geom_parser::{FragmentGeomDesc, GeomPiece, GeomLen};
use regex;
use anyhow::{Context, Result};

#[derive(Debug)]
pub struct FragmentRegexDesc {
    pub r1_re : regex::Regex,
    pub r2_re : regex::Regex
}


/// Extension methods for FragmentGeomDesc 
pub trait FragmentGeomDescExt {
    fn as_regex(&self) -> Result<FragmentRegexDesc, anyhow::Error>;
}

fn geom_piece_as_regex_string(gp: &GeomPiece) -> String {
    let mut rep = String::from("");
    match gp {
        GeomPiece::Discard(GeomLen::Bounded(x)) => {
            rep += &format!("[ACGTNacgtn]{{{}}}", x);
        }
        GeomPiece::Barcode(GeomLen::Bounded(x)) => {
            rep += &format!(r#"([ACGTNacgtn]{{{}}})"#, x);
        }
        GeomPiece::Umi(GeomLen::Bounded(x)) => {
            rep += &format!(r#"([ACGTNacgtn]{{{}}})"#, x);
        }
        GeomPiece::ReadSeq(GeomLen::Bounded(x)) => {
            rep += &format!(r#"([ACGTNacgtn]{{{}}})"#, x);
        }
        GeomPiece::Discard(GeomLen::Unbounded) => {
            rep += &format!(r#"[ACGTNacgtn]*"#);
        }
        GeomPiece::Barcode(GeomLen::Unbounded) => {
            rep += &format!(r#"([ACGTNacgtn]*)"#);
        }
        GeomPiece::Umi(GeomLen::Unbounded) => {
            rep += &format!(r#"([ACGTNacgtn]*)"#);
        }
        GeomPiece::ReadSeq(GeomLen::Unbounded) => {
            rep += &format!(r#"([ACGTNacgtn]*)"#);
        }
    }
    rep
}

impl FragmentGeomDescExt for FragmentGeomDesc {
    fn as_regex(&self) -> Result<FragmentRegexDesc, anyhow::Error> {
        let mut r1_re_str = String::from("");
        for geo_piece in &self.read1_desc {
           r1_re_str += &geom_piece_as_regex_string(geo_piece); 
        }
        let mut r2_re_str = String::from("");
        for geo_piece in &self.read2_desc {
           r2_re_str += &geom_piece_as_regex_string(geo_piece); 
        }

        let r1_re = regex::Regex::new(&r1_re_str).with_context(|| 
            format!("Could not compile {} into regex description", r1_re_str))?;
        let r2_re = regex::Regex::new(&r2_re_str).with_context(|| 
            format!("Could not compile {} into regex description", r2_re_str))?;
        Ok(FragmentRegexDesc{r1_re, r2_re})
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
