use anyhow::{bail, Context, Result};
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
        SeqPair {
            s1: String::new(),
            s2: String::new(),
        }
    }

    fn clear(&mut self) {
        self.s1.clear();
        self.s2.clear();
    }
}

impl Default for SeqPair {
    fn default() -> Self {
        Self::new()
    }
}

const BOUNDED_RANGE_LIMIT: u32 = 4;
const VAR_LEN_BC_PADDING: &[&str] = &["A", "AC", "AAG", "AAAT"];

/// Builds the parsed output string `s` given the `CaptureLocations` `clocs`,
/// the expected captured `GeomPiece`s `gpieces` and the input string `r`.  This function
/// returns true if the parse was succesful (the captured groups are what is expected)
/// and false otherwise.
#[inline(always)]
fn parse_single_read(
    clocs: &CaptureLocations,
    gpieces: &Vec<GeomPiece>,
    r: &str,
    outstr: &mut String,
) -> bool {
    // if we expected to capture the whole thing
    if gpieces.len() == 1 {
        // then just copy it over into our string
        outstr.push_str(r);
    } else {
        // otherwise, process each capture group
        for cl in 1..clocs.len() {
            if let Some(g) = clocs.get(cl) {
                outstr.push_str(r.get(g.0..g.1).unwrap());

                match gpieces.get(cl - 1) {
                    // if we captured some variable length piece of geometry
                    // then we have to apply the appropriate padding so that
                    // we can pass the result to a non-variable length parser.
                    Some(GeomPiece::Barcode(GeomLen::BoundedRange(_l, h)))
                    | Some(GeomPiece::Umi(GeomLen::BoundedRange(_l, h)))
                    | Some(GeomPiece::ReadSeq(GeomLen::BoundedRange(_l, h))) => {
                        let captured_len = (g.1 - g.0) as usize;
                        outstr.push_str(VAR_LEN_BC_PADDING[(*h as usize) - (captured_len)]);
                    }
                    _ => {
                        // fixed length, do nothing
                    }
                }
            } else {
                return false;
            }
        }
    }
    true
}


fn get_simplified_piscem_string(geo_pieces: &[GeomPiece]) -> String {
    let mut rep = String::new();
    for gp in geo_pieces {
        match gp {
            GeomPiece::Discard(GeomLen::Bounded(x)) => {
                rep += &format!("x[{}]", x);
            }
            GeomPiece::Barcode(GeomLen::Bounded(x)) => {
                rep += &format!("b[{}]", x);
            }
            GeomPiece::Umi(GeomLen::Bounded(x)) => {
                rep += &format!("u[{}]", x);
            }
            GeomPiece::ReadSeq(GeomLen::Bounded(x)) => {
                rep += &format!("r[{}]", x);
            }
            GeomPiece::Discard(GeomLen::BoundedRange(_l, h)) => {
                rep += &format!("x[{}]", h+1);
            }
            GeomPiece::Barcode(GeomLen::BoundedRange(_l, h)) => {
                rep += &format!("b[{}]", h+1);
            }
            GeomPiece::Umi(GeomLen::BoundedRange(_l, h)) => {
                rep += &format!("u[{}]", h+1);
            }
            GeomPiece::ReadSeq(GeomLen::BoundedRange(_l, h)) => {
                rep += &format!("r[{}]", h+1);
            }
            GeomPiece::Discard(GeomLen::Unbounded) => {
                rep += "x:";
            }
            GeomPiece::Barcode(GeomLen::Unbounded) => {
                rep += "b:";
            }
            GeomPiece::Umi(GeomLen::Unbounded) => {
                rep += "u:";
            }
            GeomPiece::ReadSeq(GeomLen::Unbounded) => {
                rep += "r:";
            }
            _ => { unimplemented!(); }
        }
    }
    rep
}

impl FragmentRegexDesc {
    /// Parses the read pair `r1` and `r2` in accordance with the geometry specified
    /// in `self`.  The resulting parse, if successful, is placed into the output
    /// `sp`. This function returns true if the entire *pair* of reads was parsed succesfully,
    /// and false otherwise. If the parse is not successful, nothing can be assumed about
    /// the contents of `sp`.
    pub fn parse_into(&mut self, r1: &[u8], r2: &[u8], sp: &mut SeqPair) -> bool {
        sp.clear();
        let _m1 = self.r1_re.captures_read(&mut self.r1_clocs, r1);
        let _m2 = self.r2_re.captures_read(&mut self.r2_clocs, r2);

        let s1 = unsafe { std::str::from_utf8_unchecked(r1) };
        let s2 = unsafe { std::str::from_utf8_unchecked(r2) };

        let parsed_r1 = parse_single_read(&self.r1_clocs, &self.r1_cginfo, s1, &mut sp.s1);
        if parsed_r1 {
            parse_single_read(&self.r2_clocs, &self.r2_cginfo, s2, &mut sp.s2)
        } else {
            false
        }
    }

    pub fn get_simplified_piscem_description_string(&self) -> String {
        let mut rep= String::from("");
        if !self.r1_cginfo.is_empty() {
            let d = get_simplified_piscem_string(&self.r1_cginfo);
            rep += &format!("1{{{}}}", d);
        }
        if !self.r2_cginfo.is_empty() {
            let d = get_simplified_piscem_string(&self.r2_cginfo);
            rep += &format!("2{{{}}}", d);
        }
        rep 
    }
}

/// Extension methods for FragmentGeomDesc
pub trait FragmentGeomDescExt {
    fn as_regex(&self) -> Result<FragmentRegexDesc, anyhow::Error>;
}

fn geom_piece_as_regex_string(gp: &GeomPiece) -> Result<(String, Option<GeomPiece>)> {
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
            if h - l > BOUNDED_RANGE_LIMIT {
                bail!("Bounded range can have variable width at most {} but the current element {:?} has variable width {}.",
                    BOUNDED_RANGE_LIMIT, &gp, h-l);
            }
            rep += &format!(r#"[ACGTN]{{{},{}}}"#, l, h);
            // don't need to capture
        }
        GeomPiece::Barcode(GeomLen::BoundedRange(l, h)) => {
            if h - l > BOUNDED_RANGE_LIMIT {
                bail!("Bounded range can have variable width at most {} but the current element {:?} has variable width {}.",
                    BOUNDED_RANGE_LIMIT, &gp, h-l);
            }
            rep += &format!(r#"([ACGTN]{{{},{}}})"#, l, h);
            geo = Some(gp.clone());
        }
        GeomPiece::Umi(GeomLen::BoundedRange(l, h)) => {
            if h - l > BOUNDED_RANGE_LIMIT {
                bail!("Bounded range can have variable width at most {} but the current element {:?} has variable width {}.",
                    BOUNDED_RANGE_LIMIT, &gp, h-l);
            }
            rep += &format!(r#"([ACGTN]{{{},{}}})"#, l, h);
            geo = Some(gp.clone());
        }
        GeomPiece::ReadSeq(GeomLen::BoundedRange(l, h)) => {
            if h - l > BOUNDED_RANGE_LIMIT {
                bail!("Bounded range can have variable width at most {} but the current element {:?} has variable width {}.",
                    BOUNDED_RANGE_LIMIT, &gp, h-l);
            }
            rep += &format!(r#"([ACGTN]{{{},{}}})"#, l, h);
            geo = Some(gp.clone());
        }
        // fixed sequence
        GeomPiece::Fixed(NucStr::Seq(s)) => {
            // no caputre group because no need to capture this
            // right now
            rep += s;
        }
        // unbounded pieces
        GeomPiece::Discard(GeomLen::Unbounded) => {
            rep += r#"[ACGTN]*"#;
        }
        GeomPiece::Barcode(GeomLen::Unbounded) => {
            rep += r#"([ACGTN]*)"#;
            geo = Some(gp.clone());
        }
        GeomPiece::Umi(GeomLen::Unbounded) => {
            rep += r#"([ACGTN]*)"#;
            geo = Some(gp.clone());
        }
        GeomPiece::ReadSeq(GeomLen::Unbounded) => {
            rep += r#"([ACGTN]*)"#;
            geo = Some(gp.clone());
        }
    }
    Ok((rep, geo))
}

impl FragmentGeomDescExt for FragmentGeomDesc {
    fn as_regex(&self) -> Result<FragmentRegexDesc, anyhow::Error> {
        let mut r1_re_str = String::from("");
        let mut r1_cginfo = Vec::<GeomPiece>::new();
        for geo_piece in &self.read1_desc {
            let (str_piece, geo_len) = geom_piece_as_regex_string(geo_piece)?;
            r1_re_str += &str_piece;
            if let Some(elem) = geo_len {
                r1_cginfo.push(elem);
            }
        }

        let mut r2_re_str = String::from("");
        let mut r2_cginfo = Vec::<GeomPiece>::new();
        for geo_piece in &self.read2_desc {
            let (str_piece, geo_len) = geom_piece_as_regex_string(geo_piece)?;
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
            r2_clocs: cloc2,
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
