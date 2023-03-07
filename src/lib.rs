use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::thread;

use anyhow::{bail, Context, Result};
use regex::bytes::{CaptureLocations, Regex};
use seq_geom_parser::{FragmentGeomDesc, GeomLen, GeomPiece, NucStr};

use needletail::{parse_fastx_file, Sequence};
use tracing::info;

use nix::sys::stat;
use nix::unistd;
use tempfile::tempdir;

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
                        let captured_len = g.1 - g.0;
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
            // NOTE: the + 1 in the rules below assumes we will
            // only ever have variable width of geometry pieces
            // of at most 4 bases. If we need to every move
            // beyond that, this code will have to be generalized.
            GeomPiece::Discard(GeomLen::BoundedRange(_l, h)) => {
                rep += &format!("x[{}]", h + 1);
            }
            GeomPiece::Barcode(GeomLen::BoundedRange(_l, h)) => {
                rep += &format!("b[{}]", h + 1);
            }
            GeomPiece::Umi(GeomLen::BoundedRange(_l, h)) => {
                rep += &format!("u[{}]", h + 1);
            }
            GeomPiece::ReadSeq(GeomLen::BoundedRange(_l, h)) => {
                rep += &format!("r[{}]", h + 1);
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
            _ => {
                unimplemented!();
            }
        }
    }
    rep
}

fn get_simplified_geo(gp: &GeomPiece) -> GeomPiece {
    match gp {
        // NOTE: the + 1 in the rules below assumes we will
        // only ever have variable width of geometry pieces
        // of at most 4 bases. If we need to every move
        // beyond that, this code will have to be generalized.
        GeomPiece::Discard(GeomLen::BoundedRange(_l, h)) => {
            GeomPiece::Discard(GeomLen::Bounded(h + 1))
        }
        GeomPiece::Barcode(GeomLen::BoundedRange(_l, h)) => {
            GeomPiece::Barcode(GeomLen::Bounded(h + 1))
        }
        GeomPiece::Umi(GeomLen::BoundedRange(_l, h)) => GeomPiece::Umi(GeomLen::Bounded(h + 1)),
        GeomPiece::ReadSeq(GeomLen::BoundedRange(_l, h)) => {
            GeomPiece::ReadSeq(GeomLen::Bounded(h + 1))
        }
        _ => gp.clone(),
    }
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

    pub fn get_simplified_geo_desc(&self) -> FragmentGeomDesc {
        FragmentGeomDesc {
            read1_desc: self
                .r1_cginfo
                .iter()
                .map(get_simplified_geo)
                .collect::<Vec<GeomPiece>>(),
            read2_desc: self
                .r2_cginfo
                .iter()
                .map(get_simplified_geo)
                .collect::<Vec<GeomPiece>>(),
        }
    }

    pub fn get_simplified_description_string(&self) -> String {
        let mut rep = String::from("");
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
        let mut r1_re_str = String::from("^");
        let mut r1_cginfo = Vec::<GeomPiece>::new();
        for geo_piece in &self.read1_desc {
            let (str_piece, geo_len) = geom_piece_as_regex_string(geo_piece)?;
            r1_re_str += &str_piece;
            if let Some(elem) = geo_len {
                r1_cginfo.push(elem);
            }
        }

        // This seems to lead to a slight performance improvement, but consider if
        // we really want to do this.  This checks if the last GeomPiece in the
        // current description is bounded or not.  If so (i.e. if it is bounded), then
        // we add an unbounded `Discard` GeomPiece to the end followed by the
        // end of string anchor.  This anchoring of the regex (seemingly) makes matching a
        // little bit faster.
        if let Some(geo_piece) = &self.read1_desc.last() {
            if geo_piece.is_bounded() {
                let (str_piece, _geo_len) =
                    geom_piece_as_regex_string(&GeomPiece::Discard(GeomLen::Unbounded))?;
                r1_re_str += &str_piece;
            }
        }
        r1_re_str.push('$');

        let mut r2_re_str = String::from("^");
        let mut r2_cginfo = Vec::<GeomPiece>::new();
        for geo_piece in &self.read2_desc {
            let (str_piece, geo_len) = geom_piece_as_regex_string(geo_piece)?;
            r2_re_str += &str_piece;
            if let Some(elem) = geo_len {
                r2_cginfo.push(elem);
            }
        }

        // This seems to lead to a slight performance improvement, but consider if
        // we really want to do this.  This checks if the last GeomPiece in the
        // current description is bounded or not.  If so (i.e. if it is bounded), then
        // we add an unbounded `Discard` GeomPiece to the end followed by the
        // end of string anchor.  This anchoring of the regex (seemingly) makes matching a
        // little bit faster.
        if let Some(geo_piece) = &self.read2_desc.last() {
            if geo_piece.is_bounded() {
                let (str_piece, _geo_len) =
                    geom_piece_as_regex_string(&GeomPiece::Discard(GeomLen::Unbounded))?;
                r2_re_str += &str_piece;
            }
        }
        r2_re_str.push('$');

        let r1_re = Regex::new(&r1_re_str)
            .with_context(|| format!("Could not compile {} into regex description", r1_re_str))?;
        let r2_re = Regex::new(&r2_re_str)
            .with_context(|| format!("Could not compile {} into regex description", r2_re_str))?;

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

/// The information we get back from an xform function
/// that tells us the relevant information about the
/// transformation taking place
#[derive(Debug)]
pub struct FifoXFormData {
    pub r1_fifo: PathBuf,
    pub r2_fifo: PathBuf,
    pub join_handle: thread::JoinHandle<Result<()>>,
}

pub fn xform_read_pairs_to_file(
    mut geo_re: FragmentRegexDesc,
    r1: Vec<PathBuf>,
    r2: Vec<PathBuf>,
    r1_ofile: PathBuf,
    r2_ofile: PathBuf,
) -> Result<()> {
    let f1 = File::create(r1_ofile).expect("Unable to open read 1 file");
    let f2 = File::create(r2_ofile).expect("Unable to open read 2 file");

    let mut stream1 = BufWriter::new(f1);
    let mut stream2 = BufWriter::new(f2);

    let mut parsed_records = SeqPair::new();
    for (filename1, filename2) in r1.iter().zip(r2.iter()) {
        let mut reader = parse_fastx_file(filename1).expect("valid path/file");
        let mut reader2 = parse_fastx_file(filename2).expect("valid path/file");

        while let (Some(record), Some(record2)) = (reader.next(), reader2.next()) {
            let seqrec = record.expect("invalid record");
            let seqrec2 = record2.expect("invalid record");

            if geo_re.parse_into(seqrec.sequence(), seqrec2.sequence(), &mut parsed_records) {
                unsafe {
                    std::write!(
                        &mut stream1,
                        ">{}\n{}\n",
                        std::str::from_utf8_unchecked(seqrec.id()),
                        parsed_records.s1
                    )
                    .expect("couldn't write output to file 1");
                    std::write!(
                        &mut stream2,
                        ">{}\n{}\n",
                        std::str::from_utf8_unchecked(seqrec2.id()),
                        parsed_records.s2
                    )
                    .expect("couldn't write output to file 2");
                }
            }
        }
    }
    Ok(())
}

pub fn xform_read_pairs_to_fifo(
    geo_re: FragmentRegexDesc,
    r1: Vec<PathBuf>,
    r2: Vec<PathBuf>,
) -> Result<FifoXFormData> {
    if r1.len() != r2.len() {
        bail!(
            "The number of R1 files ({}) must match the number of R2 files ({})",
            r1.len(),
            r2.len()
        );
    }

    let tmp_dir = tempdir().unwrap();
    let r1_fifo = tmp_dir.path().join("r1.pipe");
    let r2_fifo = tmp_dir.path().join("r2.pipe");

    // create new fifo and give read, write and execute rights to the owner
    match unistd::mkfifo(&r1_fifo, stat::Mode::S_IRWXU) {
        Ok(_) => {
            info!("created {:?}", r1_fifo);
            assert!(std::path::Path::new(&r1_fifo).exists());
        }
        Err(err) => bail!("Error creating read 1 fifo: {}", err),
    }
    // create new fifo and give read, write and execute rights to the owner
    match unistd::mkfifo(&r2_fifo, stat::Mode::S_IRWXU) {
        Ok(_) => {
            info!("created {:?}", r2_fifo);
            assert!(std::path::Path::new(&r2_fifo).exists());
        }
        Err(err) => bail!("Error creating read 2 fifo: {}", err),
    }

    let r1_fifo_clone = r1_fifo.clone();
    let r2_fifo_clone = r2_fifo.clone();

    let join_handle: thread::JoinHandle<Result<()>> = thread::spawn(move || {
        // this is unused, but the move is made so that the tmp_dir
        // lifetime is extended and the directory stays around for
        // the duration of this thread.
        let _local_tmpdir = tmp_dir;
        xform_read_pairs_to_file(geo_re, r1, r2, r1_fifo_clone, r2_fifo_clone)
    });

    Ok(FifoXFormData {
        r1_fifo,
        r2_fifo,
        join_handle,
    })
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
