//! # seq_geom_xform
//!
//! `seq_geom_xform` is a crate for transforming complex fragment library geometries
//! from single-cell sequencing data into simple fragment library geometries.

use std::fmt;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::thread;

use anyhow::{bail, Context, Result};
use regex::bytes::{CaptureLocations, Regex};
use seq_geom_parser::{FragmentGeomDesc, GeomLen, GeomPiece, NucStr};

use needletail::{parse_fastx_file, Sequence};
use thousands::Separable;
use tracing::info;

use nix::sys::stat;
use nix::unistd;
use tempfile::tempdir;

#[derive(Debug)]
pub struct FragmentRegexDesc {
    pub r1_cginfo: Vec<GeomPiece>,
    pub r2_cginfo: Vec<GeomPiece>,
    /// The regular expression expected to match read 1
    pub r1_re: Regex,
    /// The regular expression expected to match read 1
    pub r2_re: Regex,
    /// The CaptureLocations to store capture group information
    /// for read 1. This is re-used between parsing calls to
    /// increase performance.
    r1_clocs: CaptureLocations,
    /// The CaptureLocations to store capture group information
    /// for read 2. This is re-used between parsing calls to
    /// increase performance.
    r2_clocs: CaptureLocations,
}

#[derive(Debug)]
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

/// The maximum width for a `RangedLength` that can be handled with our
/// current padding scheme.  That is, if we have a piece of geometry like
/// Umi(BoundedRange(x, y)), we must have that y-x <= 4.
const BOUNDED_RANGE_LIMIT: u32 = 4;
/// The padding that we will append to each possible length in a variable
/// length geometry piece.
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
    // process each capture group:
    // we start at 1 here because the first group is always the match of the whole string
    // and doesn't correspond to any *explicit* capture.  That is, if we explicilty
    // asked to capture the whole string, then there will still be 2 capture groups
    // (and they will be identical).  So, here, we want to skip the "trivial"
    // match of the whole string, and iterate over the remaining capture locations.
    for cl in 1..clocs.len() {
        if let Some(g) = clocs.get(cl) {
            outstr.push_str(r.get(g.0..g.1).unwrap());

            match gpieces.get(cl - 1) {
                // if we captured some variable length piece of geometry
                // then we have to apply the appropriate padding so that
                // we can pass the result to a non-variable length parser.
                Some(GeomPiece::Barcode(GeomLen::LenRange(_l, h)))
                | Some(GeomPiece::Umi(GeomLen::LenRange(_l, h)))
                | Some(GeomPiece::ReadSeq(GeomLen::LenRange(_l, h))) => {
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
    true
}

fn get_simplified_piscem_string(geo_pieces: &[GeomPiece]) -> String {
    let mut rep = String::new();
    for gp in geo_pieces {
        match gp {
            GeomPiece::Discard(GeomLen::FixedLen(x)) => {
                rep += &format!("x[{}]", x);
            }
            GeomPiece::Barcode(GeomLen::FixedLen(x)) => {
                rep += &format!("b[{}]", x);
            }
            GeomPiece::Umi(GeomLen::FixedLen(x)) => {
                rep += &format!("u[{}]", x);
            }
            GeomPiece::ReadSeq(GeomLen::FixedLen(x)) => {
                rep += &format!("r[{}]", x);
            }
            // NOTE: the + 1 in the rules below assumes we will
            // only ever have variable width of geometry pieces
            // of at most 4 bases. If we need to every move
            // beyond that, this code will have to be generalized.
            GeomPiece::Discard(GeomLen::LenRange(_l, h)) => {
                rep += &format!("x[{}]", h + 1);
            }
            GeomPiece::Barcode(GeomLen::LenRange(_l, h)) => {
                rep += &format!("b[{}]", h + 1);
            }
            GeomPiece::Umi(GeomLen::LenRange(_l, h)) => {
                rep += &format!("u[{}]", h + 1);
            }
            GeomPiece::ReadSeq(GeomLen::LenRange(_l, h)) => {
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
        GeomPiece::Discard(GeomLen::LenRange(_l, h)) => {
            GeomPiece::Discard(GeomLen::FixedLen(h + 1))
        }
        GeomPiece::Barcode(GeomLen::LenRange(_l, h)) => {
            GeomPiece::Barcode(GeomLen::FixedLen(h + 1))
        }
        GeomPiece::Umi(GeomLen::LenRange(_l, h)) => GeomPiece::Umi(GeomLen::FixedLen(h + 1)),
        GeomPiece::ReadSeq(GeomLen::LenRange(_l, h)) => {
            GeomPiece::ReadSeq(GeomLen::FixedLen(h + 1))
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
        let m1 = self.r1_re.captures_read(&mut self.r1_clocs, r1);
        let m2 = self.r2_re.captures_read(&mut self.r2_clocs, r2);

        // if the overall match was not obtained for
        // both of the reads, then don't attempt extraction.
        if m1.or(m2).is_none() {
            return false;
        }

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
    /// Return a `FragmentRegexDesc` corresponding to the current
    /// `FragmentGeomDesc`.  This function returns a `Result` that is
    /// `Ok(FragmentRegexDesc)` if the `FragmentRegexDesc` could be
    /// succesfully created and an `Err(anyhow::Error)` otherwise.
    fn as_regex(&self) -> Result<FragmentRegexDesc, anyhow::Error>;
}

fn geom_piece_as_regex_string(gp: &GeomPiece) -> Result<(String, Option<GeomPiece>)> {
    let mut rep = String::from("");
    let mut geo = None;
    match gp {
        // single lengths
        GeomPiece::Discard(GeomLen::FixedLen(x)) => {
            rep.push_str(&format!(r#"[ACGTN]{{{}}}"#, x));
            // don't need to capture
        }
        GeomPiece::Barcode(GeomLen::FixedLen(x))
        | GeomPiece::Umi(GeomLen::FixedLen(x))
        | GeomPiece::ReadSeq(GeomLen::FixedLen(x)) => {
            rep.push_str(&format!(r#"([ACGTN]{{{}}})"#, x));
            geo = Some(gp.clone());
        }
        // length ranges
        GeomPiece::Discard(GeomLen::LenRange(l, h)) => {
            if h - l > BOUNDED_RANGE_LIMIT {
                bail!("Bounded range can have variable width at most {} but the current element {:?} has variable width {}.",
                    BOUNDED_RANGE_LIMIT, &gp, h-l);
            }
            rep.push_str(&format!(r#"[ACGTN]{{{},{}}}"#, l, h));
            // don't need to capture
        }
        GeomPiece::Barcode(GeomLen::LenRange(l, h))
        | GeomPiece::Umi(GeomLen::LenRange(l, h))
        | GeomPiece::ReadSeq(GeomLen::LenRange(l, h)) => {
            if h - l > BOUNDED_RANGE_LIMIT {
                bail!("Bounded range can have variable width at most {} but the current element {:?} has variable width {}.",
                    BOUNDED_RANGE_LIMIT, &gp, h-l);
            }
            rep.push_str(&format!(r#"([ACGTN]{{{},{}}})"#, l, h));
            geo = Some(gp.clone());
        }
        // fixed sequence
        GeomPiece::Fixed(NucStr::Seq(s)) => {
            // no caputre group because no need to capture this
            // right now
            rep.push_str(s);
        }
        // unbounded pieces
        GeomPiece::Discard(GeomLen::Unbounded) => {
            rep += r#"[ACGTN]*"#;
        }
        GeomPiece::Barcode(GeomLen::Unbounded)
        | GeomPiece::Umi(GeomLen::Unbounded)
        | GeomPiece::ReadSeq(GeomLen::Unbounded) => {
            rep += r#"([ACGTN]*)"#;
            geo = Some(gp.clone());
        }
    }
    Ok((rep, geo))
}

impl FragmentGeomDescExt for FragmentGeomDesc {
    /// Return a `FragmentRegexDesc` corresponding to the current
    /// `FragmentGeomDesc`.  This function returns a `Result` that is
    /// `Ok(FragmentRegexDesc)` if the `FragmentRegexDesc` could be
    /// succesfully created and an `Err(anyhow::Error)` otherwise.
    fn as_regex(&self) -> Result<FragmentRegexDesc, anyhow::Error> {
        let mut r1_re_str = String::from("^");
        let mut r1_cginfo = Vec::<GeomPiece>::new();
        for geo_piece in &self.read1_desc {
            let (str_piece, geo_len) = geom_piece_as_regex_string(geo_piece)?;
            r1_re_str.push_str(&str_piece);
            if let Some(elem) = geo_len {
                r1_cginfo.push(elem);
            }
        }

        // This seems to lead to a slight performance improvement, but consider if
        // we really want to do this.  This checks if the last GeomPiece in the
        // current description is of fixed length or not.  If so (i.e. if it is a fixed
        // length piece), then we add an unbounded `Discard` GeomPiece to the end followed by the
        // end of string anchor.  This anchoring of the regex (seemingly) makes matching a
        // little bit faster.
        if let Some(geo_piece) = &self.read1_desc.last() {
            if geo_piece.is_fixed_len() {
                let (str_piece, _geo_len) =
                    geom_piece_as_regex_string(&GeomPiece::Discard(GeomLen::Unbounded))?;
                r1_re_str.push_str(&str_piece);
            }
        }
        r1_re_str.push('$');

        let mut r2_re_str = String::from("^");
        let mut r2_cginfo = Vec::<GeomPiece>::new();
        for geo_piece in &self.read2_desc {
            let (str_piece, geo_len) = geom_piece_as_regex_string(geo_piece)?;
            r2_re_str.push_str(&str_piece);
            if let Some(elem) = geo_len {
                r2_cginfo.push(elem);
            }
        }

        // This seems to lead to a slight performance improvement, but consider if
        // we really want to do this.  This checks if the last GeomPiece in the
        // current description is of fixed length or not.  If so (i.e. if it is a fixed
        // length piece), then we add an unbounded `Discard` GeomPiece to the end followed by the
        // end of string anchor.  This anchoring of the regex (seemingly) makes matching a
        // little bit faster.
        if let Some(geo_piece) = &self.read2_desc.last() {
            if geo_piece.is_fixed_len() {
                let (str_piece, _geo_len) =
                    geom_piece_as_regex_string(&GeomPiece::Discard(GeomLen::Unbounded))?;
                r2_re_str.push_str(&str_piece);
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
    pub join_handle: thread::JoinHandle<Result<XformStats>>,
}

/// This struct holds some basic statistics about
/// the transformation of a stream of reads.
#[derive(Debug)]
pub struct XformStats {
    pub total_fragments: u64,
    pub failed_parsing: u64,
}

impl XformStats {
    /// Create a new (empty) XformStats
    pub fn new() -> Self {
        Self {
            total_fragments: 0u64,
            failed_parsing: 0u64,
        }
    }
}

impl Default for XformStats {
    /// Create a default (empty) XformStats
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for XformStats {
    /// Formats and returns the canonical string representation of each type of
    /// `GeomPiece`.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            r#"XformStats {{ 
    total fragments: {},
    fragments failing parsing: {},
    percentage successfully transformed fragments: {:.2},
}}"#,
            self.total_fragments.separate_with_commas(),
            self.failed_parsing.separate_with_commas(),
            if self.total_fragments > 0 {
                1_f64 - ((self.failed_parsing as f64) / (self.total_fragments as f64))
            } else {
                1_f64
            } * 100_f64
        )
    }
}

/// Given input file paths (possibly multiple sets of files) in `r1` and `r2`,
/// read sequence records from these files and transform them in accordance with
/// the `FragmentRegexDesc` provided as `geo_re`.  The transformed records are then
/// written out to `r1_ofile` and `r2_ofile`. Currently all output is written in `FASTA`
/// format, so any quality lines or comment lines (if the input is `FASTQ`) will be
/// dropped.
pub fn xform_read_pairs_to_file(
    mut geo_re: FragmentRegexDesc,
    r1: &[PathBuf],
    r2: &[PathBuf],
    r1_ofile: PathBuf,
    r2_ofile: PathBuf,
) -> Result<XformStats> {
    let f1 = File::create(r1_ofile).expect("Unable to open read 1 file");
    let f2 = File::create(r2_ofile).expect("Unable to open read 2 file");

    let mut stream1 = BufWriter::new(f1);
    let mut stream2 = BufWriter::new(f2);

    let mut xform_stats = XformStats::new();
    let mut parsed_records = SeqPair::new();
    for (filename1, filename2) in r1.iter().zip(r2.iter()) {
        let mut reader = parse_fastx_file(filename1).expect("valid path/file");
        let mut reader2 = parse_fastx_file(filename2).expect("valid path/file");

        while let (Some(record), Some(record2)) = (reader.next(), reader2.next()) {
            xform_stats.total_fragments += 1;
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
            } else {
                xform_stats.failed_parsing += 1;
            }
        }
    }
    Ok(xform_stats)
}

/// Given input file paths (possibly multiple sets of files) in `r1` and `r2`,
/// and `FragmentRegexDesc` `geo_re`, this function returns a `Result<FifoXFormData>`.
/// If succesful the `Ok(FifoXFormData)` will contain the paths to 2 fifos (1 for each
/// list of input read files) as well as a `thread::JoinHandle`.  The `thread::JoinHandle`
/// will be for a spawned thread that will read sequence records from the files in `r1` and `r2`
/// and transform these reads in accordance with the `FragmentRegexDesc` provided as `geo_re`.  
/// The transformed records are then written out to the fifos given in the `FifoXFormData` struct.
/// Currently, all output is written in `FASTA` format, so any quality lines or comment lines
/// (if the input is `FASTQ`) will be dropped.  If an error occurs up to the creation of the
/// spawned thread, then this function returns an `Err(anyhow::Error)`.  The spawned thread
/// itself returns a `Result<()>`.
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

    let tmp_dir = tempdir()?;
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

    // we clone this here because we want to move these into
    // the thread that will do the transformation but we need
    // to retain a copy to pass to the FifoXFormData that we
    // will return.
    let r1_fifo_clone = r1_fifo.clone();
    let r2_fifo_clone = r2_fifo.clone();

    let join_handle: thread::JoinHandle<Result<XformStats>> = thread::spawn(move || {
        let xform_stats = xform_read_pairs_to_file(geo_re, &r1, &r2, r1_fifo_clone, r2_fifo_clone)?;
        // Explicitly check for and propagate any errors encountered in the
        // closing and deleting of the temporary directory.  The directory
        // will be deleted when the handle goes out of scope, but without
        // calling this method, any encountered errors will be silently
        // ignored.
        // see: https://docs.rs/tempfile/latest/tempfile/struct.TempDir.html#method.close
        match tmp_dir.close() {
            Ok(_) => Ok(xform_stats),
            Err(e) => {
                bail!("When closing (deleting) the temp directory, the following error was encountered {:?}", e);
            }
        }
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

    /// This test checks that technical reads from
    /// sciseq v3 can be properly parsed.  This is a set
    /// of the first few reads from SRR7827207.  The tuple
    /// in each instance holds the read itself, whether or not
    /// it whould parse succesfully (those not containing the
    /// anchor should not parse), and the expected length of the
    /// variable-length barcode piece before the anchor.
    ///
    /// The test checks that the reads that should parse do, that
    /// those that should not don't, and that the appropriate
    /// padding characters are added to each variable-length
    /// piece based on it's length.
    #[test]
    fn sciseq3_transforms() {
        let test_technical_reads = vec![
            ("TNGCGCATTCAGAGCGCCACTTTCGGAAGATATTTT", true, 9),
            ("TNTATACCTTCAGAGCGTGAGGATGTCCTAGAGGTT", true, 10),
            ("AGAGATGAATCAGAGCTGTGCCGGGCTAACCTCATT", true, 10),
            ("TGAACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTT", false, 0),
            ("AAACTCCAATCAGAGCTCCGAGACAACCATTGGATT", true, 10),
            ("ACGAGGTTTCTGAGCCGATAAAGTGATGGCCTTTTT", false, 0),
            ("GCTCTTAGTCAGAGCCGTTTTGGGCGACGCCTTTTT", true, 9),
            ("TCCGTATGTCAGAGCGACTGATGTTATAGCAGATTT", true, 9),
            ("TCTCTCCATCAGAGCAAAAGATTCATTCAATCATTC", true, 9),
            ("AGAACTCCTCTGAGCAATGTCGCTTATTCTGAGTTT", false, 0),
            ("AAGTATTGGTCAGAGCTACGCATTACGCAACTCCTT", true, 10),
            ("TGTCCTTATTCAGAGCCCATTTACGCCACGCAGCTC", true, 10),
            ("GCGCTCAATCAGAGCCGGTGGAAAGACTCAAGCCCT", true, 9),
            ("TTCTTAACCTCAGAGCTGACTAGTACTAGAGAGTTT", true, 10),
            ("TCCTCGAGTCAGAGCGTCCTGGCTTAATTATTGTTT", true, 9),
            ("AACTGGCATCAGAGCCCTCTAATGATCGCTTCTTTT", true, 9),
            ("TATGCGATTTCAGAGCGCGGGGGGCCGAGAATCCTT", true, 10),
            ("TAGTTACCTTCAGAGCGGTTTACACATCCGACTATT", true, 10),
            ("CAAGCAACTCAGAGCTTTTTCTTTCGCGGTTGGTTT", true, 9),
            ("ACCGTAGCTCAGAGCGGCCAGTTACGCAACTCCTTT", true, 9),
            ("CCAAGGATTCAGAGCGGCGCGCCATCCATGACTTTT", true, 9),
            ("TTACTAAGTCAGAGCAGAGGGACGCCAGGATCTTTT", true, 9),
            ("GTAGCGATTCAGAGCAGGCGTGATTATAGCAGATTT", true, 9),
            ("AGCAACGATCTGAGCATATTCAGACCGCGCAACCAT", false, 0),
            ("ATTAATGCCTCAGAGCACGGTACAGATCTTACGCTT", true, 10),
        ];

        let geo = FragmentGeomDesc::try_from("1{b[9-10]f[CAGAGC]u[8]b[10]}2{r:}").unwrap();

        let mut sp = SeqPair::new();
        match geo.as_regex() {
            Ok(mut geo_re) => {
                for (tr, should_parse, pref_len) in &test_technical_reads {
                    let r = geo_re.parse_into(tr.as_bytes(), tr.as_bytes(), &mut sp);
                    assert_eq!(r, *should_parse);
                    if r {
                        dbg!("tr = {}, sp = {:?}", &tr, &sp);
                        match pref_len {
                            9 => {
                                assert_eq!(&sp.s1[9..11], VAR_LEN_BC_PADDING[1]);
                            }
                            10 => {
                                assert_eq!(&sp.s1[10..11], VAR_LEN_BC_PADDING[0]);
                            }
                            _ => {
                                panic!("shouldn't happen");
                            }
                        }
                    }
                }
            }
            Err(e) => {
                panic!("{:?} : couldn't parse valid geometry!", e);
            }
        }
    }
}
