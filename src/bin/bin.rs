use std::fs::File;
use std::io::{BufWriter, Write};

use clap::Parser;

use seq_geom_parser::{FragmentGeomDesc, PiscemGeomDesc};
use seq_geom_xform::FragmentGeomDescExt;

use anyhow::Result;
use needletail::{parse_fastx_file, Sequence};

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Expected input read geometry
    #[arg(short, long)]
    geom: String,

    /// read 1 files
    #[arg(short = '1', long)]
    read1: Vec<String>,

    /// read 2 files
    #[arg(short = '2', long)]
    read2: Vec<String>,

    /// where output r1 should be written
    #[arg(short = 'o', long)]
    out1: String,

    /// where output r2 should be written
    #[arg(short = 'w', long)]
    out2: String,
}

fn process_reads(args: Args) -> Result<()> {
    let gd = args.geom;
    let geo = FragmentGeomDesc::try_from(gd.as_str()).unwrap();

    let mut parsed_records = seq_geom_xform::SeqPair::new();

    let f1 = File::create(&args.out1).expect("Unable to create file");
    let f2 = File::create(&args.out2).expect("Unable to create file");

    let mut stream1 = BufWriter::new(f1);
    let mut stream2 = BufWriter::new(f2);

    match geo.as_regex() {
        Ok(mut geo_re) => {
            println!(
                "geometry as regex = Read1 : {:?}, Read2 : {:?}",
                geo_re.r1_re, geo_re.r2_re
            );

            let simp_desc = geo_re.get_simplified_piscem_description_string();
            println!("simplified description of geometry is {}", simp_desc);

            let simp_desc = geo_re.get_simplified_geo_desc();
            println!("simplified description again {:?}, {:?}", &simp_desc.read1_desc, &simp_desc.read2_desc);

            let pd = PiscemGeomDesc::from_geom_pieces(&simp_desc.read1_desc, &simp_desc.read2_desc);
            println!("simplified piscem description again {:?}, {:?}", &pd.read1_desc, &pd.read2_desc);


            for (filename1, filename2) in args.read1.iter().zip(args.read2.iter()) {
                let mut reader = parse_fastx_file(filename1).expect("valid path/file");
                let mut reader2 = parse_fastx_file(filename2).expect("valid path/file");

                while let (Some(record), Some(record2)) = (reader.next(), reader2.next()) {
                    let seqrec = record.expect("invalid record");
                    let seqrec2 = record2.expect("invalid record");

                    if geo_re.parse_into(seqrec.sequence(), seqrec2.sequence(), &mut parsed_records)
                    {
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
        Err(e) => Err(e),
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    process_reads(args)
}
