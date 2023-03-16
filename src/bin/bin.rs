use std::path::PathBuf;

use clap::Parser;

use seq_geom_parser::FragmentGeomDesc; // PiscemGeomDesc, SalmonSeparateGeomDesc};
use seq_geom_xform::FragmentGeomDescExt;

use anyhow::Result;

use tracing::info;
use tracing_subscriber::filter::LevelFilter;
use tracing_subscriber::{fmt, prelude::*, EnvFilter};

/// Program to convert `complex` sequencing fragment geometries
/// into a simpler (normalized) form.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Expected input read geometry specification
    #[arg(short, long)]
    geom: String,

    /// read 1 files, comma delimited
    #[arg(short = '1', long, value_delimiter = ',')]
    read1: Vec<PathBuf>,

    /// read 2 files, comma delimited
    #[arg(short = '2', long, value_delimiter = ',')]
    read2: Vec<PathBuf>,

    /// where output r1 should be written (currently uncompressed)
    #[arg(short = 'o', long)]
    out1: PathBuf,

    /// where output r2 should be written (currently uncompressed)
    #[arg(short = 'w', long)]
    out2: PathBuf,
}

fn process_reads(args: Args) -> Result<()> {
    let gd = args.geom;
    let geo = FragmentGeomDesc::try_from(gd.as_str()).unwrap();

    match geo.as_regex() {
        Ok(geo_re) => {
            info!(
                "geometry as regex = Read1 : {:?}, Read2 : {:?}",
                geo_re.r1_re, geo_re.r2_re
            );

            let simp_desc = geo_re.get_simplified_description_string();
            info!(
                "description the simplified version of this geometry is {}",
                simp_desc
            );

            /*
            let simp_desc = geo_re.get_simplified_geo_desc();

            let pd = PiscemGeomDesc::from_geom_pieces(&simp_desc.read1_desc, &simp_desc.read2_desc);
            info!(
                "piscem description of simplified geometry {:?}, {:?}",
                &pd.read1_desc, &pd.read2_desc
            );

            let sd = SalmonSeparateGeomDesc::from_geom_pieces(
                &simp_desc.read1_desc,
                &simp_desc.read2_desc,
            );
            info!("salmon description of simplified geometry {:?}", &sd);
            */

            seq_geom_xform::xform_read_pairs_to_file(
                geo_re,
                &args.read1,
                &args.read2,
                args.out1,
                args.out2,
            )?;
            Ok(())
        }
        Err(e) => Err(e),
    }
}

fn main() -> Result<()> {
    tracing_subscriber::registry()
        .with(fmt::layer())
        .with(
            EnvFilter::builder()
                .with_default_directive(LevelFilter::INFO.into())
                .from_env_lossy(),
        )
        .init();

    let args = Args::parse();
    process_reads(args)
}
