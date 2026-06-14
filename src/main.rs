use anyhow::Result;
use bioprep::breakend::vcf_to_structural_variants;
use clap::{Parser, Subcommand, ValueEnum};
use std::{
    fmt,
    path::{Path, PathBuf},
};

#[derive(Debug, Clone, Copy, ValueEnum)]
enum VcfTypes {
    Purple,
}

impl fmt::Display for VcfTypes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VcfTypes::Purple => write!(f, "Purple"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, ValueEnum)]
enum SvOutputTypes {
    Bedpe,
    BreakendTsv,
}

#[derive(Debug, Clone, PartialEq, Eq, ValueEnum)]
enum SnvOutputTypes {
    Maf,
}

#[derive(Parser)]
#[command(
    version,
    about = "Transform structural variant VCFs to other formats",
    long_about = None
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Convert svcf file to other file formats
    Svcf {
        /// Path to a structural variant vcf
        #[arg(short = 'i', long = "input", value_name = "vcf")]
        svcf: PathBuf,

        /// Input sv vcf filetype
        #[arg(long, value_enum, default_value_t = VcfTypes::Purple)]
        from: VcfTypes,

        /// Output filetype
        #[arg(long, value_enum, value_name = "filetype")]
        to: SvOutputTypes,
    },
    /// Convert SNV/INDEL mutation VCF file to other file formats
    Vcf {
        /// Path to a structural variant vcf
        #[arg(short = 'i', long = "input", value_name = "vcf")]
        vcf: PathBuf,

        /// Input sv vcf filetype
        #[arg(long, value_enum, default_value_t = VcfTypes::Purple)]
        from: VcfTypes,

        /// Output filetype
        #[arg(long, value_enum, value_name = "filetype")]
        to: SnvOutputTypes,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Svcf { svcf, from, to } => {
            let vaf_field = match from {
                VcfTypes::Purple => "PURPLE_AF",
            };

            match to {
                SvOutputTypes::Bedpe => svcf_to_bedpe(&svcf, vaf_field)?,
                SvOutputTypes::BreakendTsv => todo!(),
            };
        }
        Commands::Vcf {
            vcf: _,
            from: _,
            to: _,
        } => {
            todo!("No implementation available for converting vcf to anything");
        }
    };

    Ok(())
}

fn svcf_to_bedpe(vcf: &Path, vaf_field: &str) -> Result<()> {
    // Get serialised version fo sv VCF
    let structural_variants = vcf_to_structural_variants(vcf, vaf_field)?;

    // Print to STDOUT
    let stdout = std::io::stdout();
    structural_variants.write_bedpe_tsv(&stdout)?;

    Ok(())
}
