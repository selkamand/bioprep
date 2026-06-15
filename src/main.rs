use anyhow::Result;
use bioprep::conversions::{svcf_to_bedpe, svcf_to_breakend_tsv};
use clap::{Parser, Subcommand, ValueEnum};
use std::{fmt, path::PathBuf};

#[derive(Debug, Clone, Copy, ValueEnum)]
enum SvcfTypes {
    Purple,
}

impl fmt::Display for SvcfTypes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SvcfTypes::Purple => write!(f, "Purple"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, ValueEnum)]
enum SvOutputTypes {
    Bedpe,
    BreakendTsv,
}
#[derive(Debug, Clone, Copy, ValueEnum)]
enum SnvVcfTypes {
    Purple,
}

impl fmt::Display for SnvVcfTypes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SnvVcfTypes::Purple => write!(f, "Purple"),
        }
    }
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
        #[arg(long, value_enum, default_value_t = SvcfTypes::Purple)]
        from: SvcfTypes,

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
        #[arg(long, value_enum, default_value_t = SnvVcfTypes::Purple)]
        from: SnvVcfTypes,

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
                SvcfTypes::Purple => "PURPLE_AF",
            };

            match to {
                SvOutputTypes::Bedpe => svcf_to_bedpe(&svcf, vaf_field)?,
                SvOutputTypes::BreakendTsv => svcf_to_breakend_tsv(&svcf, vaf_field)?,
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
