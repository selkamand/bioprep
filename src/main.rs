use anyhow::Result;
use bioprep::{
    config::{SnvTool, SvTool, configure_for_snv_tool},
    conversions::{convert_snv_vcf_to_tsv, convert_svcf_to_bedpe, convert_svcf_to_breakend_tsv},
};
use clap::{Parser, Subcommand, ValueEnum, ValueHint};
use std::{fmt, path::PathBuf};

#[derive(Debug, Clone, Copy, ValueEnum)]
enum SvcfTool {
    Purple,
}

impl fmt::Display for SvcfTool {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SvcfTool::Purple => write!(f, "Purple"),
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
    Tsv,
}

#[derive(Parser)]
#[command(
    version,
    about = "Prepare biological variant files and tally common mutation classes",
    long_about = None,
    subcommand_required = true,
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
#[command()]
enum Commands {
    /// Perform Common file format conversions
    Convert {
        #[command(subcommand)]
        input: ConversionInputCommands,
    },

    /// Tally variants into mutational signature classification schemes
    Tally {
        #[command(subcommand)]
        scheme: ClassificationSchemes,
    },
}

#[derive(Subcommand)]
enum ConversionInputCommands {
    /// Convert svcf file to other file formats
    Svcf {
        /// Path to a structural variant vcf
        #[arg(short = 'i', long = "input", value_name = "vcf", value_hint = ValueHint::FilePath)]
        svcf: PathBuf,

        /// Input sv vcf filetype
        #[arg(long, value_enum, value_name = "tool")]
        from: SvTool,

        /// Output filetype
        #[arg(long, value_enum, value_name = "filetype")]
        to: SvOutputTypes,
    },
    /// Convert SNV/INDEL mutation VCF file to other file formats
    Vcf {
        /// Path to a SNV/MNV/INDEL variant vcf
        #[arg(short = 'i', long = "input", value_name = "vcf", value_hint = ValueHint::FilePath)]
        vcf: PathBuf,

        /// What type of SNV vcf was supplied
        #[arg(long, value_enum, value_name = "tool")]
        from: SnvTool,

        /// Output filetype
        #[arg(long, value_enum, value_name = "filetype")]
        to: SnvOutputTypes,
    },
}

#[derive(Subcommand)]
enum ClassificationSchemes {
    /// Tally SNVs into SBS96 trinucleotide classes
    Sbs96 {
        /// Path to a standardised bioprep SNV TSV
        #[arg(short = 'i', long = "input", value_name = "mutation tsv", value_hint = ValueHint::FilePath)]
        snv_tsv: PathBuf,

        /// Reference genome FASTA used to fetch trinucleotide context
        #[arg(short = 'r', long = "reference", value_name = "fasta", value_hint = ValueHint::FilePath)]
        reference: PathBuf,
    },
    /// Tally SNVs into SBS6 substitution classes
    Sbs6 {
        /// Path to a standardised bioprep SNV TSV
        #[arg(short = 'i', long = "input", value_name = "mutation tsv", value_hint = ValueHint::FilePath)]
        snv_tsv: PathBuf,

        /// Reference genome FASTA used to fetch trinucleotide context
        #[arg(short = 'r', long = "reference", value_name = "fasta", value_hint = ValueHint::FilePath)]
        reference: PathBuf,
    },
    /// Validate SV32 inputs. Classification rules are not implemented yet.
    Sv32 {
        /// Path to a standardised bioprep BEDPE-like TSV
        #[arg(short = 'i', long = "input", value_name = "bedpe-tsv", value_hint = ValueHint::FilePath)]
        bedpe: PathBuf,

        /// Reference genome FASTA
        #[arg(short = 'r', long = "reference", value_name = "fasta", value_hint = ValueHint::FilePath)]
        reference: PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Convert { input } => match input {
            ConversionInputCommands::Svcf { svcf, from, to } => {
                let config = bioprep::config::configure_for_sv_tool(from);

                match to {
                    SvOutputTypes::Bedpe => convert_svcf_to_bedpe(&svcf, config)?,
                    SvOutputTypes::BreakendTsv => convert_svcf_to_breakend_tsv(&svcf, config)?,
                };
            }

            ConversionInputCommands::Vcf { vcf, from, to } => {
                let config = configure_for_snv_tool(from);
                match to {
                    SnvOutputTypes::Tsv => convert_snv_vcf_to_tsv(&vcf, config)?,
                }
            }
        },
        Commands::Tally { scheme } => match scheme {
            ClassificationSchemes::Sbs96 { snv_tsv, reference } => {
                bioprep::tally::tally_sbs96(&snv_tsv, &reference)?;
            }
            ClassificationSchemes::Sbs6 { snv_tsv, reference } => {
                todo!("No implementation for SBS6 tallying yet")
            }
            ClassificationSchemes::Sv32 { bedpe, reference } => {
                todo!("No implementation for Sv32 tallying yet")
            }
        },
    };

    Ok(())
}
