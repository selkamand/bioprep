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

    /// Compute sets of statistics
    Stats {
        // schemes
        #[command(subcommand)]
        statset: StatSets,
    },

    /// Mutational Signatures
    Signatures {},

    /// Predict biological properties of a tumour
    Predict {
        #[command(subcommand)]
        model: PredictionModels,
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
    /// Tally single base changes into transitions / transversions
    Titv {
        /// Path to a standardised bioprep SNV TSV
        #[arg(short = 'i', long = "input", value_name = "mutation tsv", value_hint = ValueHint::FilePath)]
        snv_tsv: PathBuf,
    },
    /// Tally small mutation types (snv, doublet, mnv, deletion, insertion)
    Smallmuts {
        /// Path to a standardised bioprep SNV TSV
        #[arg(short = 'i', long = "input", value_name = "mutation tsv", value_hint = ValueHint::FilePath)]
        snv_tsv: PathBuf,
    },
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
    },
    /// Tally structural variants into SV32 substitution classes
    Sv32 {
        /// Path to a standardised bioprep BEDPE-like TSV
        #[arg(short = 'i', long = "input", value_name = "bedpe-tsv", value_hint = ValueHint::FilePath)]
        bedpe: PathBuf,
    },
    /// Tally breakpoints into 4 types (Trans, Del, Inv, Tds)
    Breakpoints {
        /// Path to a standardised bioprep BEDPE-like TSV
        #[arg(short = 'i', long = "input", value_name = "bedpe-tsv", value_hint = ValueHint::FilePath)]
        bedpe: PathBuf,
    },
}

#[derive(Subcommand)]
enum PredictionModels {
    /// Predicts whether patient has a P53 dysfunction
    P53detect {
        /// Path to a standardised bioprep SNV TSV
        #[arg(short = 'i', long = "input", value_name = "instability_stats.tsv", value_hint = ValueHint::FilePath)]
        instability_stats: PathBuf,

        // Model weights
        #[arg(short = 'w', long = "weights", value_name = "weights.tsv", value_hint = ValueHint::FilePath)]
        weights: PathBuf,
    },

    /// Predicts whether a haematological cancer patient has an ETV6 disruption based on mutational signature analysis
    ETV6detect {
        /// Path to a standardised bioprep SNV TSV
        #[arg(short = 'i', long = "input", value_name = "features.tsv", value_hint = ValueHint::FilePath)]
        features: PathBuf,

        // Model weights
        #[arg(short = 'w', long = "weights", value_name = "weights.tsv", value_hint = ValueHint::FilePath)]
        weights: PathBuf,
    },
}

#[derive(Subcommand)]
enum StatSets {
    /// Compute measures of genome instability including autosomal LOH and SV burden
    GenomeInstability {
        /// Path to a standardised bioprep SNV TSV
        #[arg(short = 'i', long = "input", value_name = "mutations.tsv", value_hint = ValueHint::FilePath)]
        mutations: PathBuf,
    },
    /// Estimate mitochondrial burden
    Mitochondria {
        /// Path to a idxstats tsv
        #[arg(short = 'i', long = "input", value_name = "idxstats.tsv", value_hint = ValueHint::FilePath)]
        idxstats: PathBuf,

        /// Purity
        #[arg(long = "purity", value_name = "purity")]
        purity: f32,

        /// Ploidy
        #[arg(long = "ploidy", value_name = "ploidy")]
        ploidy: f32,

        /// Mitochondrial chromosome name
        #[arg(long = "mtname", value_name = "mitochondrial contig name")]
        mtname: String,

        /// Autosomal chromosome names (comma separated)
        #[arg(long = "mtname", value_name = "chr1,chr2,...", value_delimiter=',', num_args = 1..)]
        chrnames: Vec<String>,
    },
}
fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Convert { input } => match input {
            ConversionInputCommands::Svcf { svcf, from, to } => {
                let config = bioprep::config::configure_for_sv_tool(from);

                match to {
                    SvOutputTypes::Bedpe => convert_svcf_to_bedpe(svcf.as_path(), config)?,
                    SvOutputTypes::BreakendTsv => {
                        convert_svcf_to_breakend_tsv(svcf.as_path(), config)?
                    }
                };
            }

            ConversionInputCommands::Vcf { vcf, from, to } => {
                let config = configure_for_snv_tool(from);
                match to {
                    SnvOutputTypes::Tsv => convert_snv_vcf_to_tsv(vcf.as_path(), config)?,
                }
            }
        },
        Commands::Tally { scheme } => match scheme {
            ClassificationSchemes::Sbs96 { snv_tsv, reference } => {
                bioprep::tally::tally_sbs96(snv_tsv.as_path(), reference.as_path())?;
            }
            ClassificationSchemes::Sbs6 { snv_tsv } => {
                bioprep::tally::tally_sbs6(&snv_tsv, std::io::stdout().lock())?
            }
            ClassificationSchemes::Sv32 { bedpe: _ } => {
                todo!("No implementation for Sv32 tallying yet")
            }
            ClassificationSchemes::Titv { snv_tsv } => {
                bioprep::tally::tally_titv(snv_tsv.as_path(), std::io::stdout().lock())?;
            }
            ClassificationSchemes::Smallmuts { snv_tsv } => {
                bioprep::tally::tally_small_mutation_types(
                    snv_tsv.as_path(),
                    std::io::stdout().lock(),
                )?;
            }
            ClassificationSchemes::Breakpoints { bedpe } => {
                bioprep::tally::tally_breakpoint_types(&bedpe, std::io::stdout().lock())?
            }
        },
        Commands::Stats { statset } => match statset {
            StatSets::GenomeInstability { mutations: _ } => {
                todo!("No implementation for genome instability stats")
            }
            StatSets::Mitochondria {
                idxstats: _,
                purity: _,
                ploidy: _,
                mtname: _,
                chrnames: _,
            } => {
                todo!("No implementation for mitochondrial burden stats")
            }
        },
        Commands::Predict { model } => match model {
            PredictionModels::P53detect {
                instability_stats: _,
                weights: _,
            } => todo!("P53detect model is not yet implemented"),
            PredictionModels::ETV6detect {
                features: _,
                weights: _,
            } => {
                todo!("ETV6detect model not yet implemented")
            }
        },
        Commands::Signatures {} => todo!("Signature fitting not implemented yet"),
    };

    Ok(())
}
