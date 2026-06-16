use crate::error::Error;
use crate::error::Result;

// A Small Variant (SNV / MNV / INDEL)

#[derive(serde::Deserialize, serde::Serialize, Debug, Clone, PartialEq)]
pub struct Mutation {
    pub chrom: String,

    /// A 1-based VCF style position of the variant
    pub pos: u64,

    /// Reference Sequence
    #[serde(rename = "ref")]
    pub reference: String,

    /// Alternative Sequence
    #[serde(rename = "alt")]
    pub alternative: String,

    /// A purity adjusted variant allele frequency
    pub vaf: f32,
}

/// A single structural-variant breakend parsed from a VCF record.
///
/// A breakend represents one side of a structural variant breakpoint. Paired
/// breakends are linked by `mateid`. single breakends have no mate so mateid is set to NULL
///
/// Coordinates are BED-style:
/// - `start` and `pos` are 0-based.
/// - `end` is non-inclusive (1-based).
/// - `start..end` represents the confidence interval around `pos`,
///   derived from the VCF `CIPOS` INFO field.
///
/// Strand is inferred from the alt allele content
#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct Breakend {
    /// Reference sequence or contig name for this breakend, e.g. `"chr3"`.
    pub chrom: String,

    /// Start of the breakend confidence interval.
    ///
    /// This is a 0-based BED-style coordinate. It is derived from the VCF
    /// 1-based `POS` field and the lower bound of the `CIPOS` INFO field.
    pub start: u64,

    /// End of the breakend confidence interval.
    ///
    /// This is a non-inclusive BED-style end coordinate. It is derived from the
    /// VCF 1-based `POS` field and the upper bound of the `CIPOS` INFO field.
    ///
    pub end: u64,

    /// Representative breakend position.
    ///
    /// This is the VCF `POS` field converted from 1-based to 0-based
    /// coordinates. Unlike [`Breakend::start`] and [`Breakend::end`], this is a
    /// single representative position rather than an uncertainty interval.
    pub pos: u64,

    /// VCF `ID` for this breakend record.
    ///
    /// In paired GRIDSS/PURPLE-style VCFs, the mate breakend should refer to
    /// this value in its `MATEID` INFO field.
    pub id: String,

    /// VCF `MATEID` for this breakend record.
    ///
    /// This is [`Some`] for paired breakends and [`None`] for single breakends
    /// or records without a usable `MATEID`.
    pub mateid: Option<String>,

    /// Orientation of this breakend.
    ///
    /// This is inferred from the VCF ALT allele breakend notation.
    pub strand: Strand,

    /// VCF `QUAL` score for this breakend.
    ///
    /// Missing or unparsable quality scores may be represented as `NaN`,
    /// depending on the parser configuration.
    pub qual: f32,

    /// Variant allele fraction for this breakend.
    ///
    /// This is parsed from the configured VAF INFO field, such as
    /// `PURPLE_AF`.
    pub vaf: f32,
}

impl std::fmt::Display for Breakend {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}-{} (pos={}) id={} mateid={:?} strand={} qual={} vaf={}",
            self.chrom,
            self.start,
            self.end,
            self.pos,
            self.id,
            self.mateid,
            self.strand,
            self.qual,
            self.vaf,
        )
    }
}

#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub struct SimpleBreakend {
    /// Reference sequence or contig name for this breakend, e.g. `"chr3"`.
    pub chrom: String,

    /// Representative breakend position.
    ///
    /// 1-based position from VCF POS field.
    pub pos: u64,

    /// Orientation of this breakend.
    ///
    /// This is inferred from the VCF ALT allele breakend notation.
    pub strand: Strand,

    /// Variant allele fraction for this breakend.
    ///
    /// This is parsed from the configured VAF INFO field, such as
    /// `PURPLE_AF`.
    pub vaf: f32,

    /// VCF `ID` for this breakend record.
    ///
    /// In paired GRIDSS/PURPLE-style VCFs, the mate breakend should refer to
    /// this value in its `MATEID` INFO field.
    pub id: String,

    /// VCF `MATEID` for this breakend record.
    ///
    /// This is [`Some`] for paired breakends and [`None`] for single breakends
    /// or records without a usable `MATEID`.
    pub mateid: Option<String>,

    /// VCF `QUAL` score for this breakend.
    ///
    /// Missing or unparsable quality scores may be represented as `NaN`,
    /// depending on the parser configuration.
    pub qual: f32,
}

/// TODO: We should DEFINITELY just add a new function that parses a record directly to this
/// SimpleBreakend structure as that will be much more efficient than building up a full 'Breakend'
/// struct and cloning many of the fields
pub fn breakend_to_simple_breakend(breakend: &Breakend) -> SimpleBreakend {
    SimpleBreakend {
        chrom: breakend.chrom.clone(),
        pos: breakend.pos + 1, // We don't check this because pos was origininally 1-based anyway
        strand: breakend.strand.clone(),
        vaf: breakend.vaf,
        id: breakend.id.clone(),
        mateid: breakend.mateid.clone(),
        qual: breakend.qual,
    }
}

pub fn breakend_is_single(breakend: &Breakend) -> bool {
    breakend.mateid.is_none()
}

fn _trim_trailing_lowercase_char(mut s: String) {
    if s.chars().last().is_some_and(|c| c.is_lowercase()) {
        s.pop();
    }
}

/// A paired structural-variant breakpoint represented in BEDPE-like form.
///
/// The first ten fields follow the common BEDPE column layout:
///
/// `chrom1`, `start1`, `end1`, `chrom2`, `start2`, `end2`,
/// `name`, `score`, `strand1`, `strand2`.
///
/// Additional fields record the variant allele fraction and representative
/// breakend position for each side of the breakpoint:
///
/// `vaf1`, `vaf2`, `pos1`, `pos2`.
///
/// Coordinates are inherited from [`Breakend`]:
///
/// - `start1`, `start2`, `pos1`, and `pos2` are 0-based.
/// - `end1` and `end2` are non-inclusive BED-style end coordinates.
/// - `start..end` represents the confidence interval around each breakend.
///
/// `name` is produced from the two breakend IDs using
/// [`Breakpoint::combined_identifier`].
///
/// `score` is produced using [`Breakpoint::combined_quality_score`].
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct BreakpointBedpe {
    pub chrom1: String,
    pub start1: u64,
    pub end1: u64,
    pub chrom2: String,
    pub start2: u64,
    pub end2: u64,
    pub name: String,
    pub score: f32,
    pub strand1: Strand,
    pub strand2: Strand,
    pub vaf1: f32,
    pub vaf2: f32,
    pub pos1: u64,
    pub pos2: u64,
}

pub fn breakpoint_to_breakpoint_bedpe(breakpoint: Breakpoint) -> BreakpointBedpe {
    let name = breakpoint.combined_identifier(IdMergeStrategy::DotSeparate);
    let score = breakpoint.combined_quality_score(QualityMergeStrategy::First);

    BreakpointBedpe {
        chrom1: breakpoint.first.chrom,
        start1: breakpoint.first.start,
        end1: breakpoint.first.end,
        chrom2: breakpoint.second.chrom,
        start2: breakpoint.second.start,
        end2: breakpoint.second.end,
        name,
        score,
        strand1: breakpoint.first.strand,
        strand2: breakpoint.second.strand,
        vaf1: breakpoint.first.vaf,
        vaf2: breakpoint.second.vaf,
        pos1: breakpoint.first.pos,
        pos2: breakpoint.second.pos,
    }
}
/// Genomic breaks represented by the two ends of the genome that got stitched together post-break
pub struct Breakpoint {
    pub first: Breakend,
    pub second: Breakend,
}

impl Breakpoint {
    /// Create a breakpoint level identifier by concatenating first and second breakpoint ID (separated
    /// by '.')
    pub fn combined_identifier(&self, strategy: IdMergeStrategy) -> String {
        match strategy {
            IdMergeStrategy::DotSeparate => {
                format!("{}.{}", self.first.id, self.second.id)
            }
        }
    }

    /// Create a breakpoint level identifier by concatenating first and second breakpoint ID (separated
    /// by '.')
    pub fn combined_quality_score(&self, strategy: QualityMergeStrategy) -> f32 {
        match strategy {
            QualityMergeStrategy::First => self.first.qual,
            // QualityMergeStrategy::Second => self.second.qual,
            // QualityMergeStrategy::Mean => (self.first.qual + self.second.qual) / 2.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QualityMergeStrategy {
    First,
    // Second,
    // Mean,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IdMergeStrategy {
    DotSeparate,
}

fn later_breakend_should_be_first(earlier: &Breakend, later: &Breakend) -> bool {
    earlier.chrom == later.chrom && later.pos <= earlier.pos
}

/// Build a [`Breakpoint`] from two mated breakends encountered in VCF order.
///
/// `earlier` is the breakend that appeared first in the VCF stream, and `later`
/// is its mate that appeared later in the VCF stream. These names refer only to
/// record order in the input VCF, not to genomic coordinate order.
///
/// BEDPE output order is assigned as follows:
///
/// - if both breakends are on the same chromosome, the breakend with the lower
///   position is emitted as `first`;
/// - if both breakends are on different chromosomes, the original VCF order is
///   preserved, so `earlier` is emitted as `first`.
///
/// This function takes ownership of both breakends and moves them into the
/// returned [`Breakpoint`], avoiding any additional cloning.
pub(crate) fn breakpoint_from_vcf_pair(earlier: Breakend, later: Breakend) -> Breakpoint {
    if later_breakend_should_be_first(&earlier, &later) {
        Breakpoint {
            first: later,
            second: earlier,
        }
    } else {
        Breakpoint {
            first: earlier,
            second: later,
        }
    }
}

// Basic Strand Enum
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub enum Strand {
    Plus,
    Minus,
}

// Display strand
impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Plus => write!(f, "+"),
            Strand::Minus => write!(f, "-"),
        }
    }
}
//
// /// Structural variants parsed from an SV VCF.
// ///
// /// The parsed variants are split into three groups:
// /// - `breakpoints`: complete paired breakends.
// /// - `single_breakends`: breakends with no `MATEID`.
// /// - `unmatched_breakends`: breakends with a `MATEID` whose mate was not found,
// ///   usually because the mate was filtered out or absent from the input VCF.
// pub struct StructuralVariants {
//     /// Complete paired breakpoints that can be written to BEDPE.
//     pub breakpoints: Vec<Breakpoint>,
//
//     /// Breakends with no `MATEID`, representing single breakends.
//     pub single_breakends: Vec<Breakend>,
//
//     /// Breakends with a `MATEID` whose mate was not found in the input VCF.
//     pub unmatched_breakends: Vec<Breakend>,
// }
//
// /// Breakend Type
// enum Breaktype {
//     Paired,
//     PairedWithMateMissing,
//     Single,
// }
