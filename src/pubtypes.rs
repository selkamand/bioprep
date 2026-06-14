// Public types that support common transformations
// Note we only provide the std impls. Most behaviour
// Is actually implemented as non-associated pure functions

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

/// Genomic breaks represented by the two ends of the genome that got stitched together post-break
struct Breakpoint {
    first: Breakend,
    second: Breakend,
}

/// Basic Strand Enum
#[derive(Debug, Clone, PartialEq, Eq)]
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
