//! End to end conversions of between biological file formats.

use crate::{
    config::{SnvToolConfig, SvToolConfig},
    error::{Error, Result},
    pubtypes::{self, Breakend, Breakpoint, BreakpointBedpe},
    vcfutils,
};
use std::{collections::HashMap, io::Write, path::Path};

// Define Fundamental Structs

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct MutationConversionStats {
    pub total_records: u64,
    pub pass_records: u64,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct BreakendConversionStats {
    pub total_breakends_before_filtering: u64,
    pub total_breakends_after_filtering: u64,
    pub single_breakends_after_filtering: u64,
    pub paired_breakends_after_filtering: u64,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct BreakpointConversionStats {
    pub total_breakends_before_filtering: u64,
    pub total_breakends_after_filtering: u64,
    pub proper_breakpoints: u64,
    pub single_breakends: u64,
    pub unmatched_breakends: u64,
}
impl MutationConversionStats {
    pub fn print_summary(&self, mut writer: impl Write) -> std::io::Result<()> {
        writeln!(writer, "\n=============\nSummary\n=============")?;
        writeln!(
            writer,
            "{} total mutations in input file",
            self.total_records
        )?;
        writeln!(
            writer,
            "{} mutations after PASS-filtering",
            self.pass_records
        )?;
        writeln!(
            writer,
            "\nThe {} pass mutations were converted",
            self.pass_records
        )?;
        Ok(())
    }
}

impl BreakendConversionStats {
    pub fn print_summary(&self, mut writer: impl Write) -> std::io::Result<()> {
        writeln!(writer, "\n=============\nSummary\n=============")?;
        writeln!(
            writer,
            "{} total breakends in input file",
            self.total_breakends_before_filtering
        )?;
        writeln!(
            writer,
            "{} breakends after PASS-filtering",
            self.total_breakends_after_filtering
        )?;
        writeln!(
            writer,
            "{}/{} are single breakends",
            self.single_breakends_after_filtering, self.total_breakends_after_filtering,
        )?;
        writeln!(
            writer,
            "{}/{} breakends are paired (MATEID field present)",
            self.paired_breakends_after_filtering, self.total_breakends_after_filtering,
        )?;
        writeln!(
            writer,
            "Note we have not verified that the MATEIDs match IDs present in VCF"
        )?;
        writeln!(
            writer,
            "\nThe {} pass breakends have been converted",
            self.total_breakends_after_filtering
        )?;

        Ok(())
    }
}

impl BreakpointConversionStats {
    pub fn print_summary(&self, mut writer: impl Write) -> std::io::Result<()> {
        writeln!(writer, "\n=============\nSummary\n=============")?;
        writeln!(
            writer,
            "{} total breakends in input file",
            self.total_breakends_before_filtering
        )?;
        writeln!(
            writer,
            "{} breakends after PASS-filtering",
            self.total_breakends_after_filtering
        )?;
        writeln!(
            writer,
            "{} proper breakpoints (pairs of breakends that both survived filtering)",
            self.proper_breakpoints
        )?;
        writeln!(
            writer,
            "{} single breakends (no MATEID)",
            self.single_breakends
        )?;
        writeln!(
            writer,
            "{} unmatched breakends (had MATEID but mate not in VCF after filtering)",
            self.unmatched_breakends
        )?;
        writeln!(
            writer,
            "\nThe {} proper breakpoints were converted",
            self.proper_breakpoints
        )?;

        Ok(())
    }
}

/// Convert Structural Variant VCF to a BEDPE-like TSV format file
///
/// The output contains one row per paired [`Breakpoint`]. Single breakends (no MATEID) and
/// unmatched breakends (has MATEID but mate not found in VCF after filtering) are intentionally excluded
/// from this bedpe file.
///
/// This function will filter for PASS variants and sort breakends in each breakpoint based on POS
/// (smaller first) or by order in VCF stream if breakends are on different chromosomes.
///
/// Outputs a BEDPE-like format data with one row per breakpoint including the following columns:
///
/// 1. chrom1: Chromosome of one side of first breakend in pair.
/// 2. start1: Zero-based starting position of the lower confidence interval of first breakend.
/// 3. end1: One-based end position of the upper confidence interval of first breakend.
/// 4. chrom2: Chromosome of second breakend in pair.
/// 5. start2: Zero-based starting position of the lower confidence interval of second breakend in pair.
/// 6. end2: One-based end position of the upper confidence interval of second breakend in pair.
/// 7. name: Breakpoint identifier.
/// 8. score: quality score (from first breakpoint in svcf)
/// 9. strand1: strand of the first breakend in pair
/// 10. strand2: strand for the second breakend in pair
///
/// Plus additional columns: (downstream tools like bedtools allow any number of additional columns - these will just be passed-through)
///
/// 11. vaf1: purity adjusted VAF of first breakend in pair (e.g. from PURPLE_VAF info field if `--from purple`)
/// 12. vaf2: purity adjusted VAF of second breakend in pair (e.g. from PURPLE_VAF info field if `--from purple`)
/// 13. pos1: Zero-based position of first breakend in pair (derived from POS field).
/// 14. pos2: Zero-based position of second breakend in pair (derived from POS field).
pub fn convert_svcf_to_bedpe(vcf: &Path, options: SvToolConfig) -> Result<()> {
    let stdout = std::io::stdout();
    let stderr = std::io::stderr();
    let stats = write_svcf_as_bedpe(vcf, options, stdout.lock())?;
    stats
        .print_summary(stderr.lock())
        .map_err(|source| Error::write("stats", source))?;
    Ok(())
}

/// Write a structural variant VCF as BEDPE-like TSV.
///
/// The output contains one row per paired [`Breakpoint`]. Single breakends (no MATEID) and
/// unmatched breakends (has MATEID but mate not found in VCF after filtering) are intentionally
/// excluded.
///
/// See [`convert_svcf_to_bedpe`] docs for full breakdown of output filetype
pub fn write_svcf_as_bedpe<W: Write>(
    vcf: &Path,
    options: SvToolConfig,
    writer: W,
) -> Result<BreakpointConversionStats> {
    let mut reader = vcfutils::build_vcf_reader(vcf)?;
    let header = vcfutils::read_vcf_header(&mut reader, vcf)?;

    let mut writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(writer);

    let mut stats = BreakpointConversionStats::default();
    let mut waitlist: HashMap<String, Breakend> = HashMap::new();

    for result in reader.records() {
        stats.total_breakends_before_filtering += 1;
        let record = result.map_err(|source| Error::parse_vcf_record(vcf, source))?;

        if !vcfutils::is_pass(&record, &header)? {
            continue;
        }

        stats.total_breakends_after_filtering += 1;
        let breakend = crate::vcfutils::record_to_breakend(&record, &header, &options.vaf_field)?;

        let Some(mateid) = &breakend.mateid else {
            stats.single_breakends += 1;
            continue;
        };

        let Some(mate_breakend) = waitlist.remove(mateid) else {
            waitlist.insert(breakend.id.clone(), breakend);
            continue;
        };

        let breakpoint: Breakpoint = pubtypes::breakpoint_from_vcf_pair(mate_breakend, breakend);
        let bedpe: BreakpointBedpe = pubtypes::breakpoint_to_breakpoint_bedpe(breakpoint);

        writer
            .serialize(bedpe)
            .map_err(|source| Error::write("bedpe", source))?;

        stats.proper_breakpoints += 1;
    }

    writer
        .flush()
        .map_err(|source| Error::flush("bedpe", source))?;

    stats.unmatched_breakends = waitlist.len() as u64;
    Ok(stats)
}

/// Outputs a TSV with one row per breakend
/// (each side paired breakpoints will have their own row).
///
/// Will filter for PASS variants only
///
/// Columns include:
///
/// 1. chromosome: Chromosome of breakend.
/// 2. position: 1-based position of breakend as described by POS column in vcf.
/// 3. vaf: Purity adjusted variant allele frequency supporting breakend (e.g. from PURPLE_VAF info field if `--from purple`).
/// 4. id: id of breakend.
/// 5. mateid: id of mate (set to `.` if single breakend)
/// 6. qual: quality of breakend.
pub fn convert_svcf_to_breakend_tsv(vcf: &Path, options: SvToolConfig) -> Result<()> {
    let stdout = std::io::stdout();
    let stderr = std::io::stderr();

    let stats = write_svcf_as_breakend_tsv(vcf, options, stdout.lock())?;
    stats
        .print_summary(stderr.lock())
        .map_err(|source| Error::write("stats", source))?;
    Ok(())
}

/// Write a structural variant VCF as breakend TSV.
pub fn write_svcf_as_breakend_tsv<W: Write>(
    vcf: &Path,
    options: SvToolConfig,
    writer: W,
) -> Result<BreakendConversionStats> {
    let mut reader = vcfutils::build_vcf_reader(vcf)?;
    let header = vcfutils::read_vcf_header(&mut reader, vcf)?;

    let mut writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(writer);

    let mut stats = BreakendConversionStats::default();

    for result in reader.records() {
        stats.total_breakends_before_filtering += 1;
        let record = result.map_err(|source| Error::parse_vcf_record(vcf, source))?;

        if !vcfutils::is_pass(&record, &header)? {
            continue;
        }

        stats.total_breakends_after_filtering += 1;
        let breakend = crate::vcfutils::record_to_breakend(&record, &header, &options.vaf_field)?;
        let simplebreakend = crate::pubtypes::breakend_to_simple_breakend(&breakend);

        if pubtypes::breakend_is_single(&breakend) {
            stats.single_breakends_after_filtering += 1
        } else {
            stats.paired_breakends_after_filtering += 1
        }

        writer
            .serialize(simplebreakend)
            .map_err(|source| Error::write("breakend-tsv", source))?;
    }

    writer
        .flush()
        .map_err(|source| Error::flush("breakend-tsv", source))?;

    Ok(stats)
}

/// Output a TSV with one row per mutation
///
/// vaf_field should represent the INFO field of VCF describing purity adjusted variant allelic
/// frequency
///
/// This function will return an error if input VCF has multiallelic sites - the error message will
/// tell user to normalise with bcftools norm to split these multiallelic sites.
///
///
/// 1. chrom: Chromosome of breakend.
/// 2. pos: one-based position of variant
/// 3. ref: reference sequence
/// 4. alt: alternate sequence
/// 5. vaf: variant allele frequency of tumour sample (adjusted for purity). Must be an INFO field (not FORMAT).
pub fn convert_snv_vcf_to_tsv(vcf: &Path, options: SnvToolConfig) -> Result<()> {
    let stdout = std::io::stdout();
    let stderr = std::io::stdout();
    let stats = write_snv_vcf_as_tsv(vcf, options, stdout.lock())?;
    stats
        .print_summary(stderr.lock())
        .map_err(|source| Error::write("stats", source))?;
    Ok(())
}

/// Write a SNV/MNV/INDEL VCF as mutation TSV.
///
/// This function returns an error for multiallelic sites. Normalize with `bcftools norm` before
/// conversion if multiallelic records are present.
///
/// See [`convert_snv_vcf_to_tsv`] for a full description of output filetype
pub fn write_snv_vcf_as_tsv<W: Write>(
    vcf: &Path,
    options: SnvToolConfig,
    writer: W,
) -> Result<MutationConversionStats> {
    let mut reader = vcfutils::build_vcf_reader(vcf)?;
    let header = vcfutils::read_vcf_header(&mut reader, vcf)?;

    let mut writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(writer);

    let mut stats = MutationConversionStats::default();

    for result in reader.records() {
        stats.total_records += 1;
        let record = result.map_err(|source| Error::parse_vcf_record(vcf, source))?;

        if !vcfutils::is_pass(&record, &header)? {
            continue;
        }

        stats.pass_records += 1;
        let mutation = crate::vcfutils::record_to_mutation(&record, &header, &options.vaf_field)?;

        writer
            .serialize(mutation)
            .map_err(|source| Error::write("snv-tsv", source))?;
    }

    writer
        .flush()
        .map_err(|source| Error::flush("snv-tsv", source))?;

    Ok(stats)
}
