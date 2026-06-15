// End to end conversions of file formats from one to another
// Where inputs are filetypes and outputs are writers

use crate::{
    error::{Error, Result},
    pubtypes::{self, Breakend, Breakpoint},
    vcfutils,
};
use std::{
    collections::HashMap,
    io::{BufWriter, Write},
    path::Path,
};

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
pub fn svcf_to_bedpe(vcf: &Path, vaf_field: &str) -> Result<()> {
    // Create Reader to VCF
    let mut reader = vcfutils::build_vcf_reader(vcf)?;
    let header = vcfutils::read_vcf_header(&mut reader)?;

    // Create a buffered writer to stdout
    let stdout = std::io::stdout().lock();
    let mut writer = BufWriter::new(stdout);

    // Write header line
    pubtypes::write_bedpe_header(&mut writer)?;

    // Setup iterators
    let mut n_single_breakends: u32 = 0;
    let mut n_paired_breakends: u32 = 0;

    // Setup a hashmap waitlist which breakends will be stashed in until we find their mate
    let mut waitlist: HashMap<String, Breakend> = HashMap::new();
    for result in reader.records() {
        let record = result.map_err(|source| Error::ParseVcfRecord {
            path: vcf.to_owned(),
            source,
        })?;

        // Skip non-pass variants
        if !vcfutils::is_pass(&record, &header)? {
            continue;
        }

        // Parse Breakend
        let breakend = crate::vcfutils::record_to_breakend(&record, &header, vaf_field)?;

        // Skip the rest of loop if its a single breakend
        let Some(mateid) = &breakend.mateid else {
            n_single_breakends += 1;
            continue;
        };

        // Add breakend to waitlist if mate hasn't been seen before
        // Otherwise pull mate breakend from the waitlist (deleting it from the hashmap to save
        // memory)
        let Some(mate_breakend) = waitlist.remove(mateid) else {
            waitlist.insert(breakend.id.clone(), breakend);
            continue;
        };

        // If mate has been seen on waitlist, create a breakpoint struct
        // But sort which breakend is 'first' (e.g. chrom1/start1/etc) or second in breakpoint
        // struct based on the following rules:
        // (1) if chromosomes are the same, which Pos field is smaller
        // (2) if chromosomes are not the same, 'first' breakend is based on order in VCF
        //  - note 'mate_breakend' always occurs earlier in the VCF stream than 'breakend' because of how
        //  the waitlist stashing works.
        let breakpoint: Breakpoint = pubtypes::breakpoint_from_vcf_pair(mate_breakend, breakend);

        // Write breakpoint to stdout
        crate::pubtypes::write_breakpoint_as_bedpe(&breakpoint, &mut writer)?;

        // Add tally of paired breakends
        n_paired_breakends += 1;
    }

    // Flush buffer to make sure all breapoints are written to stdout
    let _ = writer.flush();

    // Count number of breakends that were lift in waitlist (we never found their mate in the VCF)
    let n_unmatched_breakends = waitlist.len();
    // Check how many single breakends we had left
    eprintln!("Wrote {n_paired_breakends} paired breakends to BEDPE file");
    eprintln!("Found {n_single_breakends} single breakends (no MATEID)");
    eprintln!(
        "Found {n_unmatched_breakends} unmatched breakends (had MATEID but mate not in VCF after filtering)"
    );

    Ok(())
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
pub fn svcf_to_breakend_tsv(vcf: &Path, vaf_field: &str) -> Result<()> {
    // Create Reader to VCF
    let mut reader = vcfutils::build_vcf_reader(vcf)?;
    let header = vcfutils::read_vcf_header(&mut reader)?;

    // Create a buffered writer to stdout
    let stdout = std::io::stdout().lock();
    let mut writer = BufWriter::new(stdout);

    // Write header line
    pubtypes::write_breakend_tsv_header(&mut writer)?;

    for result in reader.records() {
        let record = result.map_err(|source| Error::ParseVcfRecord {
            path: vcf.to_owned(),
            source,
        })?;

        // Skip non-pass variants
        if !vcfutils::is_pass(&record, &header)? {
            continue;
        }

        // Parse Breakend
        let breakend = crate::vcfutils::record_to_breakend(&record, &header, vaf_field)?;

        // Write breakend to stdout
        crate::pubtypes::write_breakend_as_tsv(&breakend, &mut writer)?;
    }

    // Flush buffer to make sure all breapoints are written to stdout
    let _ = writer.flush();

    Ok(())
}
