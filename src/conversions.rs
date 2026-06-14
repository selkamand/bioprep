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
/// from this bedpe file
///
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

        // If mate has been seen, create a breakpoint struct
        // TODO: Choose which breakend is first or second based on Pos (if they're on the same
        // chromosome)
        let breakpoint = Breakpoint {
            first: breakend,
            second: mate_breakend,
        };

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
