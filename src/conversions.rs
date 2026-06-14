// End to end conversions of file formats from one to another
// Where inputs are filetypes and outputs are writers

use crate::{breakend::StructuralVariants, error::Result, parsers::svcf_to_structural_variants};
use std::path::Path;

/// Convert Structural Variant VCF
fn svcf_to_bedpe(vcf: &Path, vaf_field: &str) -> Result<()> {
    // Get serialised version for sv VCF
    let structural_variants: StructuralVariants = svcf_to_structural_variants(vcf, vaf_field)?;

    // Print to STDOUT
    let stdout = std::io::stdout();
    structural_variants.write_bedpe_tsv(&stdout)?;

    Ok(())
}
