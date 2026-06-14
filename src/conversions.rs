// End to end conversions of file formats from one to another
// Where inputs are filetypes and outputs are writers

use crate::error::Result;

fn svcf_to_bedpe(vcf: &Path, vaf_field: &str) -> Result<()> {
    // Get serialised version fo sv VCF
    let structural_variants = svcf_to_structural_variants(vcf, vaf_field)?;

    // Print to STDOUT
    let stdout = std::io::stdout();
    structural_variants.write_bedpe_tsv(&stdout)?;

    Ok(())
}
