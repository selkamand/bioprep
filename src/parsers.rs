// Parsing files to public types
use std::path::Path;

use noodles::vcf;
use noodles::vcf::io::CompressionMethod;

use crate::error::{Error, Result};
use crate::pubtypes::StructuralVariants;

/// VCF
pub fn svcf_to_structural_variants(vcf: &Path, vaf_field: &str) -> Result<StructuralVariants> {
    let compression_method = match vcf.extension().is_some_and(|ext| ext == "gz") {
        true => CompressionMethod::Bgzf,
        false => CompressionMethod::None,
    };

    let mut reader = vcf::io::reader::Builder::default()
        .set_compression_method(compression_method)
        .build_from_path(vcf)
        .map_err(|source| Error::ReadVcf {
            path: vcf.to_owned(),
            source,
        })?;

    let header = reader.read_header().map_err(|source| Error::ReadVcf {
        path: vcf.to_owned(),
        source,
    })?;
}
