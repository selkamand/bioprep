//! Assorted utilities for working with VCF files

use crate::error::Error;
use crate::error::Result;
use crate::pubtypes::Breakend;
use crate::pubtypes::Mutation;
use crate::pubtypes::Strand;
use noodles::vcf;
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record::Filters;
use noodles::vcf::variant::record::Ids;
use noodles::vcf::variant::record::info::field;
use noodles::vcf::variant::record::info::field::value::Array;
use std::io::BufRead;
use std::path::Path;

// Extract MATEID info field
pub fn parse_mate_id(record: &vcf::Record, header: &vcf::Header) -> Result<Option<String>> {
    let info = record.info();
    let mateid_key: &str = "MATEID";
    let Some(value_result) = info.get(header, "MATEID") else {
        return Ok(None);
    };

    let Some(value) = value_result.map_err(|error| {
        Error::invalid_info("MATEID", format!("failed to parse field: {error}"))
    })?
    else {
        return Ok(None);
    };

    let mate_id = match value {
        field::Value::Integer(n) => n.to_string(),
        field::Value::Float(n) => n.to_string(),
        field::Value::Character(c) => c.to_string(),
        field::Value::String(s) => s.to_string(),
        field::Value::Flag => {
            return Err(Error::invalid_info(
                "MATEID",
                "field is a flag type, which cannot be coerced to a string",
            ));
        }
        field::Value::Array(arr) => parse_one_string_from_array(arr, mateid_key)?,
    };

    if mate_id == "." {
        Ok(None)
    } else {
        Ok(Some(mate_id))
    }
}

// Get stock standard ID column from VCF
pub(crate) fn parse_id(record: &vcf::Record) -> Result<String> {
    let ids = record.ids();

    if ids.len() > 1 {
        return Err(Error::invalid_record(
            "multiple IDs found in a single ID column of SV VCF",
        ));
    }

    if ids.len() == 0 {
        return Err(Error::invalid_record(format!(
            "SV entry lacks a value in the ID column [{}]",
            extract_chrom_pos_ref_alt_neverfail(record),
        )));
    }

    let Some(id_str) = ids.iter().next() else {
        return Err(Error::invalid_record(format!(
            "SV entry lacks a value in the ID column [{}]",
            extract_chrom_pos_ref_alt_neverfail(record),
        )));
    };

    let id = id_str.to_owned();

    if id == "." {
        return Err(Error::invalid_record(format!(
            "SV entry has a missing ID value '.' [{}]",
            extract_chrom_pos_ref_alt_neverfail(record),
        )));
    }

    Ok(id)
}

pub(crate) fn extract_chrom_pos_ref_alt_neverfail(record: &vcf::Record) -> String {
    let pos = match record.variant_start() {
        Some(rp) => match rp {
            Ok(pos) => &pos.to_string(),
            Err(_err) => "?",
        },
        None => "?",
    };

    let refbase = record.reference_bases();
    let altbases = record.alternate_bases();
    let altbase = altbases.iter().next().unwrap_or(Ok("?")).unwrap_or("?");

    format!(
        "{}:{} {} > {}",
        record.reference_sequence_name(),
        pos,
        refbase,
        altbase,
    )
}

pub(crate) fn parse_cipos(record: &vcf::Record, header: &vcf::Header) -> Result<(i32, i32)> {
    let info = record.info();

    let Some(value_result) = info.get(header, "CIPOS") else {
        return Err(Error::invalid_info("CIPOS", "field not found"));
    };

    let Some(value) = value_result
        .map_err(|error| Error::invalid_info("CIPOS", format!("failed to parse field: {error}")))?
    else {
        return Err(Error::invalid_info("CIPOS", "field has no value"));
    };

    let array = match value {
        field::Value::Array(array) => array,
        _ => {
            return Err(Error::invalid_info(
                "CIPOS",
                format!("expected array, got {value:#?}"),
            ));
        }
    };

    let array_int = match array {
        field::value::Array::Integer(values) => values,
        _ => {
            return Err(Error::invalid_info(
                "CIPOS",
                format!("expected integer array, got {array:#?}"),
            ));
        }
    };

    if array_int.len() != 2 {
        return Err(Error::invalid_info(
            "CIPOS",
            format!("expected 2 numbers, found {}", array_int.len()),
        ));
    }

    let mut iter = array_int.iter();

    let lo = iter
        .next()
        .ok_or_else(|| Error::invalid_info("CIPOS", "expected first integer value"))?
        .map_err(|error| {
            Error::invalid_info("CIPOS", format!("failed to parse first value: {error}"))
        })?
        .ok_or_else(|| Error::invalid_info("CIPOS", "first integer value is missing"))?;

    let hi = iter
        .next()
        .ok_or_else(|| Error::invalid_info("CIPOS", "expected second integer value"))?
        .map_err(|error| {
            Error::invalid_info("CIPOS", format!("failed to parse second value: {error}"))
        })?
        .ok_or_else(|| Error::invalid_info("CIPOS", "second integer value is missing"))?;

    Ok((lo, hi))
}

fn parse_vaf(record: &vcf::Record, header: &vcf::Header, vaf_field: &str) -> Result<f32> {
    let info = record.info();

    let Some(value_result) = info.get(header, vaf_field) else {
        return Err(Error::invalid_info(vaf_field, "field not found"));
    };

    let Some(value) = value_result.map_err(|error| {
        Error::invalid_info(vaf_field, format!("failed to parse field: {error}"))
    })?
    else {
        return Err(Error::invalid_info(
            vaf_field,
            "field is present but has no value",
        ));
    };

    match value {
        field::Value::Float(vaf) => Ok(vaf),
        field::Value::Array(array) => parse_first_float_from_array(array, vaf_field),
        field::Value::Integer(_) => Err(Error::invalid_info(
            vaf_field,
            "expected Float, got Integer",
        )),
        field::Value::Flag => Err(Error::invalid_info(vaf_field, "expected Float, got Flag")),
        field::Value::Character(_) => Err(Error::invalid_info(
            vaf_field,
            "expected Float, got Character",
        )),
        field::Value::String(_) => {
            Err(Error::invalid_info(vaf_field, "expected Float, got String"))
        }
    }
}
fn parse_first_float_from_array(array: Array<'_>, field_name: &str) -> Result<f32> {
    let Array::Float(values) = array else {
        return Err(Error::invalid_info(field_name, "expected a Float array"));
    };

    let mut iter = values.iter();

    let value_result = iter
        .next()
        .ok_or_else(|| Error::invalid_info(field_name, "expected one float value"))?;

    let value = value_result
        .map_err(|error| {
            Error::invalid_info(field_name, format!("failed to parse float value: {error}"))
        })?
        .ok_or_else(|| Error::invalid_info(field_name, "contains a missing float value"))?;

    Ok(value)
}

pub(crate) fn parse_one_string_from_array(array: Array<'_>, field_name: &str) -> Result<String> {
    let Array::String(values) = array else {
        return Err(Error::invalid_info(field_name, "expected a String array"));
    };

    if values.len() != 1 {
        return Err(Error::invalid_info(
            field_name,
            format!("expected exactly one string value, got {}", values.len()),
        ));
    }

    let mut iter = values.iter();

    let value_result = iter
        .next()
        .ok_or_else(|| Error::invalid_info(field_name, "expected one string value"))?;

    let value = value_result
        .map_err(|error| {
            Error::invalid_info(field_name, format!("failed to parse string value: {error}"))
        })?
        .ok_or_else(|| Error::invalid_info(field_name, "contains a missing string value"))?;

    Ok(value.to_string())
}

// ALT pattern             local breakend strand
// ------------------------------------------------
// s[chr:pos[              +
// s]chr:pos]              +
// ]chr:pos]s              -
// [chr:pos[s              -
// s.                      +
// .s                      -

pub(crate) fn alt_to_strand(alt: String) -> Result<Strand> {
    if alt.ends_with('[') {
        return Ok(Strand::Plus);
    }

    if alt.ends_with(']') {
        return Ok(Strand::Plus);
    }

    if alt.ends_with('.') {
        return Ok(Strand::Plus);
    }

    if alt.starts_with(']') {
        return Ok(Strand::Minus);
    }

    if alt.starts_with('[') {
        return Ok(Strand::Minus);
    }

    if alt.starts_with('.') {
        return Ok(Strand::Minus);
    }

    Err(Error::InvalidAlt { alt })
}

// Get first alternate base as string. If anything goes wrong or alt is empty return None
fn get_first_alt_as_string(record: &vcf::Record) -> Option<String> {
    let alternative_bases = record.alternate_bases();
    let first = alternative_bases.iter().next()?;

    let unwrapped = match first {
        Ok(b) => b,
        Err(_) => return None,
    };

    Some(unwrapped.to_string())
}

/// Check that variant is pass
pub(crate) fn is_pass(record: &vcf::Record, header: &vcf::Header) -> Result<bool> {
    match record.filters().is_pass(header) {
        Ok(b) => Ok(b),
        Err(err) => Err(Error::FilterStatus {
            record_id: parse_id(record).unwrap_or("ID not found".to_string()),
            source: err,
        }),
    }
}

/// Infer which compressiion method to use when reading a VCF file based on its extension
///
/// If extension is 'gz' assume Bgzf compression
/// Otherwise assume its uncompressed
pub(crate) fn infer_compression_method_based_on_extension(
    vcf: &Path,
) -> vcf::io::CompressionMethod {
    match vcf.extension().is_some_and(|ext| ext == "gz") {
        true => vcf::io::CompressionMethod::Bgzf,
        false => vcf::io::CompressionMethod::None,
    }
}

pub(crate) fn ensure_vcf_exists(vcf: &Path) -> Result<()> {
    if !vcf.exists() {
        return Err(Error::FileNotFound(vcf.to_owned()));
    }
    Ok(())
}

// Build a VCF reader (inferring the appropriate)
pub(crate) fn build_vcf_reader(vcf: &Path) -> Result<vcf::io::reader::Reader<Box<dyn BufRead>>> {
    // Check files exist
    ensure_vcf_exists(vcf)?;

    // Identify compression method
    let compression_method = infer_compression_method_based_on_extension(vcf);

    // Create reader
    let reader = vcf::io::reader::Builder::default()
        .set_compression_method(compression_method)
        .build_from_path(vcf)
        .map_err(|source| Error::ReadVcf {
            path: vcf.to_owned(),
            source,
        })?;

    Ok(reader)
}

pub(crate) fn read_vcf_header(
    reader: &mut vcf::io::reader::Reader<Box<dyn BufRead>>,
    path: &Path,
) -> Result<vcf::Header> {
    match reader.read_header() {
        Ok(header) => Ok(header),
        Err(source) => Err(Error::ParseVcfHeader {
            path: path.to_owned(),
            source,
        }),
    }
}

// Convert Record to Small mutation type
pub fn record_to_mutation(
    record: &vcf::Record,
    header: &vcf::Header,
    vaf_field: &str,
) -> Result<Mutation> {
    // Grab Chrom
    let chrom = record.reference_sequence_name();
    let pos = parse_position(record)?;
    let reference = record.reference_bases();
    let alternative_alleles = parse_alternative_allele_but_error_if_multiallelic(record)?;
    let vaf = parse_vaf(record, header, vaf_field)?;

    Ok(Mutation {
        chrom: chrom.to_owned(),
        pos,
        reference: reference.to_owned(),
        alternative: alternative_alleles,
        vaf,
    })
}

/// Parse a 1-based position from a record or return an error if the position is invalid.
pub(crate) fn parse_position(record: &vcf::Record) -> Result<u64> {
    let invalid_position = |message: String| Error::InvalidVariantRecord {
        variant: extract_chrom_pos_ref_alt_neverfail(record),
        message,
    };

    let Some(pos_result) = record.variant_start() else {
        return Err(invalid_position("no valid position found".to_owned()));
    };

    let position = pos_result.map_err(|source| {
        invalid_position(format!("error trying to access position field: {source}"))
    })?;

    let pos: u64 = position
        .get()
        .try_into()
        .map_err(|_| invalid_position("position could not be converted to u64".to_owned()))?;

    Ok(pos)
}

// let invalid_position = |message: String| Error::InvalidVariantRecord {
//         variant: extract_chrom_pos_ref_alt_neverfail(record),
//         message,
//     };
//

fn parse_alternative_allele_but_error_if_multiallelic(record: &vcf::Record) -> Result<String> {
    let altbases = record.alternate_bases();
    let first_or_dot = altbases.iter().next().unwrap_or(Ok(".")).unwrap_or(".");

    if altbases.len() > 1 {
        return Err(Error::multiple_alternative_alleles(
            extract_chrom_pos_ref_alt_neverfail(record),
            "Multiallelics are not supported. Please remove by running `bcftools norm` and retry",
        ));
    }

    return Ok(first_or_dot.to_owned());
}

/// Parse a noodles SV vcf record into a Breakend
pub fn record_to_breakend(
    record: &vcf::Record,
    header: &vcf::Header,
    vaf_field: &str,
) -> Result<Breakend> {
    // Grab ID
    let breakend_id = crate::vcfutils::parse_id(record)?;

    // Fetch MATEID (if none: Single breakend
    let mateid = crate::vcfutils::parse_mate_id(record, header).map_err(|error| {
        Error::invalid_variant_record(
            breakend_id.clone(),
            format!("failed to parse MATEID: {error}"),
        )
    })?;

    // Grab Position
    let Some(pos_res) = record.variant_start() else {
        return Err(Error::invalid_variant_record(
            breakend_id.clone(),
            "failed to get position",
        ));
    };

    let pos_usize = pos_res
        .map_err(|error| {
            Error::invalid_variant_record(
                breakend_id.clone(),
                format!("failed to parse position: {error}"),
            )
        })?
        .get();

    let pos_1based: u64 = pos_usize.try_into()?; // Convert to u64
    let pos: u64 = pos_1based.saturating_sub(1); // Convert to 0-based position

    // Grab position confidence interval CIPOS or if anything goes wrong set as 0,0
    let cipos = parse_cipos(record, header).unwrap_or_default(); // default will
    // be (0_i32, 0_i32). NOTE in PURPLE VCFs CIPOS lower element will be negative

    // Get the absolute value of lower and higher end
    let cipos_low: u64 = cipos.0.unsigned_abs().into();
    let cipos_high: u64 = cipos.1.unsigned_abs().into();

    // Get Start Position by subtracting absolute value of lower CIPOS to 0-based position
    let start = pos.saturating_sub(cipos_low);

    // Get End Position by adding absolute value of upper CIPOS to 1-based position since bedpe
    // end positions are non-inclusive
    let end = pos_1based.saturating_add(cipos_high);

    // Grab VAF
    let vaf = parse_vaf(record, header, vaf_field).map_err(|error| {
        Error::invalid_variant_record(breakend_id.clone(), format!("failed to parse VAF: {error}"))
    })?;

    // Grab QUAL  (or NAN if anything went wrong)
    let qual = record
        .quality_score()
        .unwrap_or(Ok(f32::NAN))
        .unwrap_or(f32::NAN);

    // Grab CHROM
    let chrom = record.reference_sequence_name();

    // Grab ALT (used to infer strand)
    let altbases = record.alternate_bases();
    if altbases.len() != 1 {
        return Err(Error::invalid_variant_record(
            breakend_id.clone(),
            format!(
                "expected a single alternative sequence but found {}",
                altbases.len()
            ),
        ));
    }
    let alt = get_first_alt_as_string(record).ok_or_else(|| {
        Error::invalid_variant_record(
            breakend_id.clone(),
            "failed to pull a valid alternative sequence",
        )
    })?;

    // Infer strand from alt sequence Grab strand
    let strand = alt_to_strand(alt)?;

    // Create Breakend
    Ok(Breakend {
        chrom: chrom.to_string(),
        start,
        end,
        pos,
        id: breakend_id,
        mateid,
        strand,
        qual,
        vaf,
    })
}
