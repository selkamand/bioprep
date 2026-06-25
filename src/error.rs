use std::{
    collections::HashSet,
    path::{Path, PathBuf},
};

use anyhow::Error as AnyhowError;
use clap::builder::Str;
use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Failed to find file: {0}")]
    FileNotFound(PathBuf),

    #[error("failed to open or read VCF: {path}")]
    ReadVcf {
        path: PathBuf,
        #[source]
        source: AnyhowError,
    },

    #[error("failed to parse VCF header: {path}")]
    ParseVcfHeader {
        path: PathBuf,
        #[source]
        source: AnyhowError,
    },

    #[error("failed to parse VCF record: {path}")]
    ParseVcfRecord {
        path: PathBuf,
        #[source]
        source: AnyhowError,
    },

    #[error("failed to evaluate FILTER/PASS status for VCF record {record_id}")]
    FilterStatus {
        record_id: String,
        #[source]
        source: AnyhowError,
    },

    #[error("failed to write: {filetype}")]
    Write {
        filetype: String,
        #[source]
        source: AnyhowError,
    },

    #[error("failed to flush: {filetype}")]
    Flush {
        filetype: String,
        #[source]
        source: AnyhowError,
    },

    #[error("invalid breakend pairing: {0}")]
    InvalidPairing(String),

    #[error("invalid VCF record: {message}")]
    InvalidRecord { message: String },

    #[error("invalid VCF record {variant}: {message}")]
    InvalidVariantRecord { variant: String, message: String },

    #[error("invalid INFO/{field}: {message}")]
    InvalidInfo { field: String, message: String },

    #[error("invalid ALT allele {alt:?}: failed to infer breakend strand")]
    InvalidAlt { alt: String },

    #[error("Multiple ALT allels in record {variant}: {message}")]
    MultipleAlternativeAlleles { variant: String, message: String },

    #[error("integer conversion failed")]
    IntConversion(#[from] std::num::TryFromIntError),

    #[error("failed to open or read TSV: {path}")]
    ReadTsv {
        path: PathBuf,
        #[source]
        source: AnyhowError,
    },

    #[error("failed to deserialize mutation from TSV file: {path}")]
    DeserializeMutation {
        path: PathBuf,
        #[source]
        source: AnyhowError,
    },

    #[error("failed to deserialize breakpoint from bedpe TSV file: {path}")]
    DeserializeBreakpoint {
        path: PathBuf,
        #[source]
        source: AnyhowError,
    },
    #[error("failed to deserialize copynumber segment from copynumber segment TSV file: {path}")]
    DeserializeCopynumberSegment {
        path: PathBuf,
        #[source]
        source: AnyhowError,
    },
    #[error("failed to deserialize idxstats file: {path}")]
    DeserializeIdxstats {
        path: PathBuf,
        #[source]
        source: AnyhowError,
    },
    #[error("Failed to convert bioprep mutation to seqlib equivalent. Problematic field: {field}")]
    InvalidSequenceForConversion {
        field: String,
        #[source]
        source: AnyhowError,
    },

    #[error("Failed to convert bioprep mutation to seqlib equivalent. Problematic field: position")]
    InvalidPositionForConversion {
        #[source]
        source: AnyhowError,
    },

    #[error("Failed to find any mitochondrial contigs ({mtnames}) in idxstats file: {path}")]
    MissingMitochondrialContig { path: PathBuf, mtnames: String },

    #[error("Failed to find any autosomal contigs ({autosome_names}) in idxstats file: {path}")]
    MissingAutosomalContig {
        path: PathBuf,
        autosome_names: String,
    },
}

impl Error {
    pub(crate) fn invalid_record(message: impl Into<String>) -> Self {
        Self::InvalidRecord {
            message: message.into(),
        }
    }

    pub(crate) fn invalid_variant_record(
        variant: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self::InvalidVariantRecord {
            variant: variant.into(),
            message: message.into(),
        }
    }

    pub(crate) fn write(filetype: impl Into<String>, source: AnyhowError) -> Self {
        Self::Write {
            filetype: filetype.into(),
            source,
        }
    }

    pub(crate) fn invalid_info(field: impl Into<String>, message: impl Into<String>) -> Self {
        Self::InvalidInfo {
            field: field.into(),
            message: message.into(),
        }
    }

    pub(crate) fn parse_vcf_record(path: &Path, source: AnyhowError) -> Self {
        Self::ParseVcfRecord {
            path: path.to_owned(),
            source,
        }
    }

    pub(crate) fn flush(filetype: impl Into<String>, source: AnyhowError) -> Self {
        Self::Flush {
            filetype: filetype.into(),
            source,
        }
    }

    pub(crate) fn multiple_alternative_alleles(
        variant: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self::MultipleAlternativeAlleles {
            variant: variant.into(),
            message: message.into(),
        }
    }
}
