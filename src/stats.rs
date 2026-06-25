//! A module for computing stats from basic filetypes

use serde::{Deserialize, Serialize};

use crate::{
    error::{Error, Result},
    io::{
        create_versioned_tsv_writer, read_copynumber_segements_purple, read_idxstats_tsv,
        serialize_object_to_writer,
    },
    pubtypes::{CopynumberSegments, Idxstats, SegmentClass},
};
use std::{collections::HashSet, io::Write, path::Path};

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct SegmentStats {
    pub prop_autosomal_loh: f32,
    pub prop_autosomal_homdel: f32,
    pub prop_autosomal_het: f32,
    pub total_autosome_size: u64,
}

/// Compute the proportion of the autosomal genome under different types of CN event
/// (hom del VS het vs loh)
pub fn calculate_segment_stats<W: Write>(
    segment_tsv: &Path,
    writer: W,
    autosomes: Option<HashSet<String>>,
) -> Result<()> {
    // Create reader
    let mut rdr = read_copynumber_segements_purple(segment_tsv)?;

    let mut total_autosome_size = 0u64;
    let mut total_loh = 0u64;
    let mut total_homdel = 0u64;
    let mut total_het = 0u64;

    let autosomes = autosomes.unwrap_or(default_autosome_hashset());

    for result in rdr.deserialize() {
        // Deserialize into mutation class from this crate
        let segment: CopynumberSegments =
            result.map_err(|source| Error::DeserializeCopynumberSegment {
                path: segment_tsv.to_owned(),
                source: source.into(),
            })?;

        if autosomes.contains(&segment.chr) {
            let width = segment.width();
            match segment.segment_class() {
                SegmentClass::HomozygousDeletion => total_homdel += width,
                SegmentClass::LossOfHeterozygosity => total_loh += width,
                SegmentClass::Heterozygous => total_het += width,
            }
            total_autosome_size += width
        }
    }

    let stats = SegmentStats {
        prop_autosomal_loh: total_loh as f32 / total_autosome_size as f32,
        prop_autosomal_homdel: total_homdel as f32 / total_autosome_size as f32,
        prop_autosomal_het: total_het as f32 / total_autosome_size as f32,
        total_autosome_size,
    };

    // Create a versioned writer (with tool name and version number already serialized to the
    // writer)
    let writer = create_versioned_tsv_writer(writer, "Stats (Segment)")?;

    // Serialize object to writer
    serialize_object_to_writer(writer, stats, "Stats (Segment)")?;

    Ok(())
}

#[derive(Debug, Clone, PartialEq, Default, Deserialize, Serialize)]
struct MitoStats {
    mt_genomes_per_cell: f32,
    mt_reads: u64,
    mt_reads_per_kb: f32,
    mt_genome_size: u64,
    median_autosomal_reads_per_kb: f32,
    // total_autosome_size: u64,
    // total_autosome_reads: u64,
}

/// Calculate mtdna burden from idxstats file
///
/// @param idxstats_tsv Path to idxstats file produced by `samtools idxstats <bam>`
/// @param mitochondria names of mitochondrial genome in reference genome
/// @param autosomes name of autosomal sequences in reference genome
/// @param verbose verbose mode
///
///
/// # examples
/// path = system.file(package="mitochondrir", "idxstats")
/// parse_idxstats(path)
pub fn calculate_mitochondrial_cn<W: Write>(
    idxstats_tsv: &Path,
    background_ploidy: f32,
    tumour_ploidy: f32,
    tumour_purity: f32,
    writer: W,
    autosomes: Option<HashSet<String>>,
    mitochondria: Option<HashSet<String>>,
) -> Result<()> {
    // Create reader
    let mut rdr = read_idxstats_tsv(idxstats_tsv)?;

    let mut autosomal_reads_per_kb: Vec<f32> = Vec::new();

    let mut mito_contig_seen = false;
    let autosomes = autosomes.unwrap_or(default_autosome_hashset());
    let mitochondria = mitochondria.unwrap_or(default_mitochondrial_hashset());

    let mut stats = MitoStats::default();

    for result in rdr.deserialize() {
        // Deserialize into mutation class from this crate
        let idxstats: Idxstats = result.map_err(|source| Error::DeserializeIdxstats {
            path: idxstats_tsv.to_owned(),
            source: source.into(),
        })?;

        // If autosomal
        if autosomes.contains(&idxstats.contig) {
            // stats.total_autosome_size += idxstats.length;
            // stats.total_autosome_reads += idxstats.n_mapped;

            let reads_per_kb = idxstats.n_mapped as f32 * 1000_f32 / idxstats.length as f32;
            autosomal_reads_per_kb.push(reads_per_kb);
        }

        // If mitochondrial
        if mitochondria.contains(&idxstats.contig) {
            stats.mt_reads += idxstats.n_mapped;
            stats.mt_genome_size += idxstats.length;
            stats.mt_reads_per_kb = stats.mt_reads as f32 * 1000_f32 / stats.mt_genome_size as f32;
            mito_contig_seen = true;
        }
    }

    // Return error if we never saw a mitochondrial contig
    if !mito_contig_seen {
        return Err(Error::MissingMitochondrialContig {
            path: idxstats_tsv.to_owned(),
            mtnames: hashset_to_string(&mitochondria),
        });
    }

    // Calculate median autosomal depth
    let median_autosomal_depth =
        median_f32(&mut autosomal_reads_per_kb).ok_or_else(|| Error::MissingAutosomalContig {
            path: idxstats_tsv.to_owned(),
            autosome_names: hashset_to_string(&autosomes),
        })?;

    stats.median_autosomal_reads_per_kb = median_autosomal_depth;

    // Compute Mitochondrial genomes per cell
    let mt_genoems_per_cell = calculate_mtdna_tumour(
        stats.mt_reads_per_kb,
        stats.median_autosomal_reads_per_kb,
        background_ploidy,
        tumour_ploidy,
        tumour_purity,
        true,
    );
    stats.mt_genomes_per_cell = mt_genoems_per_cell;

    // Create a versioned writer (with tool name and version number already serialized to the
    // writer)
    let writer = create_versioned_tsv_writer(writer, "Stats (MtCN)")?;

    // Serialize object to writer
    serialize_object_to_writer(writer, stats, "Stats (MtCN)")?;

    Ok(())
}

pub fn default_autosome_hashset() -> HashSet<String> {
    crate::constants::DEFAULT_AUTOSOME_NAMES
        .iter()
        .flat_map(|&s| [s.to_string(), format!("chr{s}")])
        .collect()
}

pub fn default_mitochondrial_hashset() -> HashSet<String> {
    crate::constants::DEFAULT_MITOCHONDRIAL_CHROMOSOME_NAMES
        .iter()
        .flat_map(|&s| [s.to_string(), format!("chr{s}")])
        .collect()
}
/// Calculate typical number of mitochondria per cell
///
/// from a bulk whole-genome tumour sample, calculate number of mitochondrial genomes per cell.
/// this function includes options to correct for tumour ploidy and purity,
/// but by default will assume all cells in sample are diploid
///
/// @param mitochondrial_depth depth of coverage of the mitochondrial genome
/// @param autosome_depth depth of coverage of the autosomal genome (do not normalise for ploidy yet! this function performs the normalisation for you)
/// @param background_ploidy expected ploidy of healthy cells in bulk wgs sample (default: 2 for diploid background)
/// @param tumour_ploidy median ploidy of tumour cells in bulk sample (default: 2 for diploid tumours)
/// @param tumour_purity proportion of total cells in sample that are tumour cells (must be between 0 and 1). used to correct for tumour ploidy.
/// @param per_cell return mitochondrial genomes per cell. If FALSE, return mitochondrial genomes per autosome (single copy)
/// @returns number of mitochondrial genomes per cell. if \code{per_cell=FALSE} returns number of mitochondrial genomes per single copy of the autosome.
///
/// @export
///
/// @examples
/// calculate_mtdna_tumour(
///   mitochondrial_depth = 100,
///   autosome_depth = 10,
///   background_ploidy = 2,
///   tumour_ploidy = 2,
///   tumour_purity = 1
/// )
fn calculate_mtdna_tumour(
    mitochondrial_depth: f32,
    autosome_depth: f32,
    background_ploidy: f32,
    tumour_ploidy: f32,
    tumour_purity: f32,
    per_cell: bool,
) -> f32 {
    let typical_number_of_autosome_copies_per_cell =
        background_ploidy * (1_f32 - tumour_purity) + tumour_ploidy * tumour_purity;
    let mtgenomes_per_autosome =
        typical_number_of_autosome_copies_per_cell * mitochondrial_depth / autosome_depth;

    match per_cell {
        true => mtgenomes_per_autosome * typical_number_of_autosome_copies_per_cell,
        false => mtgenomes_per_autosome,
    }
}

fn hashset_to_string(set: &HashSet<String>) -> String {
    set.iter().cloned().collect::<Vec<_>>().join(", ")
}

fn median_f32(values: &mut [f32]) -> Option<f32> {
    if values.is_empty() {
        return None;
    }

    values.sort_by(|a, b| a.total_cmp(b));

    let mid = values.len() / 2;

    if values.len().is_multiple_of(2) {
        Some((values[mid - 1] + values[mid]) / 2.0)
    } else {
        Some(values[mid])
    }
}
