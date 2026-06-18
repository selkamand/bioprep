//! Classify mutations into common schemes and tally their occurance from bioprep standard formats

use crate::{
    error::{Error, Result},
    io::read_mutations_tsv,
    pubtypes::Mutation,
    seqlibutils,
    tallytypes::{TallySbs6, TallySmallMutationType, TallyTiTv},
};
use seqlib::{
    base::{Base, ChemClass, DnaBase},
    context::{ContextWindow, Orientation},
    coord::Pos,
    mutation::{DnaSmallMutation, MutationWithContext, SmallMutationType},
    sequence::{BaseSliceExt, DnaSeq},
};
use std::{
    io::{BufWriter, Write},
    path::Path,
};

/// Tally the number of transitions vs transversions
/// and print result to stdout
pub fn tally_titv<W: Write>(snv_tsv: &Path, writer: W) -> Result<()> {
    // Create reader
    let mut rdr = read_mutations_tsv(snv_tsv)?;

    let mut writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(writer);

    // Initialise counter
    let mut tally = TallyTiTv::default();

    // Tally TiTv
    for result in rdr.deserialize() {
        let mutation_bioprep: Mutation = result.map_err(|source| Error::DeserializeMutation {
            path: snv_tsv.to_owned(),
            source,
        })?;

        let mutation = seqlibutils::mutation_to_seqlib_mutation(mutation_bioprep)?;

        match mutation.titv() {
            Some(val) => match val {
                seqlib::mutation::TiTv::Transition => tally.transition += 1,
                seqlib::mutation::TiTv::Transversion => tally.transition += 1,
            },
            None => continue,
        };
    }

    // Print to stdout

    writer
        .serialize(tally)
        .map_err(|source| Error::write("tally-titv", source))?;

    Ok(())
}

/// Tally SBS6 mutation types
pub fn tally_sbs6<W: Write>(snv_tsv: &Path, writer: W) -> Result<()> {
    // Create reader
    let mut rdr = read_mutations_tsv(snv_tsv)?;

    let mut writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(writer);

    // Initialise counter
    let mut tally = TallySbs6::default();

    // Tally TiTv
    for result in rdr.deserialize() {
        let mutation_bioprep: Mutation = result.map_err(|source| Error::DeserializeMutation {
            path: snv_tsv.to_owned(),
            source,
        })?;

        let mutation = seqlibutils::mutation_to_seqlib_mutation(mutation_bioprep)?;
        //
        // //TODO: replace this with a mutaiton class

        // match mutation.titv() {
        //     Some(val) => match val {
        //         seqlib::mutation::TiTv::Transition => tally.transition += 1,
        //         seqlib::mutation::TiTv::Transversion => tally.transition += 1,
        //     },
        //     None => continue,
        // };
    }

    // Print to stdout

    writer
        .serialize(tally)
        .map_err(|source| Error::write("tally-titv", source))?;

    Ok(())
}

/// Tally the small mutation types
pub fn tally_small_mutation_types<W: Write>(snv_tsv: &Path, writer: W) -> Result<()> {
    // Create reader
    let mut rdr = read_mutations_tsv(snv_tsv)?;

    let mut writer = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(writer);

    // Initialise counter
    let mut tally = TallySmallMutationType::default();

    // Tally TiTv
    for result in rdr.deserialize() {
        let mutation_bioprep: Mutation = result.map_err(|source| Error::DeserializeMutation {
            path: snv_tsv.to_owned(),
            source,
        })?;

        let mutation = seqlibutils::mutation_to_seqlib_mutation(mutation_bioprep)?;

        match mutation.class() {
            SmallMutationType::SNV => tally.snv += 1,
            SmallMutationType::DOUBLET => tally.doublet += 1,
            SmallMutationType::MNV => tally.mnv += 1,
            SmallMutationType::INSERTION => tally.insertion += 1,
            SmallMutationType::DELETION => tally.deletion += 1,
        }
    }

    // Print to stdout
    writer
        .serialize(tally)
        .map_err(|source| Error::write("tally-titv", source))?;

    Ok(())
}

/// Classify Single base substitutions into 96 different types
/// based on base change and and trinucleotide context of mutated base.
pub fn tally_sbs96(snv_tsv: &Path, reference: &Path) -> Result<()> {
    // Opening a reader to the SNV TSV file
    let mut mutation_rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(snv_tsv)
        .map_err(|source| Error::ReadTsv {
            path: snv_tsv.to_owned(),
            source,
        })?;

    // Create writer (to stdout)
    let stdout = std::io::stdout().lock();
    let mut writer = BufWriter::new(stdout);

    // Iterate through tsv
    for result in mutation_rdr.deserialize() {
        let mutation: Mutation = result.map_err(|source| Error::DeserializeMutation {
            path: snv_tsv.to_owned(),
            source,
        })?;

        // Read mutation
    }
    Ok(())
}
