//! Classify mutations into common schemes and tally their occurance from bioprep standard formats

use crate::{
    error::{Error, Result},
    io::{
        create_tsv_writer, create_versioned_tsv_writer, read_bedpe_tsv, read_mutations_tsv,
        serialize_object_to_writer,
    },
    pubtypes::{self, BreakpointBedpe, Mutation, Strand},
    seqlibutils,
    tallytypes::{
        BreakpointType, TallyBreakpointType, TallySbs6, TallySmallMutationType, TallyTiTv,
    },
};
use seqlib::{
    base::{ConcreteBase, DnaBase},
    mutation::SmallMutationType,
    sbs::{DnaSingleBaseSubstitution, SingleBaseSubstitution},
};
use std::{
    io::{BufWriter, Write},
    path::Path,
};

/// Tally the number of transitions vs transversions.
/// Only considers Single Base Substititions (DNA)
/// and write result to writer
pub fn tally_titv<W: Write>(snv_tsv: &Path, writer: W) -> Result<()> {
    // Create reader
    let mut rdr = read_mutations_tsv(snv_tsv)?;

    // Initialise counter
    let mut tally = TallyTiTv::default();

    // Loop through mutation file
    for result in rdr.deserialize() {
        // Deserialize into mutation class from this crate
        let mutation_bioprep: Mutation = result.map_err(|source| Error::DeserializeMutation {
            path: snv_tsv.to_owned(),
            source: source.into(),
        })?;

        // Convert to a seqlib SmallMutation type to get
        // convenient mutation classification
        let mutation = seqlibutils::mutation_to_seqlib_mutation(mutation_bioprep)?;

        // Convert to a Single Base substition (DNA or RNA)
        let sbs = match DnaSingleBaseSubstitution::try_from(&mutation) {
            Ok(sbs) => sbs,
            Err(_) => continue,
        };

        match sbs.titv() {
            seqlib::mutation::TiTv::Transition => tally.transition += 1,
            seqlib::mutation::TiTv::Transversion => tally.transversion += 1,
        };
    }

    // Create Versioned writer
    let writer = create_versioned_tsv_writer(writer, "Tally (TiTV)")?;

    // Serialize object to writer
    serialize_object_to_writer(writer, tally, "Tally (TiTv)")?;

    Ok(())
}

/// Tally the small mutation types
pub fn tally_small_mutation_types<W: Write>(snv_tsv: &Path, writer: W) -> Result<()> {
    // Create reader
    let mut rdr = read_mutations_tsv(snv_tsv)?;

    // Initialise counter
    let mut tally = TallySmallMutationType::default();

    // Tally TiTv
    for result in rdr.deserialize() {
        let mutation_bioprep: Mutation = result.map_err(|source| Error::DeserializeMutation {
            path: snv_tsv.to_owned(),
            source: source.into(),
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
    // Create Versioned Writer (pre-writes tool name and version to header)
    let writer = create_versioned_tsv_writer(writer, "Tally (TiTv)")?;

    // Serialize object to writer
    serialize_object_to_writer(writer, tally, "Tally (TiTv)")?;

    Ok(())
}

/// Tally SBS6 mutation types
pub fn tally_sbs6<W: Write>(snv_tsv: &Path, writer: W) -> Result<()> {
    // Create reader
    let mut rdr = read_mutations_tsv(snv_tsv)?;

    // Initialise tally count
    let mut tally = TallySbs6::default();

    // Loop through mutation file
    for result in rdr.deserialize() {
        // Deserialize into mutation class from this crate
        let mutation_bioprep: Mutation = result.map_err(|source| Error::DeserializeMutation {
            path: snv_tsv.to_owned(),
            source: source.into(),
        })?;

        // Convert our internal mutation type to seqlib DNA mutation
        // which has much better methods for subclassifying sequences
        let mutation = seqlibutils::mutation_to_seqlib_mutation(mutation_bioprep)?;

        // Attempt conversion to a DNA  single base substitition.
        // If we can not convert, currently continue.
        // If we want error messages we can always match on the specific error types
        // and log relevant info
        let sbs = match DnaSingleBaseSubstitution::try_from(&mutation) {
            Ok(sbs) => sbs,
            Err(_) => continue,
        };

        // Pyrimidine Center
        let pyrimidine_centered = sbs.pyrimidine_center();

        match (
            pyrimidine_centered.reference(),
            pyrimidine_centered.alternative(),
        ) {
            (DnaBase::C, DnaBase::A) => tally.c_a += 1,
            (DnaBase::C, DnaBase::G) => tally.c_g += 1,
            (DnaBase::C, DnaBase::T) => tally.c_t += 1,
            (DnaBase::T, DnaBase::A) => tally.t_a += 1,
            (DnaBase::T, DnaBase::C) => tally.t_c += 1,
            (DnaBase::T, DnaBase::G) => tally.t_g += 1,
            _ => unreachable!(
                "Pyrimidine centering should not allow non-C/T reference bases. Failed for {sbs} If you see this message please create a new github issue"
            ),
        }
    }

    // Create Versioned Writer (pre-writes tool name and version to header)
    let writer = create_versioned_tsv_writer(writer, "Tally (SBS6)")?;

    // Serialize object to writer
    serialize_object_to_writer(writer, tally, "Tally (SBS6)")?;

    Ok(())
}

/// Classify Single base substitutions into 96 different types
/// based on base change and and trinucleotide context of mutated base.
pub fn tally_sbs96(snv_tsv: &Path, _reference: &Path) -> Result<()> {
    // Opening a reader to the SNV TSV file
    let mut mutation_rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(snv_tsv)
        .map_err(|source| Error::ReadTsv {
            path: snv_tsv.to_owned(),
            source: source.into(),
        })?;

    // Create writer (to stdout)
    let stdout = std::io::stdout().lock();
    let _writer = BufWriter::new(stdout);

    // Iterate through tsv
    for result in mutation_rdr.deserialize() {
        let _mutation: Mutation = result.map_err(|source| Error::DeserializeMutation {
            path: snv_tsv.to_owned(),
            source: source.into(),
        })?;

        // Read mutation
    }
    Ok(())
}

/// Classify breakpoints as Inversions, Deletions, Tandem Duplications, and Inversions
pub fn tally_breakpoint_types<W: Write>(bedpe: &Path, writer: W) -> Result<()> {
    // Create reader
    let mut rdr = read_bedpe_tsv(bedpe)?;

    // Initialise counter
    let mut tally = TallyBreakpointType::default();

    // Iterate through breakpoint bedpe (deserialising into breakpoint bedpe format)
    for result in rdr.deserialize() {
        let breakpoint: BreakpointBedpe =
            result.map_err(|source| Error::DeserializeBreakpoint {
                path: bedpe.to_owned(),
                source: source.into(),
            })?;

        // Infer breakpoint type based on bedpe
        let breakpoint_type = BreakpointType::from_breakpoint_bedpe(&breakpoint);

        match breakpoint_type {
            BreakpointType::Translocation => tally.trans += 1,
            BreakpointType::Deletion => tally.del += 1,
            BreakpointType::Inversion => tally.inv += 1,
            BreakpointType::TandemDuplication => tally.tds += 1,
        };
    }

    // Create Versioned Writer (pre-writes tool name and version to header)
    let writer = create_versioned_tsv_writer(writer, "Tally (BreakpointTypes)")?;

    // Serialize object to writer
    serialize_object_to_writer(writer, tally, "Tally (BreakpointTypes)")?;

    Ok(())
}
