//! Classify mutations into common schemes and tally their occurance from bioprep standard formats

use crate::{
    error::{Error, Result},
    pubtypes::Mutation,
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
    }
    Ok(())
}
