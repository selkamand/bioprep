use crate::error::{Error, Result};
use std::{fs::File, path::Path};

pub fn read_mutations_tsv(snv_tsv: &Path) -> Result<csv::Reader<File>> {
    csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(snv_tsv)
        .map_err(|source| Error::ReadTsv {
            path: snv_tsv.to_owned(),
            source,
        })
}
