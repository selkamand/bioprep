use crate::error::{Error, Result};
use std::{fs::File, io::Write, path::Path};

pub fn read_mutations_tsv(snv_tsv: &Path) -> Result<csv::Reader<File>> {
    csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(snv_tsv)
        .map_err(|source| Error::ReadTsv {
            path: snv_tsv.to_owned(),
            source: source.into(),
        })
}

pub fn create_tsv_writer<W: Write>(writer: W) -> csv::Writer<W> {
    csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_writer(writer)
}

/// Serializes an object as a CSV record using an existing CSV writer.
///
/// This is a small helper for writing any [`serde::Serialize`] value through a
/// configured [`csv::Writer`]. The writer is passed in by value so callers can
/// decide where the CSV data goes, such as a file, buffer, or standard output,
/// and can configure writer options before calling this function.
///
/// The `filetype` argument is used only for error reporting. It should describe
/// the kind of file being written, such as `"CSV"`, `"mutation table"`, or
/// `"metadata file"`.
///
/// # Errors
///
/// Returns an error if CSV serialization fails. The underlying [`csv::Error`]
/// is wrapped using [`Error::write`], with `filetype` included in the resulting
/// error message to identify what kind of file could not be written.
pub(crate) fn serialize_object_to_writer<W: Write, T: serde::Serialize>(
    mut writer: csv::Writer<W>,
    object: T,
    filetype: &str,
) -> Result<()> {
    writer
        .serialize(object)
        .map_err(|source| Error::write(filetype, source.into()))?;

    Ok(())
}
