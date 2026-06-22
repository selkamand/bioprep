use crate::error::Error;
use crate::error::Result;
use seqlib::{mutation::IupacDnaSmallMutation, sequence::IupacDnaSeq};

pub fn mutation_to_seqlib_mutation(
    mutation: crate::pubtypes::Mutation,
) -> Result<IupacDnaSmallMutation> {
    let alt_sequence = IupacDnaSeq::new(&mutation.alternative).map_err(|source| {
        Error::InvalidSequenceForConversion {
            field: "alternative".to_owned(),
            source: source.into(),
        }
    })?;

    let ref_sequence = IupacDnaSeq::new(&mutation.reference).map_err(|source| {
        Error::InvalidSequenceForConversion {
            field: "reference".to_owned(),
            source: source.into(),
        }
    })?;

    let pos = seqlib::coord::Pos::try_from(mutation.pos).map_err(|source| {
        Error::InvalidPositionForConversion {
            source: source.into(),
        }
    })?;

    Ok(IupacDnaSmallMutation::new(
        mutation.chrom,
        pos,
        ref_sequence,
        alt_sequence,
        Some(seqlib::coord::Strand::Positive),
        false,
        true,
    ))
}
