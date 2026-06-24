//! Structs that contain different types of tallies

use crate::pubtypes::{self, BreakpointBedpe, Strand};

/// Count of transitions vs transversions
#[derive(Default, Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub struct TallyTiTv {
    pub transition: u64,
    pub transversion: u64,
}

/// Count of pyrmidine centered single base substitions
#[derive(Default, Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub struct TallySbs6 {
    #[serde(rename = "C>A")]
    pub c_a: u64,
    #[serde(rename = "C>G")]
    pub c_g: u64,
    #[serde(rename = "C>T")]
    pub c_t: u64,
    #[serde(rename = "T>A")]
    pub t_a: u64,
    #[serde(rename = "T>C")]
    pub t_c: u64,
    #[serde(rename = "T>G")]
    pub t_g: u64,
}

/// Count each type of small mutation
#[derive(Default, Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub struct TallySmallMutationType {
    pub snv: u64,
    pub doublet: u64,
    pub mnv: u64,
    pub insertion: u64,
    pub deletion: u64,
}

/// Count of each basic SV type
#[derive(Default, Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub struct TallyBreakpointType {
    #[serde(rename = "Translocation")]
    pub trans: u64,
    #[serde(rename = "Deletion")]
    pub del: u64,
    #[serde(rename = "Inversion")]
    pub inv: u64,
    #[serde(rename = "TandemDuplication")]
    pub tds: u64,
}

pub enum BreakpointType {
    Translocation,
    Deletion,
    Inversion,
    TandemDuplication,
}

impl BreakpointType {
    /// Identify Breakpoint Type of BreakpointBedpe file
    pub fn from_breakpoint_bedpe(breakpoint: &BreakpointBedpe) -> Self {
        if breakpoint.chrom1 != breakpoint.chrom2 {
            return Self::Translocation;
        }

        match (&breakpoint.strand1, &breakpoint.strand2) {
            (Strand::Plus, Strand::Plus) => Self::Inversion,
            (Strand::Plus, Strand::Minus) => Self::Deletion,
            (Strand::Minus, Strand::Plus) => Self::TandemDuplication,
            (Strand::Minus, Strand::Minus) => Self::Inversion,
        }
    }
}

#[derive(Default, Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub struct TallyBreakpointClusterType {
    pub clustered: u64,
    pub unclustered: u64,
}

#[derive(Default, Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub struct TallyBreakpointSize {
    under_10kb: u64,
    under_100kb: u64,
    under_10mb: u64,
    over_10mb: u64,
}

/// SV32 Classification Scheme
#[derive(Default, Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub struct TallySv32 {
    #[serde(rename = "clustered_trans")]
    clustered_trans: u64,
    #[serde(rename = "clustered_del_>10Mb")]
    clustered_del_10mb: u64,
    #[serde(rename = "clustered_inv_>10Mb")]
    clustered_inv_10mb: u64,
    #[serde(rename = "clustered_tds_>10Mb")]
    clustered_tds_10mb: u64,
    #[serde(rename = "non-clustered_trans")]
    non_clustered_trans: u64,
    #[serde(rename = "clustered_del_1-10Kb")]
    clustered_del_1_10kb: u64,
    #[serde(rename = "clustered_inv_1-10Kb")]
    clustered_inv_1_10kb: u64,
    #[serde(rename = "clustered_tds_1-10Kb")]
    clustered_tds_1_10kb: u64,
    #[serde(rename = "clustered_del_10-100Kb")]
    clustered_del_10_100kb: u64,
    #[serde(rename = "clustered_del_1Mb-10Mb")]
    clustered_del_1mb_10mb: u64,
    #[serde(rename = "clustered_inv_10-100Kb")]
    clustered_inv_10_100kb: u64,
    #[serde(rename = "clustered_inv_1Mb-10Mb")]
    clustered_inv_1mb_10mb: u64,
    #[serde(rename = "clustered_tds_10-100Kb")]
    clustered_tds_10_100kb: u64,
    #[serde(rename = "clustered_tds_1Mb-10Mb")]
    clustered_tds_1mb_10mb: u64,
    #[serde(rename = "clustered_del_100Kb-1Mb")]
    clustered_del_100kb_1mb: u64,
    #[serde(rename = "clustered_inv_100Kb-1Mb")]
    clustered_inv_100kb_1mb: u64,
    #[serde(rename = "clustered_tds_100Kb-1Mb")]
    clustered_tds_100kb_1mb: u64,
    #[serde(rename = "non-clustered_del_>10Mb")]
    non_clustered_del_10mb: u64,
    #[serde(rename = "non-clustered_inv_>10Mb")]
    non_clustered_inv_10mb: u64,
    #[serde(rename = "non-clustered_tds_>10Mb")]
    non_clustered_tds_10mb: u64,
    #[serde(rename = "non-clustered_del_1-10Kb")]
    non_clustered_del_1_10kb: u64,
    #[serde(rename = "non-clustered_inv_1-10Kb")]
    non_clustered_inv_1_10kb: u64,
    #[serde(rename = "non-clustered_tds_1-10Kb")]
    non_clustered_tds_1_10kb: u64,
    #[serde(rename = "non-clustered_del_10-100Kb")]
    non_clustered_del_10_100kb: u64,
    #[serde(rename = "non-clustered_del_1Mb-10Mb")]
    non_clustered_del_1mb_10mb: u64,
    #[serde(rename = "non-clustered_inv_10-100Kb")]
    non_clustered_inv_10_100kb: u64,
    #[serde(rename = "non-clustered_inv_1Mb-10Mb")]
    non_clustered_inv_1mb_10mb: u64,
    #[serde(rename = "non-clustered_tds_10-100Kb")]
    non_clustered_tds_10_100kb: u64,
    #[serde(rename = "non-clustered_tds_1Mb-10Mb")]
    non_clustered_tds_1mb_10mb: u64,
    #[serde(rename = "non-clustered_del_100Kb-1Mb")]
    non_clustered_del_100kb_1mb: u64,
    #[serde(rename = "non-clustered_inv_100Kb-1Mb")]
    non_clustered_inv_100kb_1mb: u64,
    #[serde(rename = "non-clustered_tds_100Kb-1Mb")]
    non_clustered_tds_100kb_1mb: u64,
}
