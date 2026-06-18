//! Structs that contain different types of tallies

/// Count of transitions vs transversions
#[derive(Default, Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub struct TallyTiTv {
    pub transition: u64,
    pub transversion: u64,
}

/// Count of transitions vs transversions
#[derive(Default, Debug, Clone, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
pub struct TallySbs6 {
    pub c_a: u64,
    pub c_g: u64,
    pub c_t: u64,
    pub t_a: u64,
    pub t_c: u64,
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

