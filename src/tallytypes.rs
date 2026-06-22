//! Structs that contain different types of tallies

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
