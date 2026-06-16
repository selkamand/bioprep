//! Configure Tool Specific Settings

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SvToolConfig {
    pub vaf_field: String,
}

/// Tools that can produce SVCFs we might want to convert
#[derive(Debug, Clone, PartialEq, Eq, clap::ValueEnum)]
pub enum SvTool {
    Purple,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnvToolConfig {
    pub vaf_field: String,
    _vaf_unadjusted_field: String,
    _depth_field: String,
}

/// Tools that can produce SVCFs we might want to convert
#[derive(Debug, Clone, PartialEq, Eq, clap::ValueEnum)]
pub enum SnvTool {
    Purple,
}

/// Fetch the configuration for a specific tool that produced the svcf
/// being analysed
pub fn configure_for_sv_tool(tool: SvTool) -> SvToolConfig {
    match tool {
        SvTool::Purple => SvToolConfig {
            vaf_field: "PURPLE_AF".to_owned(),
        },
    }
}

/// Fetch the configuration for a specific tool that produced the vcf
/// being analysed
pub fn configure_for_snv_tool(tool: SnvTool) -> SnvToolConfig {
    match tool {
        SnvTool::Purple => SnvToolConfig {
            vaf_field: "PURPLE_AF".to_owned(),
            _vaf_unadjusted_field: "AF".to_owned(),
            _depth_field: "DP".to_owned(),
        },
    }
}
