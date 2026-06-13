# bioprep

> [!WARNING]
> This package is in early development and API is not yet stable 

Convert common genomic biological fileformats into basic tabular data structures.

##  Rationale
This project supports the [clinsigs](https:://github.com/selkamand/clinsigs) nextflow pipeline for mutational signature analysis. 
We focus on supporting filetype flavours produced by the oncoanalyser pipeline (e.g. purple enriched VCFs) as this pipeline.


## Installation

Binaries and installation scripts are available in [releases](https://github.com/selkamand/bioprep/releases/)

Or to compile development version from source run:

```
cargo install --git https://github.com/selkamand/bioprep
```

## Structural variant VCFs

Large genomic breaks are commonly represented in VCF files.
The VCF specification supports multiple representations this but bioprep expects:

1. One row per breakend
2. ALT alleles encode explicit breakends with square-bracket syntax (`GCA[1:100[``) as opposed to Symbolic ALT alleles (e.g. <DEL>/<INS>)
3. Because each breakend may not confidently placed (e.g. imagine a break in a homopolymer region), 
confidence intervals (bases below / above ) are should be described as INFO field (CIPOS). 
When CIPOS field is not found we will assume the break is precisely located at POS.
4. Paired breakends are linked by an INFO field (MATEID) that matches the ID field of the paired breakend.
5. For single breakends (where only 1 side of the breakpoint is known) the INFO/MATEID field is absent
<!-- 6. QUAL scores are assigned at the breakpoint level (paired breakends have the same quality). Conversions -->

### BEDPE

> [!NOTE]
> `bioprep svcf bedpe` Outputs ONLY breakpoints where both breakends are described (and have FILTER=PASS). 
> Single breakends are excluded.

```
bioprep svcf -i <svcf> --from purple --to bedpe
```

Outputs a BEDPE-like format data with one row per breakpoint including the following columns: 

1. chrom1: Chromosome of one side of first breakend in pair.
2. start1: Zero-based starting position of the lower confidence interval of first breakend.
3. end1: One-based end position of the upper confidence interval of first breakend.
4. chrom2: Chromosome of second breakend in pair.
5. start2: Zero-based starting position of the lower confidence interval of second breakend in pair.
6. end2: One-based end position of the upper confidence interval of second breakend in pair.
7. name: Breakpoint identifier.
8. score: quality score (from first breakpoint in svcf)
9. strand1: strand of the first breakend in pair
10. strand2: strand for the second breakend in pair

Plus additional columns: (downstream tools like bedtools allow any number of additional columns - these will just be passed-through)

11. vaf1: purity adjusted VAF of first breakend in pair (e.g. from PURPLE_VAF info field if `--from purple`)
12. vaf2: purity adjusted VAF of second breakend in pair (e.g. from PURPLE_VAF info field if `--from purple`)
13. pos1: Zero-based position of first breakend in pair (derived from POS field).
14. pos2: Zero-based position of second breakend in pair (derived from POS field).

> [!Warning]
> A header row is included, in contrast to the official BEDPE specification. We'll leave up to user to drop before plugging into bedtools / other tools

> [!Warning]
> Note gridds assigns a QUAL score per breakpoint (so each side of SV should have the same QUAL score - we just use the first one encountered)

> [!NOTE]
> Which breakend is first / second is determined entirely by the order of the input VCF. Presort by coordinate if you want to guarantee lower coord breakends are chrom1/start1/end1 etc. 


### Breakend TSV

> [!WARNING]
> Breakend TSV conversion not implemented yet

```
bioprep svcf -i <vcf> --from purple --to breakend-tsv
```

Outputs a TSV with one row per breakend (each side paired breakpoints will have their own row). 

Columns include: 

1. chrom: Chromosome of breakend.
2. position: 1-based position of breakend as described by POS column in vcf.
3. vaf: Purity adjusted variant allele frequency supporting breakend (e.g. from PURPLE_VAF info field if `--from purple`).
4. id: id of breakend.
5. mateid: id of mate (set to `.` if single breakend)
6. qual: quality of breakend.

## Performance

We are designing this package for ease of future extension over performance or memory considerations.

**Why?**
1. bioprep is in early development & API is unstable. Any optimisation would be premature and can wait until we have better test coverage.
2. bioprep performs filetype conversions that in precision medicine pipelines are not meaningfully rate limiting or bloating resource usage to a point meaningfully detrimental 

**Examples**

Conversion of SV VCF files to any other output involves parsing all breakend and breakpoint data into a `StructuralVariants` object whose memory footprint depends on number of PASS breakends/breakpoints & length of ID names.
This implementation makes new conversions easy to added as transformations of this struct, however no conversions should need all that information in memory at once (or even anything close).
