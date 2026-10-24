# Global reanalysis of public cheese microbiome datasets

Reproducibility repository for the manuscript:

**Global analysis of public cheese microbiomes reveals a conserved
microbial core and reproducible geographic structure**

This repository contains the computational workflow, processed data,
statistical analyses, sensitivity analyses, and figure-generation
scripts used for the global reanalysis of public cheese-associated
microbiome datasets from the NCBI Sequence Read Archive (SRA).

Two complementary reproduction routes are provided:

1.  **Full workflow from NCBI SRA** --- SRA discovery, metadata
    curation, read retrieval, taxonomic profiling, and harmonization.
2.  **Publication-level reproduction** --- reproduction of the
    statistical analyses and manuscript figures from the processed
    matrices included in this repository.

The publication-level route is recommended because complete reprocessing
of thousands of SRA runs requires substantial storage, computation, and
external reference databases.

## Dataset overview

The initial NCBI SRA search identified **14,203 SRA records from 233
BioProjects**. After metadata filtering and curation, the analytical
collection contained **8,910 SRA records from 112 BioProjects**,
including **8,100 amplicon records from 92 BioProjects** and **810 WGS
records from 21 BioProjects**. One BioProject was represented in both
sequencing-strategy categories.

Taxonomic profiling was performed using **QIIME2/DADA2** for amplicon
datasets and **Kraken2** plus **mOTUs v4** for WGS datasets.

## Repository structure

``` text
conf/             Configuration of the full Kedro/SRA workflow
metadata/         Curation, geographic audit, and run-level traceability
processed_data/   Publication-level matrices and analysis source data
scripts/          Harmonization, statistics, figures, and sensitivity analyses
src/              Kedro pipelines from SRA discovery to taxonomic profiling
supplementary/    Supplementary tables
results/          Reproduced statistical and figure outputs
```

## 1. Full workflow from NCBI SRA

### SRA discovery

The SRA discovery workflow is implemented in:

``` text
src/metareanalysis_cheese/pipelines/fetch_runinfo/
conf/base/parameters_fetch_runinfo.yml
```

The search was performed against the **NCBI Sequence Read Archive
(SRA)** using NCBI Entrez Direct. The query required `cheese` and
`metagenome`; the exact executable query and exclusion terms are
retained in `parameters_fetch_runinfo.yml`.

The metadata snapshot used for the study was obtained on **7 February
2026**.

``` bash
kedro run --pipeline=fetch_runinfo
```

### Metadata curation

Metadata filtering is implemented in:

``` text
src/metareanalysis_cheese/pipelines/process_metadata/
conf/base/parameters_process_metadata.yml
```

Only records with `LibrarySource = METAGENOMIC` and sequencing strategy
`AMPLICON` or `WGS` were retained for downstream consideration.

For BioProjects with at least five runs, projects with a whitelist
fraction \>=0.80 were automatically retained, those \<=0.20 were
automatically excluded, and intermediate cases were manually reviewed.
BioProjects with fewer than five runs were also flagged for manual
review.

Manual decisions are retained in:

``` text
metadata/curation/bioprojects_manual_FILLED.tsv
```

A total of **29 BioProjects** underwent manual review: **20 retained and
9 excluded**.

The resulting collection is therefore described as a collection of
**public cheese-associated microbiome datasets**, rather than a
collection restricted exclusively to directly sampled final cheese
matrices.

``` bash
kedro run --pipeline=process_metadata
```

### Sequence retrieval

Sequence retrieval and FASTQ preparation are implemented in:

``` text
src/metareanalysis_cheese/pipelines/sra_fetch_unpack/
conf/base/parameters_sra_fetch_unpack.yml
```

The workflow supports ENA HTTPS retrieval and SRA Toolkit-based
retrieval/unpacking, with paired-end and single-end read detection.

Raw sequencing reads are not redistributed. Temporary FASTQ data are
written to `work/sra_stage_fastq/`, which is excluded from version
control.

``` bash
kedro run --pipeline=sra_fetch_unpack
```

## 2. Amplicon profiling

Amplicon processing is implemented in:

``` text
src/metareanalysis_cheese/pipelines/sra_analyze_amplicon/
conf/base/parameters_sra_analyze_amplicon.yml
```

The workflow used QIIME2 2026.1, DADA2, automatic paired/single-end
handling, Phred+33 encoding, `denoise-paired` or `denoise-single`,
truncation lengths of 0, remaining DADA2 defaults, and the SILVA 138 99%
Naive Bayes classifier.

No separate global primer-removal step or variable-region-specific
processing branch was applied.

Successful profiles were obtained for **7,649/8,100 SRA records from
88/92 BioProjects**.

The SILVA classifier is an external resource expected at:

``` text
external_resources/silva-138-99-nb-classifier.qza
```

``` bash
kedro run --pipeline=sra_analyze_amplicon
```

## 3. WGS profiling

WGS processing is implemented in:

``` text
src/metareanalysis_cheese/pipelines/sra_analyze_wgs/
conf/base/parameters_sra_analyze_wgs.yml
```

### Kraken2

Kraken2 was applied using a locally installed reference database.
Publication matrices were constructed from Kraken2 report percentages at
genus (`G`) and species (`S`) ranks.

Successful Kraken2 profiles were obtained for **634 SRA records from 20
BioProjects**.

The database is expected at:

``` text
external_resources/kraken2_db/
```

### mOTUs

mOTUs v4 was applied independently to the WGS collection. Count profiles
were subsequently converted to relative abundance.

Successful profiles were obtained for **618 SRA records from 19
BioProjects**. These same **618 runs from 19 BioProjects** constitute
the directly paired Kraken2--mOTUs comparison.

``` bash
kedro run --pipeline=sra_analyze_wgs
```

## 4. Taxonomic harmonization

BioProject-level matrix construction is implemented in:

``` text
scripts/02_harmonization/01_build_bioproject_taxonomic_matrices.R
```

The workflow standardizes taxon labels, removes unassigned/unclassified
entries, converts profiles to relative abundance where required, and
constructs genus- and species-level matrices.

The host genera excluded during publication-level harmonization are:

``` text
Homo, Mus, Rattus, Bos, Sus, Gallus, Canis, Felis
```

Processed matrices used for the manuscript are provided in:

``` text
processed_data/bioproject_matrices/
```

These are the recommended starting point for publication-level
reproduction.

## 5. Geographic metadata

Geographic inference uses **BioSample-reported sample origin**, not
sequencing-center location.

Audit scripts:

``` text
scripts/01_metadata_and_geography/01_audit_biosample_metadata.py
scripts/01_metadata_and_geography/02_assign_sample_origin_country.py
scripts/01_metadata_and_geography/03_compare_center_and_sample_origin.py
scripts/01_metadata_and_geography/04_audit_geographic_metadata_coverage.py
```

Outputs are provided in `metadata/geography/`.

A BioProject was assigned a single sample-origin country only when
available BioSample geographic metadata consistently supported one
country. Multi-country projects were retained as multi-country records.

Among records with sample-origin metadata, the comparison with earlier
sequencing-center-derived labels showed **98.94% agreement**. This is an
SRA-record-weighted metadata comparison, not independent geographic
validation.

For all **73 amplicon BioProjects** classified as having a single
BioSample-reported sample-origin country, origin metadata were available
and concordant for 100% of constituent SRA records.

## 6. Publication-level statistical analysis

This is the recommended route for reproducing the manuscript without
downloading and reprocessing the original sequencing reads.

``` bash
Rscript scripts/03_statistical_analysis/02_geographic_and_compositional_analysis.R
```

The workflow performs Hellinger transformation, Bray-Curtis
dissimilarity calculation, geographic PERMANOVA, PERMDISP, balanced
geographic subsampling, compositional covariate analyses, removal of
dominant cheese-associated genera, and continuous genus-proportion
sensitivity analysis.

The primary geographic analysis contains **73 BioProjects from 20
countries**.

Primary PERMANOVA:

``` text
Country: Df = 19; R2 = 0.3656; F = 1.6078; P = 0.001
```

Because dispersion was significant in the primary unbalanced analysis,
balanced subsampling used three BioProjects per country across eight
countries. Across 1,000 iterations, median R2 was **0.4309** with a
**95% empirical interval of 0.3782--0.4996**.

## 7. Figure reproduction

### Figure 1

``` bash
python scripts/04_figures/01_figure1_geographic_distribution.py
```

Uses BioSample-reported sample-origin metadata.

### Figure 2

``` bash
Rscript scripts/04_figures/02_figure2_dataset_composition.R
```

Sequencing-center-derived geographic metadata in descriptive
supplementary outputs are distinct from sample-origin metadata used for
geographic inference.

### Figure 3

``` bash
Rscript scripts/04_figures/03_figure3_genus_prevalence.R
```

Genus prevalence uses a **0.0001 relative-abundance (0.01%)** detection
threshold.

The publication matrices contain 86 QIIME2 rows, 21 Kraken2 rows, and 20
mOTUs rows. These can include platform-level BioProject strata and are
not the unique-BioProject denominators used for species-level WGS
sensitivity analyses.

### Figure 4

``` bash
Rscript scripts/03_statistical_analysis/02_geographic_and_compositional_analysis.R
Rscript scripts/04_figures/04_figure4_geographic_structure.R
```

### Figure 5A --- method concordance

``` bash
Rscript scripts/04_figures/06_figure5a_method_concordance.R
```

Kraken2--mOTUs is a direct WGS comparison. Comparisons involving QIIME2
represent aggregate patterns across independently assembled amplicon and
WGS collections and are not direct paired-method validation.

### Figure 5B --- widespread species

``` bash
Rscript scripts/04_figures/07_figure5b_widespread_species.R
```

Species-level analyses use **20 unique Kraken2 BioProjects** and **19
unique mOTUs BioProjects** after collapsing platform-level duplicates.

A widespread species must independently satisfy in both pipelines:

-   detection at \>=0.01% relative abundance
-   prevalence \>=50%
-   mean abundance \>=0.01%

Twelve species satisfied these criteria; three were present in \>=75% of
BioProjects in both pipelines.

### Figure 5C/D

``` bash
Rscript scripts/04_figures/08_figure5cd_species_analysis.R
```

Kraken2 and mOTUs shared **448 species** before filtering. For panel D,
**103 species** remained after requiring mean prevalence \>=10% and mean
relative abundance \>=0.01%.

Mean species abundances showed:

``` text
Spearman rho = 0.534
P = 6.36 × 10^-9
```

This is an across-species mean-abundance comparison, not a measure of
classification accuracy.

### Assemble Figure 5

``` bash
Rscript scripts/04_figures/09_assemble_figure5.R
```

Final outputs:

``` text
results/figures/figure5/Figure5_FINAL_REVISED.png
results/figures/figure5/Figure5_FINAL_REVISED.pdf
```

## 8. Species-level sensitivity analyses

``` bash
Rscript scripts/05_sensitivity_analysis/01_species_core_threshold_sensitivity.R
```

At the primary 0.01% detection threshold, the numbers of species
retained were:

``` text
Prevalence >=10%: 67
Prevalence >=25%: 33
Prevalence >=50%: 12
Prevalence >=75%: 3
Prevalence >=90%: 0
```

Eleven of the twelve widespread species remained above the 50% criterion
in every leave-one-BioProject-out iteration.

Source data are provided under `processed_data/core_species/`.

## 9. Run-level traceability

Run-level processing information is provided in:

``` text
metadata/processing/Supplementary_Table_S11_run_level_processing.tsv
```

The table links individual SRA accessions to BioProject, BioSample,
sequencing strategy, library layout, platform, metadata-curation
decision and source, BioSample-reported sample origin, geographic
evidence/status, and QIIME2/Kraken2/mOTUs profiling status.

## 10. External resources

Large external resources are intentionally not distributed:

``` text
external_resources/
├── silva-138-99-nb-classifier.qza
└── kraken2_db/
```

Their paths are configured in:

``` text
conf/base/parameters_sra_analyze_amplicon.yml
conf/base/parameters_sra_analyze_wgs.yml
```

## 11. Software environment

The workflow uses separate software environments for orchestration and taxonomic profiling. Exported Mamba environments are provided in `environments/`:

- `kedro_cheese_meta.yml` — Kedro orchestration, metadata processing, SRA retrieval utilities, and downstream analysis dependencies.

- `qiime2-amplicon-2026.1.yml` — QIIME2 2026.1 environment used for amplicon profiling with DADA2 and taxonomic classification.

- `motus4.yml` — mOTUs 4.0.4 environment used for WGS taxonomic profiling.

Kraken2 was executed as a system-installed executable (`Kraken2 2.1.2`) rather than through a dedicated Mamba environment. The Kraken2 database itself is not distributed with this repository because of its size; the expected repository-relative database location is configured in `conf/base/parameters_sra_analyze_wgs.yml`.

The pipeline invokes the QIIME2 and mOTUs environments through `mamba run`, allowing the main Kedro workflow to remain in the `kedro_cheese_meta` environment.

Python project configuration and dependencies are provided in
`pyproject.toml` and `requirements.txt`.

Major software includes Kedro, NCBI Entrez Direct, SRA Toolkit, QIIME2
2026.1, DADA2, SILVA 138 99%, Kraken2, mOTUs v4, and R.

Downstream R analyses use packages including `data.table`, `dplyr`,
`vegan`, `ggplot2`, `ggrepel`, `MASS`, and `magick`.

## 12. Recommended reproduction order

``` bash
# Statistical analyses
Rscript scripts/03_statistical_analysis/02_geographic_and_compositional_analysis.R

# Figures 1–4
python scripts/04_figures/01_figure1_geographic_distribution.py
Rscript scripts/04_figures/02_figure2_dataset_composition.R
Rscript scripts/04_figures/03_figure3_genus_prevalence.R
Rscript scripts/04_figures/04_figure4_geographic_structure.R

# Figure 5
Rscript scripts/04_figures/06_figure5a_method_concordance.R
Rscript scripts/04_figures/07_figure5b_widespread_species.R
Rscript scripts/04_figures/08_figure5cd_species_analysis.R
Rscript scripts/04_figures/09_assemble_figure5.R

# Species sensitivity
Rscript scripts/05_sensitivity_analysis/01_species_core_threshold_sensitivity.R

# Supplementary outputs
Rscript scripts/05_sensitivity_analysis/03_generate_supplementary_outputs.R
```

## 13. Registered Kedro pipelines

``` text
fetch_runinfo
process_metadata
summarize_runinfo
by_continent
sra_fetch_unpack
sra_analyze_wgs
sra_analyze_amplicon
sra_orchestrator
```

Individual stages can be run separately, for example:

``` bash
kedro run --pipeline=fetch_runinfo
kedro run --pipeline=process_metadata
kedro run --pipeline=sra_fetch_unpack
kedro run --pipeline=sra_analyze_amplicon
kedro run --pipeline=sra_analyze_wgs
```

The complete registered workflow can also be invoked with:

``` bash
kedro run
```

A complete raw-data rerun is computationally intensive. For verification
of the manuscript analyses, the processed-data route is recommended.

## 14. Data availability

Original sequencing reads remain publicly available from the NCBI
Sequence Read Archive under the SRA accessions represented in the
run-level metadata.

Raw reads are therefore not duplicated here.

Processed BioProject-level taxonomic matrices, geographic metadata
audits, concordance source data, species-level sensitivity data, and
supplementary tables required for publication-level reproduction are
included in the repository.

## Citation

Citation information will be updated following publication of the
associated manuscript.
