#!/usr/bin/env Rscript

# ============================================================
# Reviewer 3 — master revision analysis
#
# Purpose:
#   Address Reviewer 3 general comments 1–3 concerning:
#   1) sample-origin country vs sequencing-center country
#   2) PERMANOVA / PERMDISP / imbalance / sparse countries
#   3) validation of inferred starter-associated composition types
#
# IMPORTANT:
#   - Does NOT rerun QIIME2.
#   - Does NOT overwrite original Figure 4 outputs.
#   - Uses the same BioProject-level genus matrix as original Fig. 4.
#   - All outputs are written to results/statistics/geographic_analysis/
#
# Original Figure 4 logic preserved where appropriate:
#   BioProject genus matrix
#       -> Hellinger transformation
#       -> Bray-Curtis
#       -> PERMANOVA (999 permutations)
#       -> PERMDISP
#
# Additional revision analyses:
#   - confirmed sample-origin country
#   - countries represented by >=2 BioProjects
#   - balanced subsampling
#   - k = 2:10 cluster validation
#   - cluster stability
#   - country effect after removing 5 starter-associated genera
#   - continuous starter-genus model
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(vegan)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

set.seed(12345)

# ============================================================
# 1. USER SETTINGS
# ============================================================

# Repository root (portable; independent of current working directory)
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)

if (length(file_arg) == 1) {
  script_file <- normalizePath(sub("^--file=", "", file_arg))
  repo_dir <- normalizePath(file.path(dirname(script_file), "..", ".."))
} else {
  repo_dir <- normalizePath(getwd())
}

project_dir <- repo_dir
taxa_file <- file.path(project_dir, "processed_data/bioproject_matrices/QIIME2_Amplicon_genus_BioProject_taxa_wide.csv")

origin_file <- file.path(
  project_dir,
  "metadata/geography/sra_origin_country_audit.tsv"
)

center_origin_file <- file.path(
  project_dir,
  "metadata/geography/sra_center_vs_sample_origin.tsv"
)

out_dir <- file.path(project_dir, "results/statistics/geographic_analysis")

n_perm <- 999

# Number of balanced bootstrap/subsampling iterations.
# 1000 is publication-appropriate and still computationally manageable
# for this BioProject-level dataset.
n_balanced_iter <- 1000

# Repeated k-means runs used to assess stability
n_cluster_stability_iter <- 250

# k values to evaluate
k_range <- 2:10

# Original manuscript value
original_k <- 5

starter_taxa_requested <- c(
  "Lactococcus",
  "Streptococcus",
  "Lacticaseibacillus",
  "Lactobacillus",
  "Leuconostoc"
)

# Minimum independent BioProjects per country for the
# "adequately represented country" sensitivity analysis.
min_projects_country <- 2

# ============================================================
# 2. OUTPUT DIRECTORIES
# ============================================================

dirs <- list(
  root = out_dir,
  origin = file.path(out_dir, "01_origin_audit"),
  confirmed = file.path(out_dir, "02_confirmed_origin"),
  min2 = file.path(out_dir, "03_country_min2"),
  balanced = file.path(out_dir, "04_balanced_subsampling"),
  cluster = file.path(out_dir, "05_cluster_validation"),
  starter_removed = file.path(out_dir, "06_starter_removed"),
  continuous = file.path(out_dir, "07_continuous_starter"),
  fig4 = file.path(out_dir, "08_figure4_source_tables"),
  supplementary = file.path(out_dir, "09_supplementary_tables")
)

invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# ============================================================
# 3. HELPERS
# ============================================================

msg <- function(...) {
  message(sprintf(...))
}

safe_fwrite <- function(x, path) {
  fwrite(as.data.frame(x), path, sep = "\t", na = "NA")
}

extract_adonis_term <- function(model, term) {

  x <- as.data.frame(model)
  x$Term <- rownames(x)

  hit <- x[x$Term == term, , drop = FALSE]

  if (nrow(hit) == 0) {
    return(data.frame(
      Term = term,
      Df = NA_real_,
      SumOfSqs = NA_real_,
      R2 = NA_real_,
      F = NA_real_,
      p = NA_real_
    ))
  }

  data.frame(
    Term = term,
    Df = hit$Df[1],
    SumOfSqs = hit$SumOfSqs[1],
    R2 = hit$R2[1],
    F = hit$F[1],
    p = hit$`Pr(>F)`[1]
  )
}

extract_adonis_table <- function(model) {

  x <- as.data.frame(model)
  x$Term <- rownames(x)
  rownames(x) <- NULL

  x %>%
    select(
      Term,
      everything()
    )
}

extract_permdisp <- function(distance_obj, group, variable_name,
                            permutations = n_perm) {

  bd <- betadisper(distance_obj, group)

  pt <- permutest(
    bd,
    permutations = permutations
  )

  data.frame(
    Variable = variable_name,
    n_groups = nlevels(factor(group)),
    F = unname(pt$tab[1, "F"]),
    p = unname(pt$tab[1, "Pr(>F)"])
  )
}

# ------------------------------------------------------------
# Run the exact Figure-4-style Country / Cluster analysis
# ------------------------------------------------------------

run_country_cluster_models <- function(micro_mat,
                                       meta,
                                       country_col,
                                       cluster_col = "Cluster",
                                       permutations = n_perm) {

  stopifnot(all(rownames(micro_mat) == meta$BioProject))

  micro_hel <- decostand(
    micro_mat,
    method = "hellinger"
  )

  bc <- vegdist(
    micro_hel,
    method = "bray"
  )

  model_meta <- data.frame(
    BioProject = meta$BioProject,
    Country = factor(meta[[country_col]]),
    Cluster = factor(meta[[cluster_col]])
  )

  rownames(model_meta) <- model_meta$BioProject

  ad_country <- adonis2(
    bc ~ Country,
    data = model_meta,
    permutations = permutations,
    by = "terms"
  )

  ad_cluster <- adonis2(
    bc ~ Cluster,
    data = model_meta,
    permutations = permutations,
    by = "terms"
  )

  ad_country_cluster <- adonis2(
    bc ~ Country + Cluster,
    data = model_meta,
    permutations = permutations,
    by = "terms"
  )

  ad_cluster_country <- adonis2(
    bc ~ Cluster + Country,
    data = model_meta,
    permutations = permutations,
    by = "terms"
  )

  country_total <- ad_country$R2[1]
  cluster_total <- ad_cluster$R2[1]

  cc <- as.data.frame(ad_cluster_country)
  cc$Term <- rownames(cc)

  cs <- as.data.frame(ad_country_cluster)
  cs$Term <- rownames(cs)

  country_unique <- cc$R2[cc$Term == "Country"]
  cluster_unique <- cs$R2[cs$Term == "Cluster"]

  shared <- country_total - country_unique

  residual <- 1 - (
    country_unique +
      cluster_unique +
      shared
  )

  partition <- data.frame(
    Component = c(
      "Country unique",
      "Shared country-composition type",
      "Composition type unique",
      "Residual"
    ),
    R2 = c(
      country_unique,
      shared,
      cluster_unique,
      residual
    )
  )

  permanova_summary <- bind_rows(
    extract_adonis_term(ad_country, "Country") %>%
      mutate(Model = "Country only"),

    extract_adonis_term(ad_cluster, "Cluster") %>%
      mutate(Model = "Composition type only"),

    extract_adonis_term(ad_country_cluster, "Country") %>%
      mutate(Model = "Country + Composition type"),

    extract_adonis_term(ad_country_cluster, "Cluster") %>%
      mutate(Model = "Country + Composition type"),

    extract_adonis_term(ad_cluster_country, "Cluster") %>%
      mutate(Model = "Composition type + Country"),

    extract_adonis_term(ad_cluster_country, "Country") %>%
      mutate(Model = "Composition type + Country")
  ) %>%
    select(
      Model,
      Term,
      Df,
      SumOfSqs,
      R2,
      F,
      p
    )

  dispersion <- bind_rows(
    extract_permdisp(
      bc,
      model_meta$Country,
      "Country",
      permutations
    ),
    extract_permdisp(
      bc,
      model_meta$Cluster,
      "Composition type",
      permutations
    )
  )

  list(
    distance = bc,
    meta = model_meta,
    country = ad_country,
    cluster = ad_cluster,
    country_cluster = ad_country_cluster,
    cluster_country = ad_cluster_country,
    permanova_summary = permanova_summary,
    partition = partition,
    permdisp = dispersion
  )
}

# ------------------------------------------------------------
# Country-only model
# ------------------------------------------------------------

run_country_only <- function(micro_mat,
                             meta,
                             country_col,
                             permutations = n_perm) {

  micro_hel <- decostand(
    micro_mat,
    method = "hellinger"
  )

  bc <- vegdist(
    micro_hel,
    method = "bray"
  )

  model_meta <- data.frame(
    BioProject = meta$BioProject,
    Country = factor(meta[[country_col]])
  )

  rownames(model_meta) <- model_meta$BioProject

  ad <- adonis2(
    bc ~ Country,
    data = model_meta,
    permutations = permutations,
    by = "terms"
  )

  disp <- extract_permdisp(
    bc,
    model_meta$Country,
    "Country",
    permutations
  )

  list(
    distance = bc,
    adonis = ad,
    summary = extract_adonis_term(ad, "Country"),
    permdisp = disp
  )
}

# ============================================================
# 4. CHECK INPUT FILES
# ============================================================

required_files <- c(
  taxa_file,
  origin_file,
  center_origin_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing input files:",
      paste(missing_files, collapse = "\n")
    )
  )
}

msg("All required input files found.")

# ============================================================
# 5. LOAD INPUT DATA
# ============================================================

taxa <- fread(taxa_file)
origin <- fread(origin_file)
center_origin <- fread(center_origin_file)

msg("QIIME2 BioProject matrix: %d rows x %d columns",
    nrow(taxa), ncol(taxa))

msg("Origin audit: %d SRA", nrow(origin))
msg("Center-vs-origin audit: %d SRA", nrow(center_origin))

required_taxa_meta <- c(
  "BioProject",
  "CenterName_standardized",
  "Continent",
  "Country",
  "Platform"
)

missing_meta <- setdiff(
  required_taxa_meta,
  colnames(taxa)
)

if (length(missing_meta) > 0) {
  stop(
    paste(
      "Missing required metadata columns in taxa matrix:",
      paste(missing_meta, collapse = ", ")
    )
  )
}

# ============================================================
# 6. IDENTIFY TAXON COLUMNS
# ============================================================

metadata_cols <- c(
  "BioProject",
  "CenterName_standardized",
  "Continent",
  "Country",
  "Platform"
)

tax_cols <- setdiff(
  colnames(taxa),
  metadata_cols
)

tax_cols <- tax_cols[
  sapply(taxa[, ..tax_cols], is.numeric)
]

starter_taxa <- intersect(
  starter_taxa_requested,
  tax_cols
)

if (length(starter_taxa) != length(starter_taxa_requested)) {

  warning(
    paste(
      "Not all requested starter-associated genera were found.",
      "Detected:",
      paste(starter_taxa, collapse = ", ")
    )
  )
}

if (length(starter_taxa) < 2) {
  stop("Too few starter-associated genera detected.")
}

msg(
  "Starter-associated genera detected: %s",
  paste(starter_taxa, collapse = ", ")
)

# ============================================================
# 7. REBUILD BIOPROJECT TAXON MATRIX EXACTLY AS FIGURE 4
# ============================================================

taxa_bp <- taxa %>%
  group_by(BioProject) %>%
  summarise(
    Old_CenterCountry = {
      x <- na.omit(Country)
      if (length(x) == 0) NA_character_ else as.character(x[1])
    },
    Continent_old = {
      x <- na.omit(Continent)
      if (length(x) == 0) NA_character_ else as.character(x[1])
    },
    Platform = paste(
      sort(unique(na.omit(Platform))),
      collapse = "; "
    ),
    across(
      all_of(tax_cols),
      \(x) mean(x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

msg("BioProjects in original Figure-4 matrix: %d", nrow(taxa_bp))

# ============================================================
# 8. BUILD SAMPLE-ORIGIN COUNTRY AT BIOPROJECT LEVEL
# ============================================================

origin_amp <- origin %>%
  filter(SequencingStrategy == "AMPLICON")

# Run-level coverage
origin_run_summary <- origin_amp %>%
  mutate(
    Is_confirmed = OriginStatus == "confirmed" &
      !is.na(SampleOriginCountry) &
      SampleOriginCountry != ""
  ) %>%
  group_by(BioProject) %>%
  summarise(
    n_SRA_total = n(),
    n_SRA_confirmed_origin = sum(Is_confirmed),
    pct_SRA_confirmed_origin =
      100 * n_SRA_confirmed_origin / n_SRA_total,

    n_confirmed_countries = n_distinct(
      SampleOriginCountry[Is_confirmed]
    ),

    confirmed_countries = paste(
      sort(unique(
        SampleOriginCountry[Is_confirmed]
      )),
      collapse = ";"
    ),

    .groups = "drop"
  )

# Confirmed country assignments only
origin_confirmed_bp <- origin_amp %>%
  filter(
    OriginStatus == "confirmed",
    !is.na(SampleOriginCountry),
    SampleOriginCountry != ""
  ) %>%
  group_by(BioProject) %>%
  summarise(
    n_confirmed_origin_SRA = n(),
    n_confirmed_origin_BioSamples = n_distinct(BioSample),
    n_origin_countries = n_distinct(SampleOriginCountry),
    SampleOriginCountry = ifelse(
      n_origin_countries == 1,
      first(SampleOriginCountry),
      NA_character_
    ),
    all_confirmed_origin_countries = paste(
      sort(unique(SampleOriginCountry)),
      collapse = ";"
    ),
    .groups = "drop"
  )

origin_bp <- taxa_bp %>%
  select(
    BioProject,
    Old_CenterCountry
  ) %>%
  left_join(
    origin_run_summary,
    by = "BioProject"
  ) %>%
  left_join(
    origin_confirmed_bp,
    by = "BioProject"
  ) %>%
  mutate(
    Origin_assignment = case_when(
      is.na(n_confirmed_origin_SRA) |
        n_confirmed_origin_SRA == 0 ~
        "no_confirmed_origin",

      n_origin_countries == 1 ~
        "single_confirmed_origin",

      n_origin_countries > 1 ~
        "multiple_confirmed_origins",

      TRUE ~
        "unknown"
    ),

    Center_vs_origin = case_when(
      Origin_assignment != "single_confirmed_origin" ~
        "not_comparable",

      is.na(Old_CenterCountry) ~
        "not_comparable",

      Old_CenterCountry == SampleOriginCountry ~
        "agreement",

      TRUE ~
        "disagreement"
    )
  )

safe_fwrite(
  origin_bp,
  file.path(
    dirs$origin,
    "S01_Amplicon_BioProject_origin_audit.tsv"
  )
)

# ============================================================
# 9. RUN-LEVEL CENTER VS ORIGIN AUDIT FOR AMPLICON
# ============================================================

center_amp <- center_origin %>%
  filter(SequencingStrategy == "AMPLICON")

center_origin_summary <- center_amp %>%
  count(
    Comparison,
    name = "n_SRA"
  ) %>%
  mutate(
    pct_all_Amplicon_SRA =
      100 * n_SRA / sum(n_SRA)
  )

safe_fwrite(
  center_origin_summary,
  file.path(
    dirs$origin,
    "S02_Amplicon_center_vs_origin_summary.tsv"
  )
)

center_origin_disagreements <- center_amp %>%
  filter(Comparison == "disagreement") %>%
  count(
    CenterCountry,
    SampleOriginCountry,
    BioProject,
    name = "n_SRA"
  ) %>%
  arrange(
    desc(n_SRA),
    BioProject
  )

safe_fwrite(
  center_origin_disagreements,
  file.path(
    dirs$origin,
    "S03_Amplicon_center_vs_origin_disagreements.tsv"
  )
)

# ============================================================
# 10. CONFIRMED-ORIGIN ANALYTIC DATASET
# ============================================================

confirmed_bp <- taxa_bp %>%
  inner_join(
    origin_bp %>%
      filter(
        Origin_assignment == "single_confirmed_origin"
      ) %>%
      select(
        BioProject,
        SampleOriginCountry,
        n_SRA_total,
        n_SRA_confirmed_origin,
        pct_SRA_confirmed_origin,
        Center_vs_origin
      ),
    by = "BioProject"
  )

# Remove zero-community BioProjects
confirmed_mat <- as.matrix(
  confirmed_bp[, tax_cols]
)

rownames(confirmed_mat) <- confirmed_bp$BioProject
confirmed_mat[is.na(confirmed_mat)] <- 0

keep_nonzero <- rowSums(confirmed_mat) > 0

confirmed_bp <- confirmed_bp[
  keep_nonzero,
]

confirmed_mat <- confirmed_mat[
  keep_nonzero,
  ,
  drop = FALSE
]

stopifnot(
  all(
    rownames(confirmed_mat) ==
      confirmed_bp$BioProject
  )
)

msg(
  "Confirmed-origin Figure-4 BioProjects: %d",
  nrow(confirmed_bp)
)

msg(
  "Confirmed-origin countries: %d",
  n_distinct(confirmed_bp$SampleOriginCountry)
)

# ============================================================
# 11. COUNTRY COVERAGE TABLE
# ============================================================

country_coverage <- confirmed_bp %>%
  group_by(SampleOriginCountry) %>%
  summarise(
    n_BioProjects = n_distinct(BioProject),
    n_SRA_total = sum(n_SRA_total, na.rm = TRUE),
    n_SRA_confirmed_origin =
      sum(n_SRA_confirmed_origin, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    desc(n_BioProjects),
    desc(n_SRA_total),
    SampleOriginCountry
  )

safe_fwrite(
  country_coverage,
  file.path(
    dirs$origin,
    "S04_confirmed_origin_country_coverage.tsv"
  )
)

# ============================================================
# 11B. PRIMARY GEOGRAPHIC MODEL — ALL CONFIRMED-ORIGIN PROJECTS
#
# Geography is evaluated on all BioProjects with a single
# confirmed BioSample-reported sample-origin country.
# This analysis is independent of the presence of the five
# genera used later to construct composition types.
# ============================================================

model_country_primary <- run_country_only(
  micro_mat = confirmed_mat,
  meta = confirmed_bp,
  country_col = "SampleOriginCountry"
)

safe_fwrite(
  model_country_primary$summary,
  file.path(
    dirs$confirmed,
    "PERMANOVA_Country_only_ALL_confirmed_origin.tsv"
  )
)

safe_fwrite(
  model_country_primary$permdisp,
  file.path(
    dirs$confirmed,
    "PERMDISP_Country_ALL_confirmed_origin.tsv"
  )
)

msg(
  paste0(
    "PRIMARY geography: %d BioProjects; %d countries; ",
    "Country R2=%.4f, F=%.4f, p=%.4g; ",
    "PERMDISP F=%.4f, p=%.4g"
  ),
  nrow(confirmed_bp),
  n_distinct(confirmed_bp$SampleOriginCountry),
  model_country_primary$summary$R2,
  model_country_primary$summary$F,
  model_country_primary$summary$p,
  model_country_primary$permdisp$F,
  model_country_primary$permdisp$p
)

# ============================================================
# 12. REBUILD ORIGINAL k=5 COMPOSITION TYPES
#
# Important:
# clustering is recomputed on the confirmed-origin dataset
# so all inferential analyses use exactly the same projects.
# ============================================================

starter_mat <- confirmed_mat[
  ,
  starter_taxa,
  drop = FALSE
]

starter_positive <- rowSums(starter_mat) > 0

if (any(!starter_positive)) {

  excluded_no_starter <- data.frame(
    BioProject = rownames(starter_mat)[!starter_positive],
    Reason = "zero abundance across all five starter-associated genera"
  )

  safe_fwrite(
    excluded_no_starter,
    file.path(
      dirs$cluster,
      "BioProjects_excluded_zero_starter_abundance.tsv"
    )
  )
}

confirmed_bp_cluster <- confirmed_bp[
  starter_positive,
]

confirmed_mat_cluster <- confirmed_mat[
  starter_positive,
  ,
  drop = FALSE
]

starter_mat <- starter_mat[
  starter_positive,
  ,
  drop = FALSE
]

starter_rel <- sweep(
  starter_mat,
  1,
  rowSums(starter_mat),
  "/"
)

starter_rel[is.na(starter_rel)] <- 0

set.seed(123)
km5 <- kmeans(
  starter_rel,
  centers = original_k,
  nstart = 100
)

confirmed_bp_cluster$Cluster_raw <- paste0(
  "C",
  km5$cluster
)

# Relabel exactly as in original Figure 4
cluster_profile <- as.data.frame(starter_rel)
cluster_profile$Cluster_raw <-
  confirmed_bp_cluster$Cluster_raw

cluster_profile_mean <- cluster_profile %>%
  group_by(Cluster_raw) %>%
  summarise(
    across(
      all_of(starter_taxa),
      \(x) mean(x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    Dominant_starter = starter_taxa[
      which.max(
        c_across(all_of(starter_taxa))
      )
    ]
  ) %>%
  ungroup()

desired_order <- c(
  "Lacticaseibacillus",
  "Leuconostoc",
  "Lactococcus",
  "Streptococcus",
  "Lactobacillus"
)

cluster_profile_mean <- cluster_profile_mean %>%
  mutate(
    Dominant_starter = factor(
      Dominant_starter,
      levels = desired_order
    )
  ) %>%
  arrange(
    Dominant_starter,
    Cluster_raw
  ) %>%
  mutate(
    Cluster = paste0(
      "C",
      row_number()
    )
  )

cluster_map <- cluster_profile_mean %>%
  select(
    Cluster_raw,
    Cluster,
    Dominant_starter
  )

confirmed_bp_cluster <- confirmed_bp_cluster %>%
  left_join(
    cluster_map,
    by = "Cluster_raw"
  ) %>%
  mutate(
    Cluster = factor(
      Cluster,
      levels = paste0(
        "C",
        1:original_k
      )
    )
  )

# Add continuous starter proportions
starter_rel_df <- as.data.frame(starter_rel)
starter_rel_df$BioProject <- rownames(starter_rel)

confirmed_bp_cluster <- confirmed_bp_cluster %>%
  left_join(
    starter_rel_df,
    by = "BioProject",
    suffix = c("", "_starter_rel")
  )

# Restore matrix row order after join
confirmed_bp_cluster <- confirmed_bp_cluster %>%
  arrange(
    match(
      BioProject,
      rownames(confirmed_mat_cluster)
    )
  )

stopifnot(
  all(
    confirmed_bp_cluster$BioProject ==
      rownames(confirmed_mat_cluster)
  )
)

# ============================================================
# 13. SAVE k=5 SOURCE TABLES
# ============================================================

cluster_membership <- confirmed_bp_cluster %>%
  select(
    BioProject,
    Old_CenterCountry,
    SampleOriginCountry,
    Center_vs_origin,
    Cluster,
    Cluster_raw,
    all_of(starter_taxa)
  )

safe_fwrite(
  cluster_membership,
  file.path(
    dirs$fig4,
    "Figure4_confirmed_origin_cluster_membership.tsv"
  )
)

safe_fwrite(
  cluster_map,
  file.path(
    dirs$fig4,
    "Figure4_cluster_relabeling_map.tsv"
  )
)

cluster_comp <- confirmed_bp_cluster %>%
  group_by(Cluster) %>%
  summarise(
    n_BioProjects = n(),
    across(
      all_of(starter_taxa),
      \(x) mean(x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

safe_fwrite(
  cluster_comp,
  file.path(
    dirs$fig4,
    "Figure4_cluster_composition.tsv"
  )
)

country_cluster_counts <- confirmed_bp_cluster %>%
  count(
    SampleOriginCountry,
    Cluster,
    name = "n_BioProjects"
  )

safe_fwrite(
  country_cluster_counts,
  file.path(
    dirs$fig4,
    "Figure4_country_by_composition_type_counts.tsv"
  )
)

# ============================================================
# 14. CONFIRMED-ORIGIN FULL FIGURE-4 MODEL
# ============================================================

model_confirmed <- run_country_cluster_models(
  micro_mat = confirmed_mat_cluster,
  meta = confirmed_bp_cluster,
  country_col = "SampleOriginCountry",
  cluster_col = "Cluster"
)

safe_fwrite(
  model_confirmed$permanova_summary,
  file.path(
    dirs$confirmed,
    "PERMANOVA_confirmed_origin_country_composition_type.tsv"
  )
)

safe_fwrite(
  model_confirmed$partition,
  file.path(
    dirs$confirmed,
    "Variance_partition_confirmed_origin.tsv"
  )
)

safe_fwrite(
  model_confirmed$permdisp,
  file.path(
    dirs$confirmed,
    "PERMDISP_confirmed_origin.tsv"
  )
)

capture.output(
  model_confirmed$country,
  file = file.path(
    dirs$confirmed,
    "PERMANOVA_Country_only.txt"
  )
)

capture.output(
  model_confirmed$cluster,
  file = file.path(
    dirs$confirmed,
    "PERMANOVA_CompositionType_only.txt"
  )
)

capture.output(
  model_confirmed$country_cluster,
  file = file.path(
    dirs$confirmed,
    "PERMANOVA_Country_then_CompositionType.txt"
  )
)

capture.output(
  model_confirmed$cluster_country,
  file = file.path(
    dirs$confirmed,
    "PERMANOVA_CompositionType_then_Country.txt"
  )
)

# ============================================================
# 15. COUNTRIES WITH >=2 INDEPENDENT BIOPROJECTS
#
# IMPORTANT:
# The primary geographic sensitivity analysis is constructed
# independently of composition-type availability.
#
# Two datasets are therefore retained:
#   A) geographic subset: all confirmed-origin BioProjects
#      from countries represented by >=2 BioProjects
#   B) composition-type subset: starter-positive BioProjects
#      from the same eligible countries
# ============================================================

adequate_countries <- country_coverage %>%
  filter(
    n_BioProjects >= min_projects_country
  ) %>%
  pull(SampleOriginCountry)

# ------------------------------------------------------------
# A. Geographic sensitivity dataset
# ------------------------------------------------------------

min2_geo_meta <- confirmed_bp %>%
  filter(
    SampleOriginCountry %in% adequate_countries
  )

min2_geo_mat <- confirmed_mat[
  min2_geo_meta$BioProject,
  ,
  drop = FALSE
]

msg(
  ">= %d BioProjects/country GEOGRAPHIC subset: %d BioProjects, %d countries",
  min_projects_country,
  nrow(min2_geo_meta),
  n_distinct(min2_geo_meta$SampleOriginCountry)
)

safe_fwrite(
  min2_geo_meta %>%
    count(
      SampleOriginCountry,
      name = "n_BioProjects"
    ) %>%
    arrange(
      desc(n_BioProjects),
      SampleOriginCountry
    ),
  file.path(
    dirs$min2,
    "country_group_sizes_geographic.tsv"
  )
)

model_min2_geo <- run_country_only(
  micro_mat = min2_geo_mat,
  meta = min2_geo_meta,
  country_col = "SampleOriginCountry"
)

safe_fwrite(
  model_min2_geo$summary,
  file.path(
    dirs$min2,
    "PERMANOVA_min2_Country_only_ALL_eligible.tsv"
  )
)

safe_fwrite(
  model_min2_geo$permdisp,
  file.path(
    dirs$min2,
    "PERMDISP_min2_Country_ALL_eligible.tsv"
  )
)

# ------------------------------------------------------------
# B. Composition-type subset
#
# Only analyses explicitly involving composition type require
# nonzero abundance across the five composition-associated genera.
# ------------------------------------------------------------

min2_meta <- confirmed_bp_cluster %>%
  filter(
    SampleOriginCountry %in% adequate_countries
  )

min2_mat <- confirmed_mat_cluster[
  min2_meta$BioProject,
  ,
  drop = FALSE
]

msg(
  ">= %d BioProjects/country COMPOSITION-TYPE subset: %d BioProjects, %d countries",
  min_projects_country,
  nrow(min2_meta),
  n_distinct(min2_meta$SampleOriginCountry)
)

safe_fwrite(
  min2_meta %>%
    count(
      SampleOriginCountry,
      name = "n_BioProjects"
    ) %>%
    arrange(
      desc(n_BioProjects),
      SampleOriginCountry
    ),
  file.path(
    dirs$min2,
    "country_group_sizes_composition_type.tsv"
  )
)

model_min2 <- run_country_cluster_models(
  micro_mat = min2_mat,
  meta = min2_meta,
  country_col = "SampleOriginCountry",
  cluster_col = "Cluster"
)

safe_fwrite(
  model_min2$permanova_summary,
  file.path(
    dirs$min2,
    "PERMANOVA_min2_country_composition_type.tsv"
  )
)

safe_fwrite(
  model_min2$partition,
  file.path(
    dirs$min2,
    "Variance_partition_min2.tsv"
  )
)

safe_fwrite(
  model_min2$permdisp,
  file.path(
    dirs$min2,
    "PERMDISP_min2.tsv"
  )
)

# ============================================================
# 16. BALANCED COUNTRY SUBSAMPLING
#
# We balance BioProject counts, not SRA counts.
#
# Only countries represented by >=3 independent BioProjects
# are retained. Three BioProjects are sampled without
# replacement from each retained country in every iteration.
#
# This avoids the special two-observation-per-group geometry
# that makes PERMDISP difficult to interpret for n=2.
#
# The analysis is a sensitivity analysis for unequal country
# representation and is not interpreted as an independent
# estimate of variance explained by geography.
# ============================================================

balanced_n <- 3

country_sizes_balanced <- confirmed_bp %>%
  count(
    SampleOriginCountry,
    name = "n_BioProjects"
  )

balanced_countries <- country_sizes_balanced %>%
  filter(
    n_BioProjects >= balanced_n
  ) %>%
  pull(SampleOriginCountry)

balanced_meta <- confirmed_bp %>%
  filter(
    SampleOriginCountry %in% balanced_countries
  )

balanced_mat <- confirmed_mat[
  balanced_meta$BioProject,
  ,
  drop = FALSE
]

msg(
  "Balanced subsampling: %d countries x %d BioProjects/country = %d BioProjects/iteration",
  length(balanced_countries),
  balanced_n,
  length(balanced_countries) * balanced_n
)

balanced_results <- vector(
  "list",
  n_balanced_iter
)

set.seed(20261021)

for (i in seq_len(n_balanced_iter)) {

  selected <- balanced_meta %>%
    group_by(SampleOriginCountry) %>%
    slice_sample(
      n = balanced_n,
      replace = FALSE
    ) %>%
    ungroup()

  selected <- selected %>%
    arrange(
      match(
        BioProject,
        rownames(balanced_mat)
      )
    )

  mat_i <- balanced_mat[
    selected$BioProject,
    ,
    drop = FALSE
  ]

  fit_i <- tryCatch(
    run_country_only(
      mat_i,
      selected,
      "SampleOriginCountry",
      permutations = n_perm
    ),
    error = function(e) NULL
  )

  if (is.null(fit_i)) {

    balanced_results[[i]] <- data.frame(
      iteration = i,
      n_countries = length(balanced_countries),
      n_per_country = balanced_n,
      n_BioProjects = nrow(selected),
      Df = NA_real_,
      R2 = NA_real_,
      F = NA_real_,
      p = NA_real_,
      PERMDISP_F = NA_real_,
      PERMDISP_p = NA_real_
    )

  } else {

    balanced_results[[i]] <- data.frame(
      iteration = i,
      n_countries = length(balanced_countries),
      n_per_country = balanced_n,
      n_BioProjects = nrow(selected),
      Df = fit_i$summary$Df,
      R2 = fit_i$summary$R2,
      F = fit_i$summary$F,
      p = fit_i$summary$p,
      PERMDISP_F = fit_i$permdisp$F,
      PERMDISP_p = fit_i$permdisp$p
    )
  }

  if (i %% 50 == 0) {
    msg(
      "Balanced subsampling: %d / %d",
      i,
      n_balanced_iter
    )
  }
}

balanced_results <- bind_rows(
  balanced_results
)

safe_fwrite(
  balanced_results,
  file.path(
    dirs$balanced,
    "balanced_subsampling_all_iterations.tsv"
  )
)

balanced_summary <- data.frame(
  Metric = c(
    "Number of iterations",
    "Countries",
    "BioProjects per country",
    "BioProjects per iteration",
    "Median R2",
    "R2 2.5 percentile",
    "R2 97.5 percentile",
    "Mean R2",
    "SD R2",
    "Median F",
    "Fraction PERMANOVA p<=0.05",
    "Fraction PERMANOVA p<=0.01",
    "Fraction PERMANOVA p<=0.001",
    "Median PERMDISP F",
    "Fraction PERMDISP p<=0.05"
  ),
  Value = c(
    nrow(balanced_results),
    length(balanced_countries),
    balanced_n,
    balanced_n * length(balanced_countries),

    median(
      balanced_results$R2,
      na.rm = TRUE
    ),

    quantile(
      balanced_results$R2,
      0.025,
      na.rm = TRUE
    ),

    quantile(
      balanced_results$R2,
      0.975,
      na.rm = TRUE
    ),

    mean(
      balanced_results$R2,
      na.rm = TRUE
    ),

    sd(
      balanced_results$R2,
      na.rm = TRUE
    ),

    median(
      balanced_results$F,
      na.rm = TRUE
    ),

    mean(
      balanced_results$p <= 0.05,
      na.rm = TRUE
    ),

    mean(
      balanced_results$p <= 0.01,
      na.rm = TRUE
    ),

    mean(
      balanced_results$p <= 0.001,
      na.rm = TRUE
    ),

    median(
      balanced_results$PERMDISP_F,
      na.rm = TRUE
    ),

    mean(
      balanced_results$PERMDISP_p <= 0.05,
      na.rm = TRUE
    )
  )
)

safe_fwrite(
  balanced_summary,
  file.path(
    dirs$balanced,
    "balanced_subsampling_summary.tsv"
  )
)

# ============================================================
# 17. CLUSTER VALIDATION: k = 2:10
#
# Metrics:
#   - total within-cluster SS
#   - between-cluster SS
#   - between / total SS
#   - average silhouette width
#   - minimum cluster size
#   - maximum cluster size
#
# This assesses whether k=5 is empirically supported.
# ============================================================

starter_dist <- dist(
  starter_rel,
  method = "euclidean"
)

k_validation <- vector(
  "list",
  length(k_range)
)

names(k_validation) <- k_range

for (k in k_range) {

  set.seed(1000 + k)

  km <- kmeans(
    starter_rel,
    centers = k,
    nstart = 100
  )

sil <- cluster::silhouette(
  km$cluster,
  starter_dist
)

  cluster_sizes <- table(
    km$cluster
  )

  k_validation[[as.character(k)]] <- data.frame(
    k = k,
    tot_withinss = km$tot.withinss,
    betweenss = km$betweenss,
    totss = km$totss,
    between_to_total =
      km$betweenss / km$totss,
    average_silhouette =
      mean(sil[, "sil_width"]),
    min_cluster_size =
      min(cluster_sizes),
    max_cluster_size =
      max(cluster_sizes)
  )
}

k_validation <- bind_rows(
  k_validation
)

safe_fwrite(
  k_validation,
  file.path(
    dirs$cluster,
    "k2_to_k10_validation_metrics.tsv"
  )
)

# ============================================================
# 18. K-MEANS INITIALIZATION REPRODUCIBILITY
#
# Compare repeated k-means solutions fitted to the SAME
# starter_rel matrix with a reference solution using the
# Adjusted Rand Index (ARI).
#
# IMPORTANT:
# This evaluates reproducibility across random k-means
# initializations. It is NOT BioProject-resampling stability
# and must not be interpreted as biological cluster stability.
#
# Implement ARI directly to avoid an extra package dependency.
# ============================================================

adjusted_rand_index <- function(x, y) {

  tab <- table(x, y)

  choose2 <- function(z) {
    z * (z - 1) / 2
  }

  sum_comb <- sum(
    choose2(tab)
  )

  sum_row <- sum(
    choose2(rowSums(tab))
  )

  sum_col <- sum(
    choose2(colSums(tab))
  )

  n <- sum(tab)
  total_comb <- choose2(n)

  if (total_comb == 0) {
    return(NA_real_)
  }

  expected <- (
    sum_row * sum_col
  ) / total_comb

  max_index <- 0.5 * (
    sum_row + sum_col
  )

  denom <- max_index - expected

  if (denom == 0) {
    return(NA_real_)
  }

  (
    sum_comb - expected
  ) / denom
}

cluster_stability_all <- list()

for (k in k_range) {

  set.seed(5000 + k)

  reference <- kmeans(
    starter_rel,
    centers = k,
    nstart = 250
  )$cluster

  tmp <- vector(
    "list",
    n_cluster_stability_iter
  )

  for (i in seq_len(n_cluster_stability_iter)) {

    set.seed(
      100000 +
        k * 1000 +
        i
    )

    km_i <- kmeans(
      starter_rel,
      centers = k,
      nstart = 25
    )

    ari <- adjusted_rand_index(
      reference,
      km_i$cluster
    )

    tmp[[i]] <- data.frame(
      k = k,
      iteration = i,
      ARI_vs_reference = ari
    )
  }

  cluster_stability_all[[as.character(k)]] <-
    bind_rows(tmp)

  msg(
    "Cluster stability k=%d completed",
    k
  )
}

cluster_stability_all <- bind_rows(
  cluster_stability_all
)

safe_fwrite(
  cluster_stability_all,
  file.path(
    dirs$cluster,
    "cluster_stability_all_iterations.tsv"
  )
)

cluster_stability_summary <- cluster_stability_all %>%
  group_by(k) %>%
  summarise(
    mean_ARI = mean(
      ARI_vs_reference,
      na.rm = TRUE
    ),
    median_ARI = median(
      ARI_vs_reference,
      na.rm = TRUE
    ),
    q025_ARI = quantile(
      ARI_vs_reference,
      0.025,
      na.rm = TRUE
    ),
    q975_ARI = quantile(
      ARI_vs_reference,
      0.975,
      na.rm = TRUE
    ),
    min_ARI = min(
      ARI_vs_reference,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

safe_fwrite(
  cluster_stability_summary,
  file.path(
    dirs$cluster,
    "cluster_stability_summary.tsv"
  )
)

cluster_validation_combined <- k_validation %>%
  left_join(
    cluster_stability_summary,
    by = "k"
  )

safe_fwrite(
  cluster_validation_combined,
  file.path(
    dirs$cluster,
    "cluster_validation_combined.tsv"
  )
)

# ============================================================
# 19. SENSITIVITY: REMOVE THE FIVE STARTER-ASSOCIATED GENERA
#
# Tests whether country-associated structure persists after
# removing the genera used to construct composition types.
# ============================================================

nonstarter_taxa <- setdiff(
  colnames(confirmed_mat),
  starter_taxa
)

micro_no_starter <- confirmed_mat[
  ,
  nonstarter_taxa,
  drop = FALSE
]

keep_no_starter <- rowSums(
  micro_no_starter
) > 0

micro_no_starter <- micro_no_starter[
  keep_no_starter,
  ,
  drop = FALSE
]

meta_no_starter <- confirmed_bp[
  keep_no_starter,
]

model_no_starter <- run_country_only(
  micro_no_starter,
  meta_no_starter,
  "SampleOriginCountry"
)

safe_fwrite(
  model_no_starter$summary,
  file.path(
    dirs$starter_removed,
    "PERMANOVA_country_after_removing_5_starter_genera.tsv"
  )
)

safe_fwrite(
  model_no_starter$permdisp,
  file.path(
    dirs$starter_removed,
    "PERMDISP_country_after_removing_5_starter_genera.tsv"
  )
)

removed_taxa_table <- data.frame(
  Removed_genus = starter_taxa
)

safe_fwrite(
  removed_taxa_table,
  file.path(
    dirs$starter_removed,
    "removed_starter_associated_genera.tsv"
  )
)

# Repeat on adequately represented countries
meta_no_starter_min2 <- meta_no_starter %>%
  filter(
    SampleOriginCountry %in% adequate_countries
  )

mat_no_starter_min2 <- micro_no_starter[
  meta_no_starter_min2$BioProject,
  ,
  drop = FALSE
]

model_no_starter_min2 <- run_country_only(
  mat_no_starter_min2,
  meta_no_starter_min2,
  "SampleOriginCountry"
)

safe_fwrite(
  model_no_starter_min2$summary,
  file.path(
    dirs$starter_removed,
    "PERMANOVA_country_after_removing_5_starter_genera_min2.tsv"
  )
)

safe_fwrite(
  model_no_starter_min2$permdisp,
  file.path(
    dirs$starter_removed,
    "PERMDISP_country_after_removing_5_starter_genera_min2.tsv"
  )
)

# ============================================================
# 20. CONTINUOUS STARTER-GENUS PROPORTIONS
#
# Reviewer requested examining the geographic result using
# continuous starter-genus proportions and/or community
# composition excluding those genera.
#
# The five proportions are normalized WITHIN these five genera
# for each BioProject and therefore sum to 1. Including all five
# creates exact linear dependence. We therefore use four explicit
# within-five-genus proportions and treat the fifth as the
# reference component.
#
# To avoid choosing the reference arbitrarily, select the genus
# with the lowest variance as the omitted reference and report it.
#
# Models:
#   starter continuous only
#   starter continuous + Country
#
# Country term is tested after continuous composition covariates.
# ============================================================

starter_rel_for_model <- starter_rel[
  confirmed_bp_cluster$BioProject,
  ,
  drop = FALSE
]

starter_variances <- apply(
  starter_rel_for_model,
  2,
  var,
  na.rm = TRUE
)

reference_starter <- names(
  which.min(starter_variances)
)

continuous_starters <- setdiff(
  starter_taxa,
  reference_starter
)

continuous_meta <- confirmed_bp_cluster %>%
  select(
    BioProject,
    SampleOriginCountry
  )

for (g in continuous_starters) {
  continuous_meta[[g]] <-
    starter_rel_for_model[
      continuous_meta$BioProject,
      g
    ]
}

continuous_meta$Country <- factor(
  continuous_meta$SampleOriginCountry
)

micro_hel_cont <- decostand(
  confirmed_mat_cluster,
  method = "hellinger"
)

bc_cont <- vegdist(
  micro_hel_cont,
  method = "bray"
)

starter_formula_string <- paste(
  continuous_starters,
  collapse = " + "
)

formula_starter_only <- as.formula(
  paste(
    "bc_cont ~",
    starter_formula_string
  )
)

formula_starter_country <- as.formula(
  paste(
    "bc_cont ~",
    starter_formula_string,
    "+ Country"
  )
)

set.seed(12345)

adonis_cont_starter <- adonis2(
  formula_starter_only,
  data = continuous_meta,
  permutations = n_perm,
  by = "terms"
)

set.seed(12345)

adonis_cont_country <- adonis2(
  formula_starter_country,
  data = continuous_meta,
  permutations = n_perm,
  by = "terms"
)

continuous_model_info <- data.frame(
  Item = c(
    "Starter-associated genera",
    "Reference genus omitted because proportions sum to 1",
    "Continuous genera entered into model",
    "Transformation of community matrix",
    "Distance",
    "Permutations",
    "Country definition"
  ),
  Value = c(
    paste(starter_taxa, collapse = "; "),
    reference_starter,
    paste(continuous_starters, collapse = "; "),
    "Hellinger",
    "Bray-Curtis",
    n_perm,
    "BioSample-reported sample-origin country"
  )
)

safe_fwrite(
  continuous_model_info,
  file.path(
    dirs$continuous,
    "continuous_starter_model_definition.tsv"
  )
)

safe_fwrite(
  extract_adonis_table(
    adonis_cont_starter
  ),
  file.path(
    dirs$continuous,
    "PERMANOVA_continuous_starter_only.tsv"
  )
)

safe_fwrite(
  extract_adonis_table(
    adonis_cont_country
  ),
  file.path(
    dirs$continuous,
    "PERMANOVA_continuous_starter_then_country.tsv"
  )
)

country_after_continuous <- extract_adonis_term(
  adonis_cont_country,
  "Country"
)

safe_fwrite(
  country_after_continuous,
  file.path(
    dirs$continuous,
    "Country_effect_after_continuous_starter_proportions.tsv"
  )
)

# ============================================================
# 21. CONTINUOUS MODEL ON >=2 BIOPROJECT COUNTRIES
# ============================================================

continuous_meta_min2 <- continuous_meta %>%
  filter(
    SampleOriginCountry %in% adequate_countries
  )

bc_cont_min2 <- vegdist(
  decostand(
    confirmed_mat_cluster[
      continuous_meta_min2$BioProject,
      ,
      drop = FALSE
    ],
    method = "hellinger"
  ),
  method = "bray"
)

formula_cont_min2 <- as.formula(
  paste(
    "bc_cont_min2 ~",
    starter_formula_string,
    "+ Country"
  )
)

set.seed(12345)

adonis_cont_country_min2 <- adonis2(
  formula_cont_min2,
  data = continuous_meta_min2,
  permutations = n_perm,
  by = "terms"
)

safe_fwrite(
  extract_adonis_table(
    adonis_cont_country_min2
  ),
  file.path(
    dirs$continuous,
    "PERMANOVA_continuous_starter_then_country_min2.tsv"
  )
)

country_after_continuous_min2 <-
  extract_adonis_term(
    adonis_cont_country_min2,
    "Country"
  )

safe_fwrite(
  country_after_continuous_min2,
  file.path(
    dirs$continuous,
    "Country_effect_after_continuous_starter_proportions_min2.tsv"
  )
)

# ============================================================
# 22. EXACT ANALYSIS METADATA / METHODS TABLE
# ============================================================

analysis_metadata <- data.frame(
  Parameter = c(
    "Input profiling method",
    "Taxonomic rank",
    "Independent analysis unit",
    "Original matrix",
    "Country variable",
    "Community transformation",
    "Distance metric",
    "PERMANOVA function",
    "PERMANOVA permutations",
    "PERMANOVA R2 type",
    "PERMDISP function",
    "Original composition-type genera",
    "Original clustering method",
    "Original k",
    "Original k-means nstart",
    "Adequate-country threshold",
    "Balanced-subsampling unit",
    "Balanced countries criterion",
    "Balanced BioProjects per country",
    "Balanced iterations",
    "Cluster validation k range",
    "K-means initialization reproducibility iterations"
  ),
  Value = c(
    "QIIME2 amplicon",
    "Genus",
    "BioProject",
    taxa_file,
    "BioSample-reported sample-origin country",
    "Hellinger",
    "Bray-Curtis",
    "vegan::adonis2",
    n_perm,
    "Unadjusted R2 reported by adonis2",
    "vegan::betadisper + permutest",
    paste(starter_taxa, collapse = "; "),
    "k-means on within-BioProject normalized proportions of five genera",
    original_k,
    100,
    paste0(">=", min_projects_country, " independent BioProjects per country"),
    "BioProject",
    "Countries with >=3 independent BioProjects",
    balanced_n,
    n_balanced_iter,
    paste(range(k_range), collapse = "-"),
    n_cluster_stability_iter
  )
)

safe_fwrite(
  analysis_metadata,
  file.path(
    dirs$supplementary,
    "S_Methods_exact_analysis_specification.tsv"
  )
)

# ============================================================
# 23. ANALYSIS FLOW TABLE
# ============================================================

analysis_flow <- data.frame(
  Stage = c(
    "Original Figure-4 QIIME2 matrix",
    "Single confirmed BioSample-reported sample-origin country",
    "Primary geographic analysis",
    paste0(
      "Geographic sensitivity: countries with >=",
      min_projects_country,
      " BioProjects"
    ),
    "Nonzero five-genus composition for composition-type analyses",
    "Balanced geographic sensitivity per iteration"
  ),
  n_BioProjects = c(
    nrow(taxa_bp),
    nrow(confirmed_bp),
    nrow(confirmed_bp),
    nrow(min2_geo_meta),
    nrow(confirmed_bp_cluster),
    balanced_n * length(balanced_countries)
  ),
  n_Countries = c(
    n_distinct(
      taxa_bp$Old_CenterCountry,
      na.rm = TRUE
    ),
    n_distinct(
      confirmed_bp$SampleOriginCountry
    ),
    n_distinct(
      confirmed_bp$SampleOriginCountry
    ),
    n_distinct(
      min2_geo_meta$SampleOriginCountry
    ),
    n_distinct(
      confirmed_bp_cluster$SampleOriginCountry
    ),
    length(balanced_countries)
  )
)

safe_fwrite(
  analysis_flow,
  file.path(
    dirs$supplementary,
    "S_analysis_flow.tsv"
  )
)

# ============================================================
# 24. ORIGINAL VS CONFIRMED COUNTRY TABLE
# ============================================================

old_vs_new <- confirmed_bp_cluster %>%
  select(
    BioProject,
    Old_CenterCountry,
    SampleOriginCountry,
    Center_vs_origin,
    n_SRA_total,
    n_SRA_confirmed_origin,
    pct_SRA_confirmed_origin
  ) %>%
  arrange(
    Center_vs_origin,
    BioProject
  )

safe_fwrite(
  old_vs_new,
  file.path(
    dirs$supplementary,
    "S_old_center_country_vs_confirmed_sample_origin.tsv"
  )
)

# ============================================================
# 25. COMBINED KEY STATISTICS TABLE
# ============================================================

get_country_row <- function(x) {

  z <- x %>%
    filter(Term == "Country")

  if (nrow(z) == 0) {
    return(
      data.frame(
        Df = NA,
        R2 = NA,
        F = NA,
        p = NA
      )
    )
  }

  z[1, c(
    "Df",
    "R2",
    "F",
    "p"
  )]
}

confirmed_country_only <- model_confirmed$permanova_summary %>%
  filter(
    Model == "Country only",
    Term == "Country"
  )

confirmed_country_after_cluster <-
  model_confirmed$permanova_summary %>%
  filter(
    Model == "Composition type + Country",
    Term == "Country"
  )

min2_country_only <- model_min2$permanova_summary %>%
  filter(
    Model == "Country only",
    Term == "Country"
  )

min2_country_after_cluster <-
  model_min2$permanova_summary %>%
  filter(
    Model == "Composition type + Country",
    Term == "Country"
  )

key_statistics <- bind_rows(

  data.frame(
    Analysis =
      "Confirmed origin: Country only",
    n_BioProjects =
      nrow(confirmed_bp_cluster),
    n_Countries =
      n_distinct(
        confirmed_bp_cluster$SampleOriginCountry
      ),
    Df =
      confirmed_country_only$Df,
    R2 =
      confirmed_country_only$R2,
    F =
      confirmed_country_only$F,
    p =
      confirmed_country_only$p
  ),

  data.frame(
    Analysis =
      "Confirmed origin: Country after composition type",
    n_BioProjects =
      nrow(confirmed_bp_cluster),
    n_Countries =
      n_distinct(
        confirmed_bp_cluster$SampleOriginCountry
      ),
    Df =
      confirmed_country_after_cluster$Df,
    R2 =
      confirmed_country_after_cluster$R2,
    F =
      confirmed_country_after_cluster$F,
    p =
      confirmed_country_after_cluster$p
  ),

  data.frame(
    Analysis =
      "Confirmed origin >=2 projects/country: Country only",
    n_BioProjects =
      nrow(min2_meta),
    n_Countries =
      n_distinct(
        min2_meta$SampleOriginCountry
      ),
    Df =
      min2_country_only$Df,
    R2 =
      min2_country_only$R2,
    F =
      min2_country_only$F,
    p =
      min2_country_only$p
  ),

  data.frame(
    Analysis =
      "Confirmed origin >=2 projects/country: Country after composition type",
    n_BioProjects =
      nrow(min2_meta),
    n_Countries =
      n_distinct(
        min2_meta$SampleOriginCountry
      ),
    Df =
      min2_country_after_cluster$Df,
    R2 =
      min2_country_after_cluster$R2,
    F =
      min2_country_after_cluster$F,
    p =
      min2_country_after_cluster$p
  ),

  data.frame(
    Analysis =
      "Country after removing five starter-associated genera",
    n_BioProjects =
      nrow(meta_no_starter),
    n_Countries =
      n_distinct(
        meta_no_starter$SampleOriginCountry
      ),
    Df =
      model_no_starter$summary$Df,
    R2 =
      model_no_starter$summary$R2,
    F =
      model_no_starter$summary$F,
    p =
      model_no_starter$summary$p
  ),

  data.frame(
    Analysis =
      "Country after continuous starter proportions",
    n_BioProjects =
      nrow(continuous_meta),
    n_Countries =
      n_distinct(
        continuous_meta$SampleOriginCountry
      ),
    Df =
      country_after_continuous$Df,
    R2 =
      country_after_continuous$R2,
    F =
      country_after_continuous$F,
    p =
      country_after_continuous$p
  )
)

safe_fwrite(
  key_statistics,
  file.path(
    dirs$root,
    "reviewer3_key_statistics.tsv"
  )
)

# ============================================================
# 26. FIGURE 4 SOURCE TABLES
# ============================================================

# Heatmap source
fig4_heatmap_source <- confirmed_bp_cluster %>%
  select(
    BioProject,
    SampleOriginCountry,
    Cluster,
    all_of(starter_taxa)
  )

safe_fwrite(
  fig4_heatmap_source,
  file.path(
    dirs$fig4,
    "Figure4A_heatmap_source_confirmed_origin.tsv"
  )
)

# Panel B
fig4B_source <- cluster_comp

safe_fwrite(
  fig4B_source,
  file.path(
    dirs$fig4,
    "Figure4B_composition_type_source.tsv"
  )
)

# Panel C
safe_fwrite(
  model_confirmed$partition,
  file.path(
    dirs$fig4,
    "Figure4C_variance_partition_source.tsv"
  )
)

# Panel D
fig4D_source <- bind_rows(
  model_confirmed$permanova_summary %>%
    filter(
      Model == "Country only",
      Term == "Country"
    ) %>%
    mutate(
      Effect = "Country"
    ),

  model_confirmed$permanova_summary %>%
    filter(
      Model == "Composition type only",
      Term == "Cluster"
    ) %>%
    mutate(
      Effect = "Composition type"
    ),

  model_confirmed$permanova_summary %>%
    filter(
      Model == "Composition type + Country",
      Term == "Country"
    ) %>%
    mutate(
      Effect = "Country | after composition type"
    ),

  model_confirmed$permanova_summary %>%
    filter(
      Model == "Country + Composition type",
      Term == "Cluster"
    ) %>%
    mutate(
      Effect = "Composition type | after country"
    )
) %>%
  select(
    Effect,
    Df,
    R2,
    F,
    p
  )

safe_fwrite(
  fig4D_source,
  file.path(
    dirs$fig4,
    "Figure4D_PERMANOVA_source.tsv"
  )
)

# Panel E
safe_fwrite(
  model_confirmed$permdisp,
  file.path(
    dirs$fig4,
    "Figure4E_PERMDISP_source.tsv"
  )
)

# ============================================================
# 27. SUPPLEMENTARY MASTER INDEX
# ============================================================

supplementary_index <- data.frame(
  Suggested_table = c(
    "Table S_origin",
    "Table S_country_coverage",
    "Table S_Figure4_PERMANOVA",
    "Table S_Figure4_PERMDISP",
    "Table S_balanced",
    "Table S_cluster_validation",
    "Table S_starter_removed",
    "Table S_continuous_starter",
    "Table S_methods"
  ),
  File = c(
    file.path(
      "01_origin_audit",
      "S01_Amplicon_BioProject_origin_audit.tsv"
    ),
    file.path(
      "01_origin_audit",
      "S04_confirmed_origin_country_coverage.tsv"
    ),
    file.path(
      "02_confirmed_origin",
      "PERMANOVA_confirmed_origin_country_composition_type.tsv"
    ),
    file.path(
      "02_confirmed_origin",
      "PERMDISP_confirmed_origin.tsv"
    ),
    file.path(
      "04_balanced_subsampling",
      "balanced_subsampling_summary.tsv"
    ),
    file.path(
      "05_cluster_validation",
      "cluster_validation_combined.tsv"
    ),
    file.path(
      "06_starter_removed",
      "PERMANOVA_country_after_removing_5_starter_genera.tsv"
    ),
    file.path(
      "07_continuous_starter",
      "Country_effect_after_continuous_starter_proportions.tsv"
    ),
    file.path(
      "09_supplementary_tables",
      "S_Methods_exact_analysis_specification.tsv"
    )
  ),
  Purpose = c(
    "Audit of sequencing-center country versus confirmed sample-origin country",
    "Country representation after restricting to confirmed sample origins",
    "Main revised Figure 4 PERMANOVA models",
    "Main revised Figure 4 dispersion diagnostics",
    "Sensitivity to unequal numbers of BioProjects among countries",
    "Silhouette evaluation and reproducibility across repeated k-means initializations",
    "Geographic association after removing the five genera used to define composition types",
    "Geographic association after accounting for continuous starter-associated genus proportions",
    "Exact statistical specification requested by Reviewer 3"
  )
)

safe_fwrite(
  supplementary_index,
  file.path(
    dirs$supplementary,
    "SUPPLEMENTARY_TABLE_INDEX.tsv"
  )
)

# ============================================================
# 28. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    dirs$root,
    "R_sessionInfo.txt"
  )
)

# ============================================================
# 29. HUMAN-READABLE KEY RESULTS
# ============================================================

fmt <- function(x, digits = 4) {
  ifelse(
    is.na(x),
    "NA",
    format(
      round(x, digits),
      nsmall = digits,
      scientific = FALSE
    )
  )
}

get_perm_line <- function(df, model, term) {

  x <- df %>%
    filter(
      Model == model,
      Term == term
    )

  if (nrow(x) == 0) {
    return("not available")
  }

  paste0(
    "Df=",
    x$Df[1],
    ", R2=",
    fmt(x$R2[1]),
    ", F=",
    fmt(x$F[1]),
    ", p=",
    signif(x$p[1], 4)
  )
}

best_sil_k <- k_validation$k[
  which.max(
    k_validation$average_silhouette
  )
]

best_stability_k <-
  cluster_stability_summary$k[
    which.max(
      cluster_stability_summary$median_ARI
    )
  ]

k5_validation <- cluster_validation_combined %>%
  filter(k == original_k)

primary_disp_country <-
  model_country_primary$permdisp %>%
  filter(Variable == "Country")

confirmed_disp_country <-
  model_confirmed$permdisp %>%
  filter(Variable == "Country")

min2_geo_disp_country <-
  model_min2_geo$permdisp %>%
  filter(Variable == "Country")

min2_disp_country <-
  model_min2$permdisp %>%
  filter(Variable == "Country")

summary_lines <- c(

  "REVIEWER 3 — MASTER REVISION ANALYSIS",
  "=====================================",
  "",

  "INPUT",
  paste0(
    "Original QIIME2 BioProject matrix: ",
    nrow(taxa_bp),
    " BioProjects"
  ),
  "Taxonomic rank: genus",
  paste0(
    "Composition-associated genera: ",
    paste(starter_taxa, collapse = ", ")
  ),
  "",

  "GEOGRAPHIC AUDIT",
  paste0(
    "BioProjects with a single BioSample-reported sample-origin country: ",
    nrow(confirmed_bp)
  ),
  paste0(
    "Countries represented: ",
    n_distinct(confirmed_bp$SampleOriginCountry)
  ),
  paste0(
    "BioProjects entering k=5 composition-type analysis: ",
    nrow(confirmed_bp_cluster)
  ),
  "",

  "PRIMARY GEOGRAPHIC ANALYSIS",
  paste0(
    "BioProjects: ",
    nrow(confirmed_bp)
  ),
  paste0(
    "Countries: ",
    n_distinct(confirmed_bp$SampleOriginCountry)
  ),
  paste0(
    "Country only: Df=",
    model_country_primary$summary$Df[1],
    ", R2=",
    fmt(model_country_primary$summary$R2[1]),
    ", F=",
    fmt(model_country_primary$summary$F[1]),
    ", p=",
    signif(model_country_primary$summary$p[1], 4)
  ),
  paste0(
    "Country PERMDISP: F=",
    fmt(primary_disp_country$F[1]),
    ", p=",
    signif(primary_disp_country$p[1], 4)
  ),
  "",

  paste0(
    "REPLICATED-COUNTRY SENSITIVITY (>=",
    min_projects_country,
    " BIOPROJECTS/COUNTRY)"
  ),
  paste0(
    "BioProjects: ",
    nrow(min2_geo_meta)
  ),
  paste0(
    "Countries: ",
    n_distinct(min2_geo_meta$SampleOriginCountry)
  ),
  paste0(
    "Country only: Df=",
    model_min2_geo$summary$Df[1],
    ", R2=",
    fmt(model_min2_geo$summary$R2[1]),
    ", F=",
    fmt(model_min2_geo$summary$F[1]),
    ", p=",
    signif(model_min2_geo$summary$p[1], 4)
  ),
  paste0(
    "Country PERMDISP: F=",
    fmt(min2_geo_disp_country$F[1]),
    ", p=",
    signif(min2_geo_disp_country$p[1], 4)
  ),
  "",

  "COMPOSITION-TYPE ANALYSIS",
  paste0(
    "BioProjects: ",
    nrow(confirmed_bp_cluster)
  ),
  paste0(
    "Composition type only: ",
    get_perm_line(
      model_confirmed$permanova_summary,
      "Composition type only",
      "Cluster"
    )
  ),
  paste0(
    "Country after composition type: ",
    get_perm_line(
      model_confirmed$permanova_summary,
      "Composition type + Country",
      "Country"
    )
  ),
  paste0(
    "Composition type after country: ",
    get_perm_line(
      model_confirmed$permanova_summary,
      "Country + Composition type",
      "Cluster"
    )
  ),
  paste0(
    "Country PERMDISP in composition-type subset: F=",
    fmt(confirmed_disp_country$F[1]),
    ", p=",
    signif(confirmed_disp_country$p[1], 4)
  ),
  "",

  "BALANCED SUBSAMPLING",
  paste0(
    "Iterations: ",
    n_balanced_iter
  ),
  paste0(
    "Countries: ",
    length(balanced_countries)
  ),
  paste0(
    "BioProjects/country: ",
    balanced_n
  ),
  paste0(
    "Median R2: ",
    fmt(
      median(
        balanced_results$R2,
        na.rm = TRUE
      )
    )
  ),
  paste0(
    "95% empirical interval R2: ",
    fmt(
      quantile(
        balanced_results$R2,
        0.025,
        na.rm = TRUE
      )
    ),
    " - ",
    fmt(
      quantile(
        balanced_results$R2,
        0.975,
        na.rm = TRUE
      )
    )
  ),
  paste0(
    "Fraction PERMANOVA p<=0.05: ",
    fmt(
      mean(
        balanced_results$p <= 0.05,
        na.rm = TRUE
      )
    )
  ),
  paste0(
    "Fraction PERMDISP p<=0.05: ",
    fmt(
      mean(
        balanced_results$PERMDISP_p <= 0.05,
        na.rm = TRUE
      )
    )
  ),
  "",

  "CLUSTER VALIDATION",
  paste0(
    "Best k by mean silhouette: ",
    best_sil_k
  ),
  paste0(
    "Best k by median ARI stability: ",
    best_stability_k
  ),
  paste0(
    "k=5 mean silhouette: ",
    fmt(
      k5_validation$average_silhouette
    )
  ),
  paste0(
    "k=5 median ARI: ",
    fmt(
      k5_validation$median_ARI
    )
  ),
  paste0(
    "k=5 95% ARI interval: ",
    fmt(
      k5_validation$q025_ARI
    ),
    " - ",
    fmt(
      k5_validation$q975_ARI
    )
  ),
  "",

  "STARTER-GENUS REMOVAL SENSITIVITY",
  paste0(
    "Removed genera: ",
    paste(
      starter_taxa,
      collapse = ", "
    )
  ),
  paste0(
    "Country after removal: Df=",
    model_no_starter$summary$Df,
    ", R2=",
    fmt(
      model_no_starter$summary$R2
    ),
    ", F=",
    fmt(
      model_no_starter$summary$F
    ),
    ", p=",
    signif(
      model_no_starter$summary$p,
      4
    )
  ),
  "",

  "CONTINUOUS STARTER-PROPORTION SENSITIVITY",
  paste0(
    "Reference proportion omitted: ",
    reference_starter
  ),
  paste0(
    "Country after continuous starter proportions: Df=",
    country_after_continuous$Df,
    ", R2=",
    fmt(
      country_after_continuous$R2
    ),
    ", F=",
    fmt(
      country_after_continuous$F
    ),
    ", p=",
    signif(
      country_after_continuous$p,
      4
    )
  ),
  "",

  "METHODS",
  "Community transformation: Hellinger",
  "Distance: Bray-Curtis",
  paste0(
    "PERMANOVA permutations: ",
    n_perm
  ),
  "PERMANOVA R2: unadjusted R2 returned by vegan::adonis2",
  "Dispersion: vegan::betadisper followed by permutest",
  "",

  paste0(
    "All outputs: ",
    out_dir
  )
)

writeLines(
  summary_lines,
  file.path(
    dirs$root,
    "REVISION_KEY_RESULTS.txt"
  )
)

# ============================================================
# 30. DONE
# ============================================================

message("")
message(
  paste(
    summary_lines,
    collapse = "\n"
  )
)
message("")
message("====================================================")
message("REVIEWER 3 MASTER ANALYSIS COMPLETED")
message("====================================================")
message("Output directory:")
message(out_dir)
message("")
message("Most important file:")
message(
  file.path(
    out_dir,
    "REVISION_KEY_RESULTS.txt"
  )
)
message("====================================================")