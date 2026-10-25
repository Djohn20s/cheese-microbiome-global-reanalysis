#!/usr/bin/env Rscript

# ============================================================
# Figure 3. Consensus cheese microbiota across analytical pipelines
# Version: v4_all_and_plotted_concordance
#
# Input:
# - QIIME2_Amplicon_genus_BioProject_prevalence_summary.csv
# - Kraken2_WGS_genus_BioProject_prevalence_summary.csv
# - mOTUs_WGS_genus_BioProject_prevalence_summary.csv
#
# Output:
# - grouped prevalence heatmap
# - long source table for plotted genera
# - prevalence matrix for plotted genera
# - mean abundance matrix for plotted genera
# - separate matrices for biological groups
# - Supplementary table with ALL harmonised genera across pipelines
# - Supplementary prevalence and abundance matrices for ALL genera
# - pairwise concordance statistics for:
#     1) all harmonised genera
#     2) selected top 30 genera
#     3) plotted biological-group genera only
# - compact concordance summary for Results section
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(grid)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

# =====================================================
# PATHS
# =====================================================

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
audit_dir <- file.path(
  project_dir,
  "processed_data/bioproject_matrices"
)

out_dir <- file.path(
  project_dir,
  "results/figures/figure3"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =====================================================
# PARAMETERS
# =====================================================

top_n_genera <- 30

# Detection threshold for Jaccard presence/absence.
# Values are already in percent after standardisation.
detection_prevalence_threshold <- 0

# Filter used only for filtered concordance statistics.
# A genus is retained when:
# (prevalence >= 10% in at least one pipeline) AND
# (mean abundance >= 0.01% in at least one pipeline).
concordance_min_prevalence_percent <- 10
concordance_min_abundance_percent <- 0.01

pipeline_levels <- c(
  "QIIME2",
  "Kraken2",
  "mOTUs"
)

excluded_taxa <- c(
  # host taxa
  "Homo",
  "Mus",
  "Rattus",
  "Bos",
  "Sus",
  "Gallus",
  "Canis",
  "Felis",

  # unknown / uninformative labels
  "Unknown",
  "Unclassified",
  "unclassified",

  # non-bacterial / artefacts
  "Mitochondria",
  "Skunavirus",

  # taxonomic artefacts / not directly comparable
  "Escherichia-Shigella",

  # manually excluded pending taxonomic harmonisation
  "Macrococcus",
  "Mammaliicoccus"
)

# =====================================================
# BIOLOGICAL GROUPS
# =====================================================

core_cheese <- c(
  "Lactococcus",
  "Streptococcus",
  "Staphylococcus",
  "Lacticaseibacillus",
  "Leuconostoc",
  "Lactobacillus",
  "Enterococcus",
  "Lactiplantibacillus",
  "Latilactobacillus",
  "Lentilactobacillus",
  "Limosilactobacillus"
)

rind_microbiota <- c(
  "Brevibacterium",
  "Corynebacterium",
  "Brachybacterium",
  "Microbacterium",
  "Psychrobacter",
  "Marinilactibacillus",
  "Chryseobacterium"
)

environmental_taxa <- c(
  "Pseudomonas",
  "Acinetobacter",
  "Vibrio",
  "Enterobacter",
  "Klebsiella",
  "Stenotrophomonas",
  "Bacillus",
  "Clostridium",
  "Streptomyces"
)

group_order <- c(
  "Cheese-associated genera",
  "Rind-associated microbiota",
  "Environmental microbiota"
)

# =====================================================
# FUNCTIONS
# =====================================================

normalize_genus_name <- function(x) {

  x <- as.character(x)
  x <- stringr::str_squish(x)

  # GTDB/mOTUs-style suffixes:
  # Pseudomonas E -> Pseudomonas
  # Lactococcus A -> Lactococcus
  # Enterococcus B -> Enterococcus
  # Serratia J -> Serratia
  x <- stringr::str_replace(
    x,
    "^([A-Za-z-]+)\\s+[A-Z]$",
    "\\1"
  )

  x
}

assign_biological_group <- function(x) {

  dplyr::case_when(
    x %in% core_cheese ~ "Cheese-associated genera",
    x %in% rind_microbiota ~ "Rind-associated microbiota",
    x %in% environmental_taxa ~ "Environmental microbiota",
    TRUE ~ "Other"
  )
}

standardize_prevalence_table <- function(dt, pipeline_name) {

  required_cols <- c(
    "Tax",
    "prevalence_bioproject",
    "mean_abundance_bioproject"
  )

  missing_cols <- setdiff(required_cols, colnames(dt))

  if (length(missing_cols) > 0) {
    stop(
      pipeline_name,
      " table is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dt %>%
    transmute(
      Tax = normalize_genus_name(Tax),
      Pipeline = pipeline_name,
      prevalence = prevalence_bioproject * 100,
      mean_abundance = mean_abundance_bioproject * 100
    )
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(rho = NA_real_, p_value = NA_real_, n = length(x)))
  }

  ct <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE))

  list(
    rho = unname(ct$estimate),
    p_value = ct$p.value,
    n = length(x)
  )
}

calc_jaccard <- function(a_detected, b_detected) {
  union_n <- sum(a_detected | b_detected, na.rm = TRUE)
  intersect_n <- sum(a_detected & b_detected, na.rm = TRUE)
  if (union_n == 0) return(NA_real_)
  intersect_n / union_n
}

make_wide_metric <- function(df, metric_col) {
  df %>%
    select(Tax, Pipeline, all_of(metric_col)) %>%
    pivot_wider(
      names_from = Pipeline,
      values_from = all_of(metric_col),
      values_fill = 0
    ) %>%
    {
      missing_pipelines <- setdiff(pipeline_levels, colnames(.))
      if (length(missing_pipelines) > 0) {
        for (pp in missing_pipelines) .[[pp]] <- 0
      }
      .[, c("Tax", pipeline_levels)]
    }
}

pairwise_concordance <- function(
    df,
    taxa_keep = NULL,
    label = "all_harmonised_genera"
) {

  if (!is.null(taxa_keep)) {
    df <- df %>% filter(Tax %in% taxa_keep)
  }

  pairs <- combn(pipeline_levels, 2, simplify = FALSE)

  res <- lapply(pairs, function(pp) {

    a <- pp[1]
    b <- pp[2]

    # Restrict each comparison to taxa represented in both pipelines.
    # This reproduces the pairwise inner-join logic used in the central
    # harmonisation/concordance analysis.
    taxa_a <- df %>%
      filter(Pipeline == a) %>%
      pull(Tax) %>%
      unique()

    taxa_b <- df %>%
      filter(Pipeline == b) %>%
      pull(Tax) %>%
      unique()

    common_taxa <- intersect(taxa_a, taxa_b)

    pair_df <- df %>%
      filter(
        Pipeline %in% c(a, b),
        Tax %in% common_taxa
      )

    prev_wide <- pair_df %>%
      select(Tax, Pipeline, prevalence) %>%
      pivot_wider(
        names_from = Pipeline,
        values_from = prevalence
      )

    abund_wide <- pair_df %>%
      select(Tax, Pipeline, mean_abundance) %>%
      pivot_wider(
        names_from = Pipeline,
        values_from = mean_abundance
      )

    # Ensure identical taxon order for prevalence and abundance.
    common_order <- intersect(prev_wide$Tax, abund_wide$Tax)

    prev_wide <- prev_wide %>%
      filter(Tax %in% common_order) %>%
      arrange(match(Tax, common_order))

    abund_wide <- abund_wide %>%
      filter(Tax %in% common_order) %>%
      arrange(match(Tax, common_order))

    prev_a <- prev_wide[[a]]
    prev_b <- prev_wide[[b]]
    abund_a <- abund_wide[[a]]
    abund_b <- abund_wide[[b]]

    # Same filter as the central harmonisation script:
    # (prevalence >= 10% in at least one pipeline) AND
    # (mean abundance >= 0.01% in at least one pipeline).
    keep_for_filtered <- (
      (
        prev_a >= concordance_min_prevalence_percent |
          prev_b >= concordance_min_prevalence_percent
      ) &
        (
          abund_a >= concordance_min_abundance_percent |
            abund_b >= concordance_min_abundance_percent
        )
    )

    prev_cor_all <- safe_cor(prev_a, prev_b)
    abund_cor_all <- safe_cor(abund_a, abund_b)

    prev_cor_filtered <- safe_cor(
      prev_a[keep_for_filtered],
      prev_b[keep_for_filtered]
    )

    abund_cor_filtered <- safe_cor(
      abund_a[keep_for_filtered],
      abund_b[keep_for_filtered]
    )

    detected_a <- prev_a > detection_prevalence_threshold
    detected_b <- prev_b > detection_prevalence_threshold

    tibble(
      comparison = paste(a, "vs", b),
      analysis_set = label,
      n_common_genera = length(common_order),
      n_filtered_genera = sum(keep_for_filtered, na.rm = TRUE),
      prevalence_spearman_rho_all = prev_cor_all$rho,
      prevalence_p_value_all = prev_cor_all$p_value,
      prevalence_n_all = prev_cor_all$n,
      abundance_spearman_rho_all = abund_cor_all$rho,
      abundance_p_value_all = abund_cor_all$p_value,
      abundance_n_all = abund_cor_all$n,
      prevalence_spearman_rho_filtered = prev_cor_filtered$rho,
      prevalence_p_value_filtered = prev_cor_filtered$p_value,
      prevalence_n_filtered = prev_cor_filtered$n,
      abundance_spearman_rho_filtered = abund_cor_filtered$rho,
      abundance_p_value_filtered = abund_cor_filtered$p_value,
      abundance_n_filtered = abund_cor_filtered$n,
      detected_jaccard = calc_jaccard(detected_a, detected_b),
      n_detected_in_both = sum(detected_a & detected_b, na.rm = TRUE),
      n_detected_in_any = sum(detected_a | detected_b, na.rm = TRUE)
    )
  })

  bind_rows(res)
}

summarise_concordance_set <- function(concordance_df, set_name) {

  concordance_df %>%
    filter(analysis_set == set_name) %>%
    summarise(
      analysis_set = set_name,
      n_pairwise_comparisons = n(),
      n_common_genera = unique(n_common_genera)[1],
      mean_n_filtered_genera = mean(n_filtered_genera, na.rm = TRUE),
      mean_prevalence_rho_all = mean(prevalence_spearman_rho_all, na.rm = TRUE),
      min_prevalence_rho_all = min(prevalence_spearman_rho_all, na.rm = TRUE),
      max_prevalence_rho_all = max(prevalence_spearman_rho_all, na.rm = TRUE),
      mean_abundance_rho_all = mean(abundance_spearman_rho_all, na.rm = TRUE),
      min_abundance_rho_all = min(abundance_spearman_rho_all, na.rm = TRUE),
      max_abundance_rho_all = max(abundance_spearman_rho_all, na.rm = TRUE),
      mean_prevalence_rho_filtered = mean(prevalence_spearman_rho_filtered, na.rm = TRUE),
      min_prevalence_rho_filtered = min(prevalence_spearman_rho_filtered, na.rm = TRUE),
      max_prevalence_rho_filtered = max(prevalence_spearman_rho_filtered, na.rm = TRUE),
      mean_abundance_rho_filtered = mean(abundance_spearman_rho_filtered, na.rm = TRUE),
      min_abundance_rho_filtered = min(abundance_spearman_rho_filtered, na.rm = TRUE),
      max_abundance_rho_filtered = max(abundance_spearman_rho_filtered, na.rm = TRUE),
      mean_detected_jaccard = mean(detected_jaccard, na.rm = TRUE),
      min_detected_jaccard = min(detected_jaccard, na.rm = TRUE),
      max_detected_jaccard = max(detected_jaccard, na.rm = TRUE),
      .groups = "drop"
    )
}

format_p <- function(x) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "e", digits = 2))
}

# =====================================================
# LOAD PREVALENCE TABLES
# =====================================================

message("Building prevalence summaries from public BioProject matrices...")

build_prevalence_summary <- function(path, pipeline_name) {

  dt <- fread(path)

  metadata_cols <- c(
    "BioProject",
    "CenterName_standardized",
    "Continent",
    "Country",
    "Platform"
  )

  taxa_cols <- setdiff(names(dt), metadata_cols)

  if (!"BioProject" %in% names(dt)) {
    stop(pipeline_name, ": BioProject column missing.")
  }

  # Taxonomic abundance columns must be numeric.
  dt[, (taxa_cols) := lapply(.SD, as.numeric), .SDcols = taxa_cols]

  # Reproduce the prevalence definition used in the publication-level
  # BioProject audit. Each row of the harmonised matrix is retained as
  # generated by the upstream workflow, including platform-level strata.
  # A taxon is considered detected when relative abundance > 0.0001
  # (0.01%).
  detection_threshold <- 0.0001

  mat <- as.matrix(dt[, ..taxa_cols])
  storage.mode(mat) <- "numeric"

  n_bp_rows <- nrow(mat)

  message(
    pipeline_name,
    ": ",
    n_bp_rows,
    " BioProject-level matrix rows; detection threshold = ",
    detection_threshold * 100,
    "%"
  )

  data.table(
    Tax = taxa_cols,
    prevalence_bioproject =
      colSums(mat > detection_threshold, na.rm = TRUE) / n_bp_rows,
    mean_abundance_bioproject =
      colMeans(mat, na.rm = TRUE)
  )
}

qiime_raw <- build_prevalence_summary(
  file.path(
    audit_dir,
    "QIIME2_Amplicon_genus_BioProject_taxa_wide.csv"
  ),
  "QIIME2"
)

kraken_raw <- build_prevalence_summary(
  file.path(
    audit_dir,
    "Kraken2_WGS_genus_BioProject_taxa_wide.csv"
  ),
  "Kraken2"
)

motus_raw <- build_prevalence_summary(
  file.path(
    audit_dir,
    "mOTUs_WGS_genus_BioProject_taxa_wide.csv"
  ),
  "mOTUs"
)

# =====================================================
# STANDARDIZE
# =====================================================

qiime <- standardize_prevalence_table(qiime_raw, "QIIME2")
kraken <- standardize_prevalence_table(kraken_raw, "Kraken2")
motus <- standardize_prevalence_table(motus_raw, "mOTUs")

# =====================================================
# MERGE PIPELINES
# =====================================================

all_prev <- bind_rows(qiime, kraken, motus) %>%
  filter(
    !is.na(Tax),
    Tax != "",
    !Tax %in% excluded_taxa
  ) %>%
  group_by(Tax, Pipeline) %>%
  summarise(
    # After genus-name harmonisation, duplicate genera can occur.
    # Use max prevalence/abundance to preserve the detected signal.
    prevalence = max(prevalence, na.rm = TRUE),
    mean_abundance = max(mean_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Group = assign_biological_group(Tax),
    Pipeline = factor(Pipeline, levels = pipeline_levels)
  )

message("All harmonised genera: ", n_distinct(all_prev$Tax))

# =====================================================
# SUPPLEMENTARY TABLES: ALL GENERA
# =====================================================

supp_all_long <- all_prev %>%
  mutate(
    Pipeline = as.character(Pipeline),
    Group = as.character(Group)
  ) %>%
  arrange(Group, Tax, Pipeline)

fwrite(
  supp_all_long,
  file.path(out_dir, "Supplementary_Table_S1_all_harmonised_genera_long.csv")
)

supp_prev_matrix <- supp_all_long %>%
  select(Group, Tax, Pipeline, prevalence) %>%
  pivot_wider(names_from = Pipeline, values_from = prevalence, values_fill = 0) %>%
  {
    missing_pipelines <- setdiff(pipeline_levels, colnames(.))
    if (length(missing_pipelines) > 0) {
      for (pp in missing_pipelines) .[[pp]] <- 0
    }
    .
  } %>%
  mutate(
    max_prevalence = pmax(!!!rlang::syms(pipeline_levels), na.rm = TRUE),
    mean_prevalence = rowMeans(across(all_of(pipeline_levels)), na.rm = TRUE),
    n_pipelines_detected = rowSums(across(all_of(pipeline_levels), ~ .x > detection_prevalence_threshold), na.rm = TRUE)
  ) %>%
  arrange(desc(n_pipelines_detected), desc(max_prevalence), Tax)

fwrite(
  supp_prev_matrix,
  file.path(out_dir, "Supplementary_Table_S1_all_genera_prevalence_matrix.csv")
)

supp_abund_matrix <- supp_all_long %>%
  select(Group, Tax, Pipeline, mean_abundance) %>%
  pivot_wider(names_from = Pipeline, values_from = mean_abundance, values_fill = 0) %>%
  {
    missing_pipelines <- setdiff(pipeline_levels, colnames(.))
    if (length(missing_pipelines) > 0) {
      for (pp in missing_pipelines) .[[pp]] <- 0
    }
    .
  } %>%
  mutate(
    max_mean_abundance = pmax(!!!rlang::syms(pipeline_levels), na.rm = TRUE),
    mean_mean_abundance = rowMeans(across(all_of(pipeline_levels)), na.rm = TRUE),
    n_pipelines_detected = rowSums(across(all_of(pipeline_levels), ~ .x > 0), na.rm = TRUE)
  ) %>%
  arrange(desc(n_pipelines_detected), desc(max_mean_abundance), Tax)

fwrite(
  supp_abund_matrix,
  file.path(out_dir, "Supplementary_Table_S1_all_genera_mean_abundance_matrix.csv")
)

# =====================================================
# TOP GENERA
# selection by maximum prevalence across pipelines
# =====================================================

top_genera <- all_prev %>%
  group_by(Tax, Group) %>%
  summarise(
    max_prev = max(prevalence, na.rm = TRUE),
    mean_prev = mean(prevalence, na.rm = TRUE),
    max_mean_abundance = max(mean_abundance, na.rm = TRUE),
    mean_mean_abundance = mean(mean_abundance, na.rm = TRUE),
    n_pipelines_detected = sum(prevalence > detection_prevalence_threshold, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    desc(max_prev),
    desc(mean_prev),
    desc(max_mean_abundance)
  ) %>%
  slice_head(n = top_n_genera)

fwrite(
  top_genera,
  file.path(out_dir, "Figure3_selected_top30_genera.csv")
)

top_genera_other <- top_genera %>%
  filter(Group == "Other")

fwrite(
  top_genera_other,
  file.path(out_dir, "Figure3_selected_top30_genera_other_not_plotted.csv")
)

# =====================================================
# HEATMAP DATA
# only selected top genera that belong to predefined biological groups are plotted
# =====================================================

heatmap_df <- all_prev %>%
  filter(
    Tax %in% top_genera$Tax,
    Group %in% group_order
  ) %>%
  mutate(
    Group = factor(Group, levels = group_order)
  )

plotted_taxa <- heatmap_df %>%
  pull(Tax) %>%
  unique()

tax_order <- heatmap_df %>%
  group_by(Group, Tax) %>%
  summarise(mean_prev = mean(prevalence, na.rm = TRUE), .groups = "drop") %>%
  arrange(Group, mean_prev) %>%
  pull(Tax)

heatmap_df <- heatmap_df %>%
  mutate(
    Tax = factor(Tax, levels = unique(tax_order))
  )

message("Selected top genera: ", n_distinct(top_genera$Tax))
message("Plotted biological-group genera: ", n_distinct(heatmap_df$Tax))

# =====================================================
# CONCORDANCE STATISTICS ACROSS PIPELINES
# =====================================================

# Use the canonical BioProject-level prevalence summaries generated by
# scripts/02_harmonization/01_build_bioproject_taxonomic_matrices.R.
# This ensures that Figure 3 and the central concordance analysis use
# exactly the same harmonised genus sets and BioProject summaries.

canonical_prevalence_files <- c(
  QIIME2 = file.path(
    audit_dir,
    "QIIME2_Amplicon_genus_BioProject_prevalence_summary.csv"
  ),
  Kraken2 = file.path(
    audit_dir,
    "Kraken2_WGS_genus_BioProject_prevalence_summary.csv"
  ),
  mOTUs = file.path(
    audit_dir,
    "mOTUs_WGS_genus_BioProject_prevalence_summary.csv"
  )
)

missing_canonical <- canonical_prevalence_files[
  !file.exists(canonical_prevalence_files)
]

if (length(missing_canonical) > 0) {
  stop(
    "Missing canonical BioProject prevalence summaries: ",
    paste(missing_canonical, collapse = ", ")
  )
}

concordance_input <- purrr::imap_dfr(
  canonical_prevalence_files,
  function(path, pipeline_name) {

    fread(path) %>%
      as_tibble() %>%
      transmute(
        Tax = Tax,
        Pipeline = pipeline_name,

        # pairwise_concordance() uses percentage units
        prevalence = prevalence_bioproject * 100,
        mean_abundance = mean_abundance_bioproject * 100
      )
  }
)

message(
  "Canonical concordance input: ",
  n_distinct(concordance_input$Tax),
  " union genera across pipelines."
)

# ------------------------------------------------------------
# Concordance analysis
#
# QIIME2-vs-WGS comparisons use the complete method-specific
# BioProject collections and therefore describe aggregate
# cross-collection taxonomic-pattern concordance.
#
# Kraken2-vs-mOTUs is calculated separately on the identical
# set of 19 WGS BioProjects processed successfully by both
# classifiers (paired method comparison).
# ------------------------------------------------------------

paired_prevalence_files <- c(
  Kraken2 = file.path(
    audit_dir,
    "Kraken2_WGS_genus_paired19_BioProject_prevalence_summary.csv"
  ),
  mOTUs = file.path(
    audit_dir,
    "mOTUs_WGS_genus_paired19_BioProject_prevalence_summary.csv"
  )
)

missing_paired <- paired_prevalence_files[
  !file.exists(paired_prevalence_files)
]

if (length(missing_paired) > 0) {
  stop(
    "Missing paired Kraken2-mOTUs prevalence summaries: ",
    paste(missing_paired, collapse = ", ")
  )
}

paired_concordance_input <- purrr::imap_dfr(
  paired_prevalence_files,
  function(path, pipeline_name) {

    fread(path) %>%
      as_tibble() %>%
      transmute(
        Tax = Tax,
        Pipeline = pipeline_name,
        prevalence = prevalence_bioproject * 100,
        mean_abundance = mean_abundance_bioproject * 100
      )
  }
)

replace_k2_motus_with_paired <- function(
    standard_df,
    taxa_keep = NULL,
    label
) {

  paired_df <- pairwise_concordance(
    paired_concordance_input,
    taxa_keep = taxa_keep,
    label = label
  )

  bind_rows(
    standard_df %>%
      filter(comparison != "Kraken2 vs mOTUs"),
    paired_df %>%
      filter(comparison == "Kraken2 vs mOTUs")
  )
}


# ---- all harmonised genera ----

concordance_all_standard <- pairwise_concordance(
  concordance_input,
  taxa_keep = NULL,
  label = "all_harmonised_genera"
)

concordance_all <- replace_k2_motus_with_paired(
  concordance_all_standard,
  taxa_keep = NULL,
  label = "all_harmonised_genera"
)


# ---- selected top 30 genera ----

concordance_top30_standard <- pairwise_concordance(
  concordance_input,
  taxa_keep = top_genera$Tax,
  label = "selected_top30_genera"
)

concordance_top30 <- replace_k2_motus_with_paired(
  concordance_top30_standard,
  taxa_keep = top_genera$Tax,
  label = "selected_top30_genera"
)


# ---- plotted biological-group genera ----

concordance_plotted_standard <- pairwise_concordance(
  concordance_input,
  taxa_keep = plotted_taxa,
  label = "plotted_biological_group_genera"
)

concordance_plotted <- replace_k2_motus_with_paired(
  concordance_plotted_standard,
  taxa_keep = plotted_taxa,
  label = "plotted_biological_group_genera"
)

concordance_combined <- bind_rows(
  concordance_all,
  concordance_top30,
  concordance_plotted
) %>%
  mutate(
    prevalence_p_value_all_scientific = format_p(prevalence_p_value_all),
    abundance_p_value_all_scientific = format_p(abundance_p_value_all),
    prevalence_p_value_filtered_scientific = format_p(prevalence_p_value_filtered),
    abundance_p_value_filtered_scientific = format_p(abundance_p_value_filtered)
  )

fwrite(
  concordance_all,
  file.path(out_dir, "Figure3_pipeline_concordance_statistics_all_genera.csv")
)

fwrite(
  concordance_top30,
  file.path(out_dir, "Figure3_pipeline_concordance_statistics_selected_top30_genera.csv")
)

fwrite(
  concordance_plotted,
  file.path(out_dir, "Figure3_pipeline_concordance_statistics_plotted_biological_group_genera.csv")
)

fwrite(
  concordance_combined,
  file.path(out_dir, "Figure3_pipeline_concordance_statistics.csv")
)

concordance_summary <- bind_rows(
  summarise_concordance_set(concordance_combined, "all_harmonised_genera"),
  summarise_concordance_set(concordance_combined, "selected_top30_genera"),
  summarise_concordance_set(concordance_combined, "plotted_biological_group_genera")
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

fwrite(
  concordance_summary,
  file.path(out_dir, "Supplementary_Table_S2_pipeline_concordance_summary.csv")
)

# Long metric form, convenient for manuscript tables.
concordance_summary_long <- concordance_summary %>%
  pivot_longer(
    cols = -analysis_set,
    names_to = "metric",
    values_to = "value"
  )

fwrite(
  concordance_summary_long,
  file.path(out_dir, "Supplementary_Table_S2_pipeline_concordance_summary_long.csv")
)

# =====================================================
# HEATMAP
# =====================================================

p <- ggplot(
  heatmap_df,
  aes(
    x = Pipeline,
    y = Tax,
    fill = prevalence
  )
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f", prevalence)), size = 3) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradientn(
    colours = c("white", "#d9f0a3", "#78c679", "#238443"),
    limits = c(0, 100),
    name = "Prevalence (%)"
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(face = "italic", size = 10),
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.8, "cm"),
    strip.text.y = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey90", colour = NA),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# =====================================================
# SAVE FIGURE
# =====================================================

ggsave(
  file.path(out_dir, "Figure3_consensus_genera_heatmap_grouped.png"),
  p,
  width = 10,
  height = 13,
  dpi = 600
)

ggsave(
  file.path(out_dir, "Figure3_consensus_genera_heatmap_grouped.pdf"),
  p,
  width = 10,
  height = 13
)

# Keep legacy filenames too.
ggsave(
  file.path(out_dir, "Figure3_consensus_genera_heatmap.png"),
  p,
  width = 10,
  height = 13,
  dpi = 600
)

ggsave(
  file.path(out_dir, "Figure3_consensus_genera_heatmap.pdf"),
  p,
  width = 8.5,
  height = 11
)

# =====================================================
# HEATMAP MATRIX
# =====================================================

heatmap_matrix <- heatmap_df %>%
  mutate(
    Tax = as.character(Tax),
    Pipeline = as.character(Pipeline),
    Group = as.character(Group)
  ) %>%
  select(Group, Tax, Pipeline, prevalence) %>%
  pivot_wider(names_from = Pipeline, values_from = prevalence, values_fill = 0) %>%
  arrange(factor(Group, levels = group_order), Tax)

fwrite(
  heatmap_matrix,
  file.path(out_dir, "Figure3_heatmap_matrix.csv")
)

# =====================================================
# MEAN ABUNDANCE MATRIX
# =====================================================

abundance_matrix <- heatmap_df %>%
  mutate(
    Tax = as.character(Tax),
    Pipeline = as.character(Pipeline),
    Group = as.character(Group)
  ) %>%
  select(Group, Tax, Pipeline, mean_abundance) %>%
  pivot_wider(names_from = Pipeline, values_from = mean_abundance, values_fill = 0) %>%
  arrange(factor(Group, levels = group_order), Tax)

fwrite(
  abundance_matrix,
  file.path(out_dir, "Figure3_mean_abundance_matrix.csv")
)

# =====================================================
# EXPORT BIOLOGICAL GROUPS
# =====================================================

core_matrix <- heatmap_matrix %>% filter(Group == "Cheese-associated genera")
rind_matrix <- heatmap_matrix %>% filter(Group == "Rind-associated microbiota")
environment_matrix <- heatmap_matrix %>% filter(Group == "Environmental microbiota")

fwrite(core_matrix, file.path(out_dir, "Figure3_core_cheese_microbiota.csv"))
fwrite(rind_matrix, file.path(out_dir, "Figure3_rind_associated_microbiota.csv"))
fwrite(environment_matrix, file.path(out_dir, "Figure3_environmental_microbiota.csv"))

# =====================================================
# EXPORT LONG TABLE FOR PLOTTED FIGURE
# =====================================================

heatmap_table <- heatmap_df %>%
  mutate(
    Tax = as.character(Tax),
    Pipeline = as.character(Pipeline),
    Group = as.character(Group)
  ) %>%
  select(Group, Tax, Pipeline, prevalence, mean_abundance) %>%
  arrange(factor(Group, levels = group_order), Tax, Pipeline)

fwrite(
  heatmap_table,
  file.path(out_dir, "Figure3_heatmap_source_data.csv")
)

# =====================================================
# SHORT TEXT SUMMARY FOR RESULTS
# =====================================================

get_summary_value <- function(set_name, metric_name) {
  concordance_summary %>%
    filter(analysis_set == set_name) %>%
    pull(all_of(metric_name))
}

summary_table <- tibble(
  metric = c(
    "All harmonised genera",
    "Selected top genera",
    "Plotted biological-group genera",
    "All genera mean filtered prevalence Spearman rho",
    "All genera minimum filtered prevalence Spearman rho",
    "All genera mean filtered abundance Spearman rho",
    "All genera mean detected-genera Jaccard",
    "Top30 mean filtered prevalence Spearman rho",
    "Top30 minimum filtered prevalence Spearman rho",
    "Top30 mean filtered abundance Spearman rho",
    "Top30 mean detected-genera Jaccard",
    "Plotted genera mean filtered prevalence Spearman rho",
    "Plotted genera minimum filtered prevalence Spearman rho",
    "Plotted genera mean filtered abundance Spearman rho",
    "Plotted genera mean detected-genera Jaccard"
  ),
  value = c(
    n_distinct(all_prev$Tax),
    n_distinct(top_genera$Tax),
    n_distinct(heatmap_df$Tax),
    get_summary_value("all_harmonised_genera", "mean_prevalence_rho_filtered"),
    get_summary_value("all_harmonised_genera", "min_prevalence_rho_filtered"),
    get_summary_value("all_harmonised_genera", "mean_abundance_rho_filtered"),
    get_summary_value("all_harmonised_genera", "mean_detected_jaccard"),
    get_summary_value("selected_top30_genera", "mean_prevalence_rho_filtered"),
    get_summary_value("selected_top30_genera", "min_prevalence_rho_filtered"),
    get_summary_value("selected_top30_genera", "mean_abundance_rho_filtered"),
    get_summary_value("selected_top30_genera", "mean_detected_jaccard"),
    get_summary_value("plotted_biological_group_genera", "mean_prevalence_rho_filtered"),
    get_summary_value("plotted_biological_group_genera", "min_prevalence_rho_filtered"),
    get_summary_value("plotted_biological_group_genera", "mean_abundance_rho_filtered"),
    get_summary_value("plotted_biological_group_genera", "mean_detected_jaccard")
  )
) %>%
  mutate(value = ifelse(is.numeric(value), round(as.numeric(value), 4), value))

fwrite(
  summary_table,
  file.path(out_dir, "Figure3_RESULTS_SUMMARY.tsv"),
  sep = "\t"
)

# Human-readable one-page summary.
results_text <- c(
  "Figure 3 results summary",
  "========================",
  "",
  paste0("All harmonised genera: ", n_distinct(all_prev$Tax)),
  paste0("Selected top genera: ", n_distinct(top_genera$Tax)),
  paste0("Plotted biological-group genera: ", n_distinct(heatmap_df$Tax)),
  "",
  "Concordance summary:",
  capture.output(print(concordance_summary)),
  "",
  "Pairwise concordance statistics:",
  capture.output(print(concordance_combined))
)

writeLines(
  results_text,
  file.path(out_dir, "Figure3_RESULTS_TEXT.txt")
)

# =====================================================
# README
# =====================================================

summary_lines <- c(
  "Figure 3 completed.",
  "Version: v4_all_and_plotted_concordance",
  paste0("Output directory: ", out_dir),
  paste0("Audit directory: ", audit_dir),
  paste0("Top genera requested: ", top_n_genera),
  paste0("All harmonised genera exported: ", n_distinct(all_prev$Tax)),
  paste0("Selected top genera: ", n_distinct(top_genera$Tax)),
  paste0("Plotted genera: ", n_distinct(heatmap_df$Tax)),
  "",
  "Biological groups plotted:",
  "- Cheese-associated genera",
  "- Rind-associated microbiota",
  "- Environmental microbiota",
  "",
  "Concordance statistics:",
  "- pairwise Spearman correlation of BioProject prevalence",
  "- pairwise Spearman correlation of mean relative abundance",
  "- pairwise Jaccard index of detected genera",
  "- statistics exported for all genera, selected top 30 genera, and plotted biological-group genera",
  paste0("- filtered concordance threshold: prevalence >= ", concordance_min_prevalence_percent, "% in at least one pipeline AND mean abundance >= ", concordance_min_abundance_percent, "% in at least one pipeline"),
  "",
  "Excluded taxa:",
  paste(excluded_taxa, collapse = ", "),
  "",
  "Main output files:",
  "- Figure3_consensus_genera_heatmap_grouped.png",
  "- Figure3_consensus_genera_heatmap_grouped.pdf",
  "- Figure3_heatmap_matrix.csv",
  "- Figure3_mean_abundance_matrix.csv",
  "- Figure3_heatmap_source_data.csv",
  "- Figure3_core_cheese_microbiota.csv",
  "- Figure3_rind_associated_microbiota.csv",
  "- Figure3_environmental_microbiota.csv",
  "- Figure3_selected_top30_genera.csv",
  "- Figure3_selected_top30_genera_other_not_plotted.csv",
  "- Figure3_pipeline_concordance_statistics.csv",
  "- Figure3_pipeline_concordance_statistics_all_genera.csv",
  "- Figure3_pipeline_concordance_statistics_selected_top30_genera.csv",
  "- Figure3_pipeline_concordance_statistics_plotted_biological_group_genera.csv",
  "- Supplementary_Table_S2_pipeline_concordance_summary.csv",
  "- Supplementary_Table_S2_pipeline_concordance_summary_long.csv",
  "- Figure3_RESULTS_SUMMARY.tsv",
  "- Figure3_RESULTS_TEXT.txt",
  "",
  "Supplementary output files:",
  "- Supplementary_Table_S1_all_harmonised_genera_long.csv",
  "- Supplementary_Table_S1_all_genera_prevalence_matrix.csv",
  "- Supplementary_Table_S1_all_genera_mean_abundance_matrix.csv"
)

writeLines(
  summary_lines,
  file.path(out_dir, "README_Figure3.txt")
)

message(paste(summary_lines, collapse = "\n"))
message("DONE")
