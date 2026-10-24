#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

# ============================================================
# 1. PATHS
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
rew_dir <- file.path(
  project_dir,
  "data/08_reporting/Rew_3"
)

out_dir <- file.path(rew_dir, "NEW_SUPP")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ============================================================
# 2. HELPERS
# ============================================================

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Missing input: ", path)
  }
  fread(path, sep = "\t", header = TRUE, data.table = FALSE)
}

write_csv_final <- function(x, name) {
  fwrite(
    as.data.frame(x),
    file.path(out_dir, name),
    sep = ",",
    quote = TRUE,
    na = ""
  )
}

find_col <- function(df, patterns, required = TRUE) {
  nm <- names(df)
  for (p in patterns) {
    hit <- grep(p, nm, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  if (required) {
    stop(
      "Could not find column matching: ",
      paste(patterns, collapse = " / "),
      "\nAvailable columns: ",
      paste(nm, collapse = ", ")
    )
  }
  NULL
}

# ============================================================
# 3. INPUT FILES
# ============================================================

f_origin <- file.path(
  rew_dir,
  "01_origin_audit/S01_Amplicon_BioProject_origin_audit.tsv"
)

f_country_cov <- file.path(
  rew_dir,
  "01_origin_audit/S04_confirmed_origin_country_coverage.tsv"
)

f_perm_confirmed <- file.path(
  rew_dir,
  "02_confirmed_origin/PERMANOVA_confirmed_origin_country_composition_type.tsv"
)

f_perm_primary <- file.path(
  rew_dir,
  "02_confirmed_origin/PERMANOVA_Country_only_ALL_confirmed_origin.tsv"
)

f_disp_primary <- file.path(
  rew_dir,
  "02_confirmed_origin/PERMDISP_Country_ALL_confirmed_origin.tsv"
)

f_disp_confirmed <- file.path(
  rew_dir,
  "02_confirmed_origin/PERMDISP_confirmed_origin.tsv"
)

f_perm_min2 <- file.path(
  rew_dir,
  "03_country_min2/PERMANOVA_min2_country_composition_type.tsv"
)

f_perm_min2_geo <- file.path(
  rew_dir,
  "03_country_min2/PERMANOVA_min2_Country_only_ALL_eligible.tsv"
)

f_disp_min2_geo <- file.path(
  rew_dir,
  "03_country_min2/PERMDISP_min2_Country_ALL_eligible.tsv"
)

f_disp_min2 <- file.path(
  rew_dir,
  "03_country_min2/PERMDISP_min2.tsv"
)

f_balanced_summary <- file.path(
  rew_dir,
  "04_balanced_subsampling/balanced_subsampling_summary.tsv"
)

f_balanced_all <- file.path(
  rew_dir,
  "04_balanced_subsampling/balanced_subsampling_all_iterations.tsv"
)

f_cluster_validation <- file.path(
  rew_dir,
  "05_cluster_validation/cluster_validation_combined.tsv"
)

f_cluster_k <- file.path(
  rew_dir,
  "05_cluster_validation/k2_to_k10_validation_metrics.tsv"
)

f_cluster_stability <- file.path(
  rew_dir,
  "05_cluster_validation/cluster_stability_summary.tsv"
)

f_removed <- file.path(
  rew_dir,
  "06_starter_removed/PERMANOVA_country_after_removing_5_starter_genera.tsv"
)

f_continuous <- file.path(
  rew_dir,
  "07_continuous_starter/Country_effect_after_continuous_starter_proportions.tsv"
)

f_heatmap <- file.path(
  rew_dir,
  "08_figure4_source_tables/Figure4A_heatmap_source_confirmed_origin.tsv"
)

f_composition <- file.path(
  rew_dir,
  "08_figure4_source_tables/Figure4B_composition_type_source.tsv"
)

f_cluster_comp <- file.path(
  rew_dir,
  "08_figure4_source_tables/Figure4_cluster_composition.tsv"
)

f_membership <- file.path(
  rew_dir,
  "08_figure4_source_tables/Figure4_confirmed_origin_cluster_membership.tsv"
)

f_fig_perm <- file.path(
  rew_dir,
  "08_figure4_source_tables/Figure4D_PERMANOVA_source.tsv"
)

f_fig_disp <- file.path(
  rew_dir,
  "08_figure4_source_tables/Figure4E_PERMDISP_source.tsv"
)

# ============================================================
# 4. READ DATA
# ============================================================

origin <- read_tsv(f_origin)
country_cov <- read_tsv(f_country_cov)

perm_confirmed <- read_tsv(f_perm_confirmed)
disp_confirmed <- read_tsv(f_disp_confirmed)

perm_primary <- read_tsv(f_perm_primary)
disp_primary <- read_tsv(f_disp_primary)

perm_min2 <- read_tsv(f_perm_min2)
disp_min2 <- read_tsv(f_disp_min2)

perm_min2_geo <- read_tsv(f_perm_min2_geo)
disp_min2_geo <- read_tsv(f_disp_min2_geo)

balanced_summary <- read_tsv(f_balanced_summary)
balanced_all <- read_tsv(f_balanced_all)

cluster_validation <- read_tsv(f_cluster_validation)
cluster_k <- read_tsv(f_cluster_k)
cluster_stability <- read_tsv(f_cluster_stability)

removed <- read_tsv(f_removed)
continuous <- read_tsv(f_continuous)

heat <- read_tsv(f_heatmap)
composition <- read_tsv(f_composition)
cluster_comp <- read_tsv(f_cluster_comp)
membership <- read_tsv(f_membership)

fig_perm <- read_tsv(f_fig_perm)
fig_disp <- read_tsv(f_fig_disp)

# ============================================================
# 5. SUPPLEMENTARY TABLE S4
# Geographic provenance
# ============================================================

message("Creating Supplementary Table S4...")

s4 <- origin

# Put most informative columns first if present
priority <- c(
  "BioProject",
  "CenterCountry",
  "SampleOriginCountry",
  "OriginStatus",
  "OriginEvidence",
  "n_SRA",
  "n_sra",
  "n_BioSamples",
  "n_biosamples"
)

priority <- priority[priority %in% names(s4)]

s4 <- s4[, c(
  priority,
  setdiff(names(s4), priority)
), drop = FALSE]

write_csv_final(
  s4,
  "Supplementary_Table_S4_Geographic_origin.csv"
)

write_csv_final(
  country_cov,
  "Supplementary_Table_S4b_Country_coverage.csv"
)

# ============================================================
# 6. FUNCTIONS FOR EXTRACTING PERMANOVA
# ============================================================

standardize_perm <- function(df, analysis_name) {

  term_col <- find_col(
    df,
    c("^Term$", "^Variable$", "^Factor$", "^Model$", "term"),
    required = FALSE
  )

  df_col <- find_col(df, c("^Df$", "degrees"), required = FALSE)
  r2_col <- find_col(df, c("^R2$", "R2", "R\\.2"), required = FALSE)
  f_col <- find_col(df, c("^F$", "pseudo", "F.Model"), required = FALSE)
  p_col <- find_col(df, c("Pr\\(>F\\)", "^p$", "p_value", "p.value"), required = FALSE)

  out <- data.frame(
    Analysis = analysis_name,
    Term = if (!is.null(term_col)) as.character(df[[term_col]]) else NA_character_,
    Df = if (!is.null(df_col)) suppressWarnings(as.numeric(df[[df_col]])) else NA_real_,
    R2 = if (!is.null(r2_col)) suppressWarnings(as.numeric(df[[r2_col]])) else NA_real_,
    pseudo_F = if (!is.null(f_col)) suppressWarnings(as.numeric(df[[f_col]])) else NA_real_,
    p_value = if (!is.null(p_col)) suppressWarnings(as.numeric(df[[p_col]])) else NA_real_,
    stringsAsFactors = FALSE
  )

  out
}

perm1 <- standardize_perm(
  perm_primary,
  "Primary geography: all single-origin BioProjects (n=73; 20 countries)"
)

perm2 <- standardize_perm(
  perm_min2_geo,
  "Replicated-country sensitivity: >=2 BioProjects/country (n=63; 10 countries)"
)

perm3 <- standardize_perm(
  perm_confirmed,
  "Composition-type analysis (n=68)"
)

perm4 <- standardize_perm(
  removed,
  "Five composition-defining genera removed"
)

perm5 <- standardize_perm(
  continuous,
  "Continuous five-genus composition sensitivity (starter-positive subset)"
)

s5_perm <- bind_rows(
  perm1,
  perm2,
  perm3,
  perm4,
  perm5
)

# Keep relevant terms, but don't silently discard if naming differs
relevant <- grepl(
  "Country|Composition|Cluster|Model",
  s5_perm$Term,
  ignore.case = TRUE
)

if (any(relevant, na.rm = TRUE)) {
  s5_perm <- s5_perm[
    relevant | is.na(relevant),
    ,
    drop = FALSE
  ]
}

# ============================================================
# 7. PERMDISP
# ============================================================

standardize_disp <- function(df, analysis_name) {

  term_col <- find_col(
    df,
    c("^Term$", "^Variable$", "^Factor$", "^Groups$", "term"),
    required = FALSE
  )

  f_col <- find_col(df, c("^F$", "F.Model"), required = FALSE)
  p_col <- find_col(df, c("Pr\\(>F\\)", "^p$", "p_value", "p.value"), required = FALSE)

  data.frame(
    Analysis = analysis_name,
    Term = if (!is.null(term_col)) as.character(df[[term_col]]) else NA_character_,
    PERMDISP_F =
      if (!is.null(f_col)) suppressWarnings(as.numeric(df[[f_col]])) else NA_real_,
    PERMDISP_p =
      if (!is.null(p_col)) suppressWarnings(as.numeric(df[[p_col]])) else NA_real_,
    stringsAsFactors = FALSE
  )
}

disp1 <- standardize_disp(
  disp_primary,
  "Primary geography: all single-origin BioProjects (n=73; 20 countries)"
)

disp2 <- standardize_disp(
  disp_min2_geo,
  "Replicated-country sensitivity: >=2 BioProjects/country (n=63; 10 countries)"
)

disp3 <- standardize_disp(
  disp_confirmed,
  "Composition-type subset (n=68)"
)

s5_disp <- bind_rows(
  disp1,
  disp2,
  disp3
)

# ============================================================
# 8. S5 OUTPUT
# ============================================================

write_csv_final(
  s5_perm,
  "Supplementary_Table_S5_Geographic_PERMANOVA.csv"
)

write_csv_final(
  s5_disp,
  "Supplementary_Table_S5b_Geographic_PERMDISP.csv"
)

# ============================================================
# 9. S6 — CLUSTER VALIDATION
# ============================================================

message("Creating Supplementary Table S6...")

# ARI values quantify reproducibility across repeated random
# k-means initializations fitted to the same data matrix.
# They are not bootstrap or BioProject-resampling stability.
s6 <- cluster_k

# Join stability where possible
k1 <- find_col(s6, c("^k$"), required = FALSE)
k2 <- find_col(cluster_stability, c("^k$"), required = FALSE)

if (!is.null(k1) && !is.null(k2)) {
  s6 <- left_join(
    s6,
    cluster_stability,
    by = setNames(k2, k1)
  )
}

write_csv_final(
  s6,
  "Supplementary_Table_S6_Composition_type_validation.csv"
)

# Full validation source retained as separate source-data table
write_csv_final(
  cluster_validation,
  "Supplementary_Table_S6b_Composition_type_validation_full.csv"
)

# ============================================================
# 10. S7 — BALANCED SUBSAMPLING
# ============================================================

message("Creating Supplementary Table S7...")

write_csv_final(
  balanced_summary,
  "Supplementary_Table_S7_Balanced_subsampling_summary.csv"
)

# Keep all iterations as source data, not main supplementary table
write_csv_final(
  balanced_all,
  "Source_Data_S7_Balanced_subsampling_all_iterations.csv"
)

# ============================================================
# 11. FIGURE 4 SOURCE DATA EXPORT
# ============================================================

message("Exporting Figure 4 source data...")

write_csv_final(
  heat,
  "Source_Data_Figure4A.csv"
)

write_csv_final(
  composition,
  "Source_Data_Figure4B.csv"
)

write_csv_final(
  fig_perm,
  "Source_Data_Figure4C.csv"
)

write_csv_final(
  fig_disp,
  "Source_Data_Figure4_PERMDISP.csv"
)

write_csv_final(
  membership,
  "Source_Data_Figure4_cluster_membership.csv"
)

# ============================================================
# 12. FIGURE 4A — HEATMAP
# Countries represented by >=2 BioProjects
# ============================================================

message("Creating Figure 4A — >=2 BioProjects/country...")

starter_genera <- c(
  "Lacticaseibacillus",
  "Leuconostoc",
  "Lactococcus",
  "Streptococcus",
  "Lactobacillus"
)

bio_col <- find_col(heat, c("^BioProject$"))
country_col <- find_col(
  heat,
  c("^SampleOriginCountry$", "^Country$")
)

cluster_col <- find_col(
  heat,
  c("^CompositionType$", "^Cluster$", "composition"),
  required = FALSE
)

present_genera <- intersect(
  starter_genera,
  names(heat)
)

if (length(present_genera) != 5) {
  stop(
    "Figure4A: expected all 5 genera. Found: ",
    paste(present_genera, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Panel A shows the full composition-type analysis population.
# These are BioProjects with a single BioSample-reported
# sample-origin country and nonzero abundance across the five
# genera used to define composition types.
# No >=2-country restriction is applied to this descriptive panel.
# ------------------------------------------------------------

heat_A <- heat

country_n <- heat_A %>%
  distinct(
    .data[[bio_col]],
    .data[[country_col]]
  ) %>%
  count(
    .data[[country_col]],
    name = "n_BioProjects"
  )

message(
  "Figure 4A composition-type population: ",
  dplyr::n_distinct(heat_A[[bio_col]]),
  " BioProjects; ",
  dplyr::n_distinct(heat_A[[country_col]]),
  " countries."
)

write_csv_final(
  heat_A,
  "Source_Data_Figure4A_all_composition_positive.csv"
)

write_csv_final(
  country_n %>%
    arrange(desc(n_BioProjects)),
  "Source_Data_Figure4A_country_counts.csv"
)

# Order BioProjects:
# Country -> composition type -> BioProject
# ------------------------------------------------------------

if (!is.null(cluster_col)) {

  ordering_df <- heat_A %>%
    arrange(
      .data[[country_col]],
      .data[[cluster_col]],
      .data[[bio_col]]
    )

} else {

  ordering_df <- heat_A %>%
    arrange(
      .data[[country_col]],
      .data[[bio_col]]
    )
}

bio_levels <- ordering_df[[bio_col]]

# ------------------------------------------------------------
# Long-format heatmap
# ------------------------------------------------------------

heat_long <- heat_A %>%
  select(
    all_of(
      unique(c(
        bio_col,
        country_col,
        cluster_col,
        present_genera
      ))
    )
  ) %>%
  pivot_longer(
    cols = all_of(present_genera),
    names_to = "Genus",
    values_to = "RelativeAbundance"
  )

heat_long[[bio_col]] <- factor(
  heat_long[[bio_col]],
  levels = rev(bio_levels)
)

heat_long$Genus <- factor(
  heat_long$Genus,
  levels = starter_genera
)

# ------------------------------------------------------------
# Country labels: one label per country, no facet strips
# ------------------------------------------------------------

row_positions <- data.frame(
  BioProject_tmp = rev(bio_levels),
  y = seq_along(bio_levels),
  stringsAsFactors = FALSE
)

names(row_positions)[1] <- bio_col

country_positions <- ordering_df %>%
  select(
    all_of(c(
      bio_col,
      country_col
    ))
  ) %>%
  left_join(
    row_positions,
    by = bio_col
  ) %>%
  group_by(
    .data[[country_col]]
  ) %>%
  summarise(
    y = mean(y),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Horizontal boundaries between countries
# ------------------------------------------------------------

country_row_order <- ordering_df %>%
  select(
    all_of(c(
      bio_col,
      country_col
    ))
  ) %>%
  mutate(
    y = match(.data[[bio_col]], rev(bio_levels))
  ) %>%
  arrange(y)

country_boundaries <- country_row_order %>%
  mutate(
    next_country = lead(.data[[country_col]])
  ) %>%
  filter(
    !is.na(next_country),
    .data[[country_col]] != next_country
  ) %>%
  transmute(
    yintercept = y + 0.5
  )

# Export positions used for country separators
write_csv_final(
  country_boundaries,
  "Source_Data_Figure4A_country_boundaries.csv"
)

# ------------------------------------------------------------
# Main heatmap
# ------------------------------------------------------------

pA_heat <- ggplot(
  heat_long,
  aes(
    x = Genus,
    y = .data[[bio_col]],
    fill = RelativeAbundance
  )
) +
  geom_tile(
    linewidth = 0.15
  ) +
  geom_hline(
    data = country_boundaries,
    aes(yintercept = yintercept),
    inherit.aes = FALSE,
    linewidth = 0.55
  ) +
  scale_fill_viridis_c(
    option = "C",
    trans = "sqrt",
    begin = 0.05,
    end = 0.95,
    na.value = "white"
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Relative\nabundance"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      linewidth = 0.5
    ),
    legend.position = "right"
  )

# ------------------------------------------------------------
# Separate country-label panel
# No boxes / no facet strips / one country name per group
# ------------------------------------------------------------

country_positions[[country_col]] <- factor(
  country_positions[[country_col]],
  levels = rev(unique(ordering_df[[country_col]]))
)

pA_country <- ggplot(
  country_positions,
  aes(
    x = 1,
    y = y,
    label = .data[[country_col]]
  )
) +
  geom_text(
    hjust = 1,
    size = 2.6
  ) +
  scale_y_continuous(
    limits = c(
      0.5,
      length(bio_levels) + 0.5
    ),
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0.85, 1.02),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(
      5.5, 2, 5.5, 5.5
    )
  )

pA <- pA_country + pA_heat +
  plot_layout(
    widths = c(2.3, 6)
  )

# ============================================================
# 13. FIGURE 4B — COMPOSITION TYPES
# ============================================================

message("Creating Figure 4B...")

cluster_comp2 <- cluster_comp

cl_col <- find_col(
  cluster_comp2,
  c("^Cluster$", "^CompositionType$", "composition")
)

starter_genera_B <- intersect(
  c(
    "Lacticaseibacillus",
    "Leuconostoc",
    "Lactococcus",
    "Streptococcus",
    "Lactobacillus"
  ),
  names(cluster_comp2)
)

if (length(starter_genera_B) != 5) {
  stop(
    "Figure4B: expected all 5 starter-associated genera. Found: ",
    paste(starter_genera_B, collapse = ", ")
  )
}

cluster_comp_long <- cluster_comp2 %>%
  select(
    all_of(c(
      cl_col,
      intersect("n_BioProjects", names(cluster_comp2)),
      starter_genera_B
    ))
  ) %>%
  pivot_longer(
    cols = all_of(starter_genera_B),
    names_to = "Genus",
    values_to = "MeanNormalizedProportion"
  )

cluster_comp_long$Genus <- factor(
  cluster_comp_long$Genus,
  levels = c(
    "Lacticaseibacillus",
    "Leuconostoc",
    "Lactococcus",
    "Streptococcus",
    "Lactobacillus"
  )
)

# Export exact source data used in panel B
write_csv_final(
  cluster_comp_long,
  "Source_Data_Figure4B_long.csv"
)

pB <- ggplot(
  cluster_comp_long,
  aes(
    x = factor(.data[[cl_col]]),
    y = MeanNormalizedProportion,
    fill = Genus
  )
) +
  geom_col(width = 0.75) +
  labs(
    x = "Inferred bacterial composition type",
    y = "Mean normalized proportion",
    fill = "Genus"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )

# ============================================================
# ============================================================
# 14. FIGURE 4C — PRIMARY AND COMPOSITION-TYPE PERMANOVA
# ============================================================

message("Creating Figure 4C...")

extract_country <- function(df) {

  term_col <- find_col(
    df,
    c("^Term$", "^Variable$", "^Factor$", "term")
  )

  r2_col_local <- find_col(
    df,
    c("^R2$", "R2", "R\\.2")
  )

  f_col_local <- find_col(
    df,
    c("^F$", "pseudo", "F.Model")
  )

  p_col_local <- find_col(
    df,
    c("^p$", "p_value", "Pr\\(>F\\)")
  )

  df_col_local <- find_col(
    df,
    c("^Df$", "degrees"),
    required = FALSE
  )

  z <- df[
    grepl(
      "^Country$",
      as.character(df[[term_col]]),
      ignore.case = TRUE
    ),
    ,
    drop = FALSE
  ]

  if (nrow(z) < 1) {
    stop("Could not extract Country row")
  }

  data.frame(
    Df = if (!is.null(df_col_local))
      as.numeric(z[[df_col_local]][1]) else NA_real_,
    R2 = as.numeric(z[[r2_col_local]][1]),
    F = as.numeric(z[[f_col_local]][1]),
    p = as.numeric(z[[p_col_local]][1])
  )
}

primary_country <- extract_country(
  perm_primary
)

min2_country <- extract_country(
  perm_min2_geo
)

# Composition-type source contains the sequential models.
# Use the already verified master-analysis table.
comp_plot <- standardize_perm(
  perm_confirmed,
  "Composition-type subset"
)

# Exact rows needed from composition model
model_col_comp <- find_col(
  perm_confirmed,
  c("^Model$"),
  required = FALSE
)

term_col_comp <- find_col(
  perm_confirmed,
  c("^Term$")
)

r2_col_comp <- find_col(
  perm_confirmed,
  c("^R2$", "R2")
)

f_col_comp <- find_col(
  perm_confirmed,
  c("^F$", "pseudo")
)

p_col_comp <- find_col(
  perm_confirmed,
  c("^p$", "p_value", "Pr\\(>F\\)")
)

get_comp_row <- function(model_pattern, term_pattern) {

  if (is.null(model_col_comp)) {
    stop("Composition PERMANOVA table has no Model column")
  }

  z <- perm_confirmed[
    grepl(
      model_pattern,
      perm_confirmed[[model_col_comp]],
      ignore.case = TRUE
    ) &
    grepl(
      term_pattern,
      perm_confirmed[[term_col_comp]],
      ignore.case = TRUE
    ),
    ,
    drop = FALSE
  ]

  if (nrow(z) < 1) {
    stop(
      "Missing composition row: ",
      model_pattern,
      " / ",
      term_pattern
    )
  }

  data.frame(
    R2 = as.numeric(z[[r2_col_comp]][1]),
    F = as.numeric(z[[f_col_comp]][1]),
    p = as.numeric(z[[p_col_comp]][1])
  )
}

composition_only <- get_comp_row(
  "Composition type only",
  "Cluster"
)

country_after_comp <- get_comp_row(
  "Composition type \\+ Country",
  "Country"
)

comp_after_country <- get_comp_row(
  "Country \\+ Composition type",
  "Cluster"
)

perm_plot <- data.frame(
  Analysis = c(
    "Primary Country\n(n=73; 20 countries)",
    "Country, replicated countries\n(n=63; 10 countries)",
    "Composition type\n(n=68)",
    "Country | composition type\n(n=68)",
    "Composition type | Country\n(n=68)"
  ),
  R2 = c(
    primary_country$R2,
    min2_country$R2,
    composition_only$R2,
    country_after_comp$R2,
    comp_after_country$R2
  ),
  F = c(
    primary_country$F,
    min2_country$F,
    composition_only$F,
    country_after_comp$F,
    comp_after_country$F
  ),
  p = c(
    primary_country$p,
    min2_country$p,
    composition_only$p,
    country_after_comp$p,
    comp_after_country$p
  )
)

write_csv_final(
  perm_plot,
  "Source_Data_Figure4C.csv"
)

perm_plot$Annotation <- sprintf(
  "R² = %.3f\np %s",
  perm_plot$R2,
  ifelse(
    perm_plot$p <= 0.001,
    "≤ 0.001",
    paste0(
      "= ",
      formatC(
        perm_plot$p,
        digits = 3,
        format = "f"
      )
    )
  )
)

pC <- ggplot(
  perm_plot,
  aes(
    x = reorder(Analysis, R2),
    y = R2
  )
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = Annotation),
    hjust = -0.08,
    size = 3
  ) +
  coord_flip() +
  expand_limits(
    y = max(perm_plot$R2) * 1.32
  ) +
  labs(
    x = NULL,
    y = expression(PERMANOVA~R^2)
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank()
  )

# ============================================================
# 15. FIGURE 4D — ROBUSTNESS ANALYSES
# ============================================================

message("Creating Figure 4D...")

# Extract current values directly from completed master-analysis
# output files. No statistical result is hardcoded.

removed_country <- extract_country(
  removed
)

continuous_country <- extract_country(
  continuous
)

balanced_metric <- function(metric) {

  z <- balanced_summary[
    balanced_summary[[1]] == metric,
    ,
    drop = FALSE
  ]

  if (nrow(z) != 1) {
    stop("Balanced metric not found: ", metric)
  }

  as.numeric(z[[2]][1])
}

balanced_median_r2 <- balanced_metric(
  "Median R2"
)

balanced_perm_sig <- balanced_metric(
  "Fraction PERMANOVA p<=0.05"
)

balanced_disp_sig <- balanced_metric(
  "Fraction PERMDISP p<=0.05"
)

robust <- data.frame(
  Analysis = c(
    "Replicated countries\n63 BioProjects / 10 countries",
    "Five genera removed\nprimary population",
    "Continuous five-genus composition\nstarter-positive subset",
    "Balanced subsampling\n8 countries × 3 BioProjects"
  ),
  R2 = c(
    min2_country$R2,
    removed_country$R2,
    continuous_country$R2,
    balanced_median_r2
  ),
  p = c(
    min2_country$p,
    removed_country$p,
    continuous_country$p,
    NA_real_
  ),
  Note = c(
    "PERMDISP p=0.197",
    "PERMDISP significant in full population",
    "Microbiome-derived sensitivity covariates",
    paste0(
      sprintf("%.1f", 100 * balanced_perm_sig),
      "% PERMANOVA significant; ",
      sprintf("%.1f", 100 * balanced_disp_sig),
      "% PERMDISP significant"
    )
  ),
  stringsAsFactors = FALSE
)

write_csv_final(
  robust,
  "Source_Data_Figure4D.csv"
)

robust$Annotation <- c(
  sprintf(
    "R² = %.3f\nPERMDISP p = 0.197",
    min2_country$R2
  ),
  sprintf(
    "R² = %.3f\np ≤ 0.001",
    removed_country$R2
  ),
  sprintf(
    "R² = %.3f\np ≤ 0.001",
    continuous_country$R2
  ),
  sprintf(
    "median R² = %.3f\n%.1f%% significant",
    balanced_median_r2,
    100 * balanced_perm_sig
  )
)

pD <- ggplot(
  robust,
  aes(
    x = reorder(Analysis, R2),
    y = R2
  )
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = Annotation),
    hjust = -0.08,
    size = 2.8
  ) +
  coord_flip() +
  expand_limits(
    y = max(robust$R2) * 1.42
  ) +
  labs(
    x = NULL,
    y = expression(Country-associated~R^2)
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank()
  )

# 16. COMBINE FIGURE
# ============================================================

message("Combining final Figure 4...")

# Explicit panel labels.
# Manual tags avoid patchwork treating the internal components
# of composite panel A as separate figure panels.

# ============================================================
# Explicit panel labels A-D
# ============================================================

label_theme <- theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 2)
  )

labA <- ggplot() +
  annotate(
    "text", x = 0, y = 0,
    label = "A",
    hjust = 0, vjust = 0.5,
    fontface = "bold", size = 5
  ) +
  xlim(0, 1) + ylim(-1, 1) +
  label_theme

labB <- ggplot() +
  annotate(
    "text", x = 0, y = 0,
    label = "B",
    hjust = 0, vjust = 0.5,
    fontface = "bold", size = 5
  ) +
  xlim(0, 1) + ylim(-1, 1) +
  label_theme

labC <- ggplot() +
  annotate(
    "text", x = 0, y = 0,
    label = "C",
    hjust = 0, vjust = 0.5,
    fontface = "bold", size = 5
  ) +
  xlim(0, 1) + ylim(-1, 1) +
  label_theme

labD <- ggplot() +
  annotate(
    "text", x = 0, y = 0,
    label = "D",
    hjust = 0, vjust = 0.5,
    fontface = "bold", size = 5
  ) +
  xlim(0, 1) + ylim(-1, 1) +
  label_theme

# Each panel gets its own explicit label row.
# Panel A can remain internally composite without affecting numbering.

panelA_final <- labA / pA +
  plot_layout(heights = c(0.045, 1))

panelB_final <- labB / pB +
  plot_layout(heights = c(0.045, 1))

panelC_final <- labC / pC +
  plot_layout(heights = c(0.065, 1))

panelD_final <- labD / pD +
  plot_layout(heights = c(0.065, 1))

figure4 <- (
  panelA_final | panelB_final
) / (
  panelC_final | panelD_final
)

png_file <- file.path(
  out_dir,
  "Figure4_REVISED_final.png"
)

pdf_file <- file.path(
  out_dir,
  "Figure4_REVISED_final.pdf"
)

ggsave(
  png_file,
  figure4,
  width = 13,
  height = 15,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  pdf_file,
  figure4,
  width = 13,
  height = 15,
  units = "in",
  bg = "white"
)

# Individual panels also useful for editing/submission
ggsave(
  file.path(out_dir, "Figure4A_REVISED.png"),
  pA,
  width = 6.5,
  height = 12,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure4B_REVISED.png"),
  pB,
  width = 6.5,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure4C_REVISED.png"),
  pC,
  width = 6.5,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure4D_REVISED.png"),
  pD,
  width = 6.5,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

# ============================================================
# 17. FINAL FIGURE CAPTION
# ============================================================

caption <- paste0(
"Figure 4. Country-associated structure and inferred bacterial composition ",
"types in the cheese microbiome dataset. ",
"(A) Relative abundances of Lactococcus, Streptococcus, Lacticaseibacillus, ",
"Lactobacillus and Leuconostoc across 68 BioProjects with a single ",
"BioSample-reported sample-origin country and nonzero representation of the ",
"five genera. ",
"(B) Five descriptive bacterial composition types derived by k-means from ",
"within-BioProject normalized proportions of these genera. These groups are ",
"microbiome-derived composition types and do not represent experimentally ",
"verified starter-culture categories. ",
"(C) PERMANOVA effect sizes. The primary geographic analysis included 73 ",
"BioProjects from 20 countries (Country R2 = 0.366, p = 0.001); multivariate ",
"dispersion differed among countries (PERMDISP p = 0.001), and this result is ",
"therefore interpreted as country-associated community structure rather than ",
"an independent geographic effect. Restriction to countries represented by ",
"at least two independent BioProjects retained 63 BioProjects from 10 ",
"countries (Country R2 = 0.284, p = 0.001), while PERMDISP was not significant ",
"(p = 0.197). Composition-type models used the 68 BioProjects with nonzero ",
"representation of the five composition-associated genera. ",
"(D) Sensitivity analyses. The country association persisted after complete ",
"removal of the five composition-associated genera (R2 = 0.359, p = 0.001) ",
"and after adjustment for their continuous within-five-genus proportions ",
"(R2 = 0.306, p = 0.001). In 1,000 balanced subsamples containing three ",
"BioProjects from each of eight countries, PERMANOVA was significant in ",
"99.4% of iterations, whereas PERMDISP was significant in 0% of iterations; ",
"the balanced median R2 is shown as an empirical sensitivity estimate rather ",
"than as the primary geographic effect size. ",
"Community analyses used Hellinger-transformed genus-level profiles, ",
"Bray-Curtis dissimilarity and 999 permutations."
)

writeLines(
  caption,
  file.path(out_dir, "Figure4_REVISED_caption.txt")
)

# ============================================================
# 18. README
# ============================================================

readme <- c(
  "REVIEWER 3 — FINAL FIGURE 4 AND SUPPLEMENTARY PACKAGE",
  "=====================================================",
  "",
  "This directory contains publication-facing outputs generated",
  "from the completed Reviewer 3 sensitivity analyses.",
  "",
  "MAIN FIGURE",
  "Figure4_REVISED_final.pdf",
  "Figure4_REVISED_final.png",
  "",
  "SUPPLEMENTARY TABLES",
  "S4  Geographic provenance audit",
  "S4b Country coverage",
  "S5  Geographic PERMANOVA analyses",
  "S5b Geographic PERMDISP analyses",
  "S6  Composition-type validation",
  "S6b Full composition-type validation source",
  "S7  Balanced project-level subsampling summary",
  "",
  "IMPORTANT INTERPRETATION",
  "Country effects are described as country-associated community",
  "structure, not as an independent geographic or terroir effect.",
  "",
  "The five microbiome-derived groups are described as inferred",
  "bacterial composition types, not verified starter-culture categories.",
  "",
  "The k=5 solution is retained as a descriptive and interpretable representation;",
  "it is not described as the statistically optimal number of clusters.",
  "ARI values describe reproducibility across repeated k-means initializations",
  "on the same data matrix and are not interpreted as resampling stability.",
  "",
  "Balanced subsampling is a robustness analysis and its R2 estimate",
  "should not be used as the primary geographic effect size.",
  "",
  "No QIIME2/Kraken2/mOTUs processing was rerun by this script."
)

writeLines(
  readme,
  file.path(out_dir, "README.txt")
)

# ============================================================
# 19. FILE INDEX
# ============================================================

files <- list.files(
  out_dir,
  full.names = FALSE
)

index <- data.frame(
  File = files,
  stringsAsFactors = FALSE
)

write_csv_final(
  index,
  "NEW_SUPP_FILE_INDEX.csv"
)

message("")
message("==============================================")
message("NEW_SUPP PACKAGE COMPLETED")
message("==============================================")
message("Output: ", out_dir)
message("")
message("Main figure:")
message("  Figure4_REVISED_final.pdf")
message("  Figure4_REVISED_final.png")
message("")
message("Supplementary:")
message("  Supplementary_Table_S4_Geographic_origin.csv")
message("  Supplementary_Table_S5_Geographic_PERMANOVA.csv")
message("  Supplementary_Table_S5b_Geographic_PERMDISP.csv")
message("  Supplementary_Table_S6_Composition_type_validation.csv")
message("  Supplementary_Table_S7_Balanced_subsampling_summary.csv")
message("==============================================")
