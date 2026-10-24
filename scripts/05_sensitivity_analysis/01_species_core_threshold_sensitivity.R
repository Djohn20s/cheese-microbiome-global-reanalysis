#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

# ============================================================
# PATHS
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
audit_dir <- file.path(
  project_dir,
  "data/08_reporting/bioproject_level_audit_amplicon_wgs_v8"
)

out_dir <- file.path(
  project_dir,
  "data/08_reporting/Rew_3/Comment5_core_sensitivity"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

kraken_file <- file.path(
  audit_dir,
  "Kraken2_WGS_species_BioProject_taxa_wide.csv"
)

motus_file <- file.path(
  audit_dir,
  "mOTUs_WGS_species_BioProject_taxa_wide.csv"
)

# ============================================================
# PARAMETERS
# ============================================================

# Current detection threshold used in Figure 5
primary_detection_threshold <- 0.0001

# Sensitivity around detection threshold
detection_thresholds <- c(
  0,
  0.00001,
  0.0001,
  0.001
)

# Prevalence sensitivity
prevalence_thresholds <- c(
  0.10,
  0.25,
  0.50,
  0.75,
  0.90
)

# Existing mean-abundance criterion
mean_abundance_threshold <- 0.0001

# ============================================================
# LOAD
# ============================================================

K <- fread(kraken_file)
M <- fread(motus_file)

meta_candidates <- c(
  "BioProject",
  "CenterName_standardized",
  "CenterName",
  "Continent",
  "Country",
  "Platform"
)

get_tax_cols <- function(dt) {
  setdiff(names(dt), intersect(names(dt), meta_candidates))
}

# ============================================================
# COLLAPSE DUPLICATE BIOPROJECT ROWS
#
# Important because matrices contain:
# Kraken2: 21 rows / 20 unique BioProjects
# mOTUs:   20 rows / 19 unique BioProjects
#
# Taxonomic profiles are averaged within duplicated BioProjects.
# ============================================================

collapse_bp <- function(dt) {

  tax_cols <- get_tax_cols(dt)

  dt2 <- dt[
    ,
    lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE)),
    by = BioProject,
    .SDcols = tax_cols
  ]

  dt2
}

K <- collapse_bp(K)
M <- collapse_bp(M)

cat("Kraken2 unique BioProjects:", nrow(K), "\n")
cat("mOTUs unique BioProjects:", nrow(M), "\n")

# ============================================================
# COMMON SPECIES
# ============================================================

k_species <- setdiff(names(K), "BioProject")
m_species <- setdiff(names(M), "BioProject")

common_species <- intersect(k_species, m_species)

cat("Common species:", length(common_species), "\n")

# ============================================================
# BASIC SUMMARY FUNCTION
# ============================================================

summarise_method <- function(dt, species, detection_threshold) {

  n_bp <- nrow(dt)

  mat <- as.matrix(dt[, ..species])
  storage.mode(mat) <- "numeric"

  prevalence <- colMeans(
    mat > detection_threshold,
    na.rm = TRUE
  )

  mean_abundance <- colMeans(
    mat,
    na.rm = TRUE
  )

  detected_n <- colSums(
    mat > detection_threshold,
    na.rm = TRUE
  )

  data.table(
    Species = species,
    denominator = n_bp,
    detected_n = detected_n,
    prevalence = prevalence,
    mean_abundance = mean_abundance
  )
}

# ============================================================
# PRIMARY SPECIES TABLE
# Both methods evaluated independently
# ============================================================

K_primary <- summarise_method(
  K,
  common_species,
  primary_detection_threshold
)

M_primary <- summarise_method(
  M,
  common_species,
  primary_detection_threshold
)

setnames(
  K_primary,
  c("denominator", "detected_n", "prevalence", "mean_abundance"),
  c(
    "Kraken2_denominator",
    "Kraken2_detected_n",
    "Kraken2_prevalence",
    "Kraken2_mean_abundance"
  )
)

setnames(
  M_primary,
  c("denominator", "detected_n", "prevalence", "mean_abundance"),
  c(
    "mOTUs_denominator",
    "mOTUs_detected_n",
    "mOTUs_prevalence",
    "mOTUs_mean_abundance"
  )
)

primary <- merge(
  K_primary,
  M_primary,
  by = "Species",
  all = FALSE
)

primary[, prevalence_min :=
          pmin(Kraken2_prevalence, mOTUs_prevalence)]

primary[, mean_abundance_min :=
          pmin(Kraken2_mean_abundance, mOTUs_mean_abundance)]

# ============================================================
# PREVALENCE SENSITIVITY
#
# IMPORTANT:
# A species qualifies only when BOTH pipelines independently
# satisfy BOTH prevalence and abundance thresholds.
# ============================================================

prev_sensitivity <- rbindlist(
  lapply(prevalence_thresholds, function(p_thr) {

    tmp <- copy(primary)

    tmp[, qualifies :=
          Kraken2_prevalence >= p_thr &
          mOTUs_prevalence >= p_thr &
          Kraken2_mean_abundance >= mean_abundance_threshold &
          mOTUs_mean_abundance >= mean_abundance_threshold]

    data.table(
      prevalence_threshold = p_thr,
      mean_abundance_threshold = mean_abundance_threshold,
      n_common_species = nrow(tmp),
      n_species_meeting_both_methods = sum(tmp$qualifies)
    )
  })
)

fwrite(
  prev_sensitivity,
  file.path(
    out_dir,
    "Core_species_prevalence_sensitivity.tsv"
  ),
  sep = "\t"
)

# ============================================================
# MEMBERSHIP AT EACH PREVALENCE THRESHOLD
# ============================================================

membership <- rbindlist(
  lapply(prevalence_thresholds, function(p_thr) {

    tmp <- primary[
      Kraken2_prevalence >= p_thr &
      mOTUs_prevalence >= p_thr &
      Kraken2_mean_abundance >= mean_abundance_threshold &
      mOTUs_mean_abundance >= mean_abundance_threshold
    ]

    if (nrow(tmp) == 0) return(NULL)

    tmp[, prevalence_threshold := p_thr]

    tmp
  }),
  fill = TRUE
)

fwrite(
  membership,
  file.path(
    out_dir,
    "Core_species_membership_by_prevalence.tsv"
  ),
  sep = "\t"
)

# ============================================================
# DETECTION-THRESHOLD SENSITIVITY
#
# For every detection threshold x prevalence threshold.
# ============================================================

detection_sensitivity <- rbindlist(
  lapply(detection_thresholds, function(det_thr) {

    ks <- summarise_method(K, common_species, det_thr)
    ms <- summarise_method(M, common_species, det_thr)

    setnames(
      ks,
      c("denominator", "detected_n", "prevalence", "mean_abundance"),
      c("K_den", "K_n", "K_prev", "K_mean")
    )

    setnames(
      ms,
      c("denominator", "detected_n", "prevalence", "mean_abundance"),
      c("M_den", "M_n", "M_prev", "M_mean")
    )

    z <- merge(ks, ms, by = "Species")

    rbindlist(
      lapply(prevalence_thresholds, function(p_thr) {

        qualifies <-
          z$K_prev >= p_thr &
          z$M_prev >= p_thr &
          z$K_mean >= mean_abundance_threshold &
          z$M_mean >= mean_abundance_threshold

        data.table(
          detection_threshold = det_thr,
          prevalence_threshold = p_thr,
          mean_abundance_threshold = mean_abundance_threshold,
          n_species_meeting_both_methods = sum(qualifies)
        )
      })
    )
  })
)

fwrite(
  detection_sensitivity,
  file.path(
    out_dir,
    "Core_species_detection_prevalence_sensitivity.tsv"
  ),
  sep = "\t"
)

# ============================================================
# LEAVE-ONE-BIOPROJECT-OUT
#
# Main robustness definition:
# prevalence >=50% independently in both methods
# abundance >=0.0001 independently in both methods
#
# Each BioProject is removed in turn from whichever matrix
# contains it.
# ============================================================

loo_prevalence <- 0.50

all_bp <- sort(unique(c(K$BioProject, M$BioProject)))

loo_results <- list()
loo_membership <- list()

for (bp_remove in all_bp) {

  K_sub <- K[BioProject != bp_remove]
  M_sub <- M[BioProject != bp_remove]

  ks <- summarise_method(
    K_sub,
    common_species,
    primary_detection_threshold
  )

  ms <- summarise_method(
    M_sub,
    common_species,
    primary_detection_threshold
  )

  setnames(
    ks,
    c("denominator", "detected_n", "prevalence", "mean_abundance"),
    c("K_den", "K_n", "K_prev", "K_mean")
  )

  setnames(
    ms,
    c("denominator", "detected_n", "prevalence", "mean_abundance"),
    c("M_den", "M_n", "M_prev", "M_mean")
  )

  z <- merge(ks, ms, by = "Species")

  z[, qualifies :=
      K_prev >= loo_prevalence &
      M_prev >= loo_prevalence &
      K_mean >= mean_abundance_threshold &
      M_mean >= mean_abundance_threshold]

  loo_results[[bp_remove]] <- data.table(
    excluded_BioProject = bp_remove,
    Kraken2_denominator = nrow(K_sub),
    mOTUs_denominator = nrow(M_sub),
    n_widespread_species = sum(z$qualifies)
  )

  loo_membership[[bp_remove]] <- z[
    qualifies == TRUE,
    .(
      excluded_BioProject = bp_remove,
      Species
    )
  ]
}

loo_results <- rbindlist(loo_results)
loo_membership <- rbindlist(loo_membership)

fwrite(
  loo_results,
  file.path(
    out_dir,
    "Widespread_core_leave_one_BioProject_out.tsv"
  ),
  sep = "\t"
)

fwrite(
  loo_membership,
  file.path(
    out_dir,
    "Widespread_core_leave_one_BioProject_out_membership.tsv"
  ),
  sep = "\t"
)

# ============================================================
# STABILITY OF EACH SPECIES ACROSS LOO ITERATIONS
# ============================================================

if (nrow(loo_membership) > 0) {

  loo_stability <- loo_membership[
    ,
    .(
      n_iterations_core = .N
    ),
    by = Species
  ]

  loo_stability[, n_iterations_total := length(all_bp)]

  loo_stability[
    ,
    fraction_iterations_core :=
      n_iterations_core / n_iterations_total
  ]

  setorder(
    loo_stability,
    -fraction_iterations_core,
    Species
  )

  fwrite(
    loo_stability,
    file.path(
      out_dir,
      "Widespread_core_species_LOO_stability.tsv"
    ),
    sep = "\t"
  )
}

# ============================================================
# FINAL PRIMARY CATEGORIES
# ============================================================

primary[
  ,
  category := fifelse(
    Kraken2_prevalence >= 0.75 &
      mOTUs_prevalence >= 0.75 &
      Kraken2_mean_abundance >= mean_abundance_threshold &
      mOTUs_mean_abundance >= mean_abundance_threshold,
    "Highly widespread (>=75% in both)",
    fifelse(
      Kraken2_prevalence >= 0.50 &
        mOTUs_prevalence >= 0.50 &
        Kraken2_mean_abundance >= mean_abundance_threshold &
        mOTUs_mean_abundance >= mean_abundance_threshold,
      "Widespread (>=50% in both)",
      fifelse(
        Kraken2_prevalence >= 0.10 &
          mOTUs_prevalence >= 0.10 &
          Kraken2_mean_abundance >= mean_abundance_threshold &
          mOTUs_mean_abundance >= mean_abundance_threshold,
        "Recurrent (10-<50% in at least one method)",
        "Not retained"
      )
    )
  )
]

setorder(
  primary,
  -prevalence_min,
  -mean_abundance_min
)

fwrite(
  primary,
  file.path(
    out_dir,
    "Species_primary_classification.tsv"
  ),
  sep = "\t"
)

# ============================================================
# SUMMARY
# ============================================================

cat("\n============================================\n")
cat("REVIEWER 3 COMMENT 5 — CORE SENSITIVITY\n")
cat("============================================\n\n")

cat("Kraken2 denominator:", nrow(K), "BioProjects\n")
cat("mOTUs denominator:", nrow(M), "BioProjects\n")
cat("Common species:", length(common_species), "\n\n")

cat("Primary detection threshold:",
    primary_detection_threshold, "\n")
cat("Mean abundance threshold:",
    mean_abundance_threshold, "\n\n")

cat("Species meeting BOTH methods independently:\n")
print(prev_sensitivity)

cat("\nPrimary classification:\n")
print(primary[, .N, by = category][order(category)])

cat("\nLeave-one-BioProject-out at >=50% prevalence:\n")
cat(
  "range =",
  min(loo_results$n_widespread_species),
  "-",
  max(loo_results$n_widespread_species),
  "\n"
)

cat("\nOutputs:\n", out_dir, "\n")