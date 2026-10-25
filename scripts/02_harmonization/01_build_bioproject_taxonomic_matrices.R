#!/usr/bin/env Rscript

# ============================================================
# BioProject-level statistical audit: QIIME2 Amplicon + Kraken2 WGS + mOTUs WGS
# Version: v8
# - output in new folder under data/08_reporting
# - optional Kraken2 matrix construction from long TSV files
# - fixed duplicated sample columns such as SRR_x / SRR_y / SRR
# - curated institution/country/continent metadata support
# - BioProject is analytical unit
# - v7: exports SRA and BioProject counts per Country and per Continent
# - v8: uses final stage5 metadata as analytical input and exports publication-ready filtering/analysis tables
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(vegan)
  library(ggrepel)
  library(parallel)
  library(stringi)
})

# ============================================================
# SETTINGS
# ============================================================

n_cores <- 28
setDTthreads(n_cores)
options(stringsAsFactors = FALSE)
options(scipen = 999)

# TRUE only for first run, when Kraken2 matrices do not exist yet.
# After successful first run, set FALSE.
build_kraken_matrix <- FALSE

# Repository root (portable; independent of current working directory)
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)

if (length(file_arg) == 1) {
  script_file <- normalizePath(sub("^--file=", "", file_arg))
  repo_dir <- normalizePath(file.path(dirname(script_file), "..", ".."))
} else {
  repo_dir <- normalizePath(getwd())
}

# Source project containing the large SRA-level input matrices.
# These files are intentionally not duplicated in the GitHub repository.
project_dir <- normalizePath(file.path(repo_dir, ".."))

# Publication-ready BioProject-level outputs are written to the repository.
out_dir <- file.path(
  repo_dir,
  "processed_data/bioproject_matrices"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# HOST / CONTAMINANT TAXA TO REMOVE
# ============================================================

excluded_host_taxa <- c(
  "Homo",
  "Mus",
  "Rattus",
  "Bos",
  "Sus",
  "Gallus",
  "Canis",
  "Felis"
)

# ============================================================
# INPUT PATHS
# ============================================================

# Initial metadata table used only for independent count checks of the initial NCBI query.
initial_metadata_path <- file.path(
  project_dir,
  "data/03_primary/ncbi/cheese_metagenome/v2/runinfo_with_continent.csv"
)

# Final curated stage5 metadata used as the analytical metadata input.
stage5_amplicon_metadata_path <- file.path(
  project_dir,
  "data/03_primary/ncbi/cheese_metagenome/v2/stage5/runinfo_AMPLICON.csv"
)

stage5_wgs_metadata_path <- file.path(
  project_dir,
  "data/03_primary/ncbi/cheese_metagenome/v2/stage5/runinfo_WGS.csv"
)

metadata_path <- paste(
  stage5_amplicon_metadata_path,
  stage5_wgs_metadata_path,
  sep = " + "
)

# Optional curated map with columns:
# CenterName, CenterName_standardized, Country, Continent
center_country_map_path <- file.path(
  project_dir,
  "data/03_primary/ncbi/cheese_metagenome/v2/center_country_map_curated.csv"
)

qiime2_amplicon_matrix_path <- file.path(
  project_dir,
  "data/08_reporting/sra_profile_amplicon/matrices/tax_matrix_1-2marge_curated_list.tsv"
)

motus_matrix_path <- file.path(
  project_dir,
  "data/08_reporting/sra_profile_wgs/matrices/motus_matrix.tsv"
)

kraken_taxa_long_dir <- file.path(
  project_dir,
  "data/08_reporting/sra_profile_wgs/matrices/kraken_taxa_long"
)

kraken_genus_matrix_path <- file.path(
  project_dir,
  "data/08_reporting/sra_profile_wgs/matrices/kraken_genus_matrix.tsv"
)

kraken_species_matrix_path <- file.path(
  project_dir,
  "data/08_reporting/sra_profile_wgs/matrices/kraken_species_matrix.tsv"
)

# ============================================================
# THRESHOLDS
# ============================================================

min_bioprojects_for_core <- 0.50
min_mean_abundance_core <- 0.001
species_detection_bioprojects <- 3
species_detection_mean_abundance <- 0.0001
top_n_taxa <- 30

# ============================================================
# PREVALENCE DETECTION THRESHOLDS (relative abundance)
# Change only this value for the whole project:
detection_threshold <- 0.0001   # 0.01%

# Audit thresholds automatically reported
prevalence_audit_thresholds <- c(
  gt0   = 0,
  `0.01%` = 0.0001,
  `0.1%`  = 0.001,
  `1%`    = 0.01
)

# Filters used only for cross-pipeline concordance plots/statistics.
# Raw concordance tables are still saved before filtering.
min_prevalence_for_correlation <- 0.10
min_mean_abundance_for_correlation <- 0.0001

starter_genera <- c(
  "Lactococcus",
  "Streptococcus",
  "Lactobacillus"
)

nonstarter_genera <- c(
  "Leuconostoc",
  "Brevibacterium",
  "Corynebacterium"
)

# ============================================================
# README
# ============================================================

readme_text <- c(
  "# BioProject-level audit of public cheese microbiome datasets",
  "",
  "Version: v8",
  "",
  "This folder contains BioProject-level analyses for QIIME2 Amplicon, Kraken2 WGS, and mOTUs WGS datasets.",
  "",
  "## Main analyses",
  "",
  "1. Aggregation of SRA-level taxonomic profiles to BioProject-level mean relative abundance.",
  "2. Prevalence of bacterial genera across BioProjects.",
  "3. Genus-level concordance among QIIME2, Kraken2, and mOTUs.",
  "4. WGS species richness within dominant genera.",
  "5. Distribution of common starter and non-starter genera across BioProjects.",
  "6. BioProject-level Bray-Curtis PCoA and PERMANOVA.",
  "7. Dataset audit plots: SRA counts per BioProject.",
  "8. Curated CenterName/Country/Continent metadata integration.",
  "9. Country-level and continent-level SRA/BioProject summaries.",
  "10. Publication-ready filtering workflow and final analysis summary tables.",
  "",
  "## Detection thresholds",
  "",
  paste0("- Core genus: detected in >= ", min_bioprojects_for_core * 100, "% BioProjects and mean abundance >= ", min_mean_abundance_core),
  paste0("- WGS species retained: detected in >= ", species_detection_bioprojects, " BioProjects and mean abundance >= ", species_detection_mean_abundance),
  "",
  "## Notes",
  "",
  "BioProject is treated as the independent analytical unit to reduce pseudoreplication caused by multiple SRA runs within the same public project.",
  "Geographic metadata are based on the submitting institution unless true cheese origin metadata are available.",
  "If center_country_map_curated.csv exists, curated country and continent assignments are used.",
  "Kraken2 matrices can be built from long-format TSV files using build_kraken_matrix <- TRUE."
)

writeLines(readme_text, file.path(out_dir, "README.md"))

# ============================================================
# BASIC HELPERS
# ============================================================

normalize_sample_names <- function(x) {
  x <- as.character(x)
  x <- gsub("(_x|_y)$", "", x)
  x
}

normalize_center_name <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- stringr::str_replace_all(x, '"', "")
  x <- stringr::str_replace_all(x, "'", "")
  x <- stringr::str_replace_all(x, "\\.", "")
  x <- stringr::str_replace_all(x, ",", "")
  x <- stringr::str_replace_all(x, "\\(", "")
  x <- stringr::str_replace_all(x, "\\)", "")
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- stringr::str_to_upper(x)
  x <- stringr::str_squish(x)
  x
}

clean_tax <- function(x) {
  x <- as.character(x)
  x <- gsub("^.*__", "", x)
  x <- gsub("_", " ", x)
  x <- trimws(x)
  x[x == "" | is.na(x)] <- NA_character_
  x
}

extract_rank <- function(tax, rank = "genus") {

  tax <- as.character(tax)

  out <- vapply(tax, function(x) {

    if (is.na(x) || x == "") return(NA_character_)

    parts <- unlist(strsplit(x, ";", fixed = TRUE))
    parts <- trimws(parts)

    if (rank == "genus") {
      hit <- parts[grepl("^g__", parts)]
      if (length(hit) > 0) {
        y <- gsub("^g__", "", hit[length(hit)])
        y <- gsub("_", " ", y)
        y <- trimws(y)
        if (!is.na(y) && y != "" && y != "__") return(y)
      }

      sp <- parts[length(parts)]
      sp <- gsub("^s__", "", sp)
      sp <- gsub("_", " ", sp)
      sp <- trimws(sp)
      if (grepl(" ", sp)) return(strsplit(sp, " ")[[1]][1])
      return(NA_character_)
    }

    if (rank == "species") {
      hit <- parts[grepl("^s__", parts)]
      if (length(hit) > 0) {
        y <- gsub("^s__", "", hit[length(hit)])
        y <- gsub("_", " ", y)
        y <- trimws(y)
        if (!is.na(y) && y != "" && y != "__") return(y)
      }

      sp <- parts[length(parts)]
      sp <- gsub("^s__", "", sp)
      sp <- gsub("_", " ", sp)
      sp <- trimws(sp)
      if (grepl(" ", sp)) return(sp)
      return(NA_character_)
    }

    NA_character_

  }, character(1))

  out
}

# ============================================================
# MATRIX LOADERS
# ============================================================

read_tax_matrix <- function(path) {

  message("Reading matrix: ", path)

  if (!file.exists(path)) {
    stop("Matrix file does not exist: ", path)
  }

  dt <- fread(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    data.table = TRUE,
    showProgress = TRUE
  )

  setnames(dt, names(dt)[1], "Tax")
  dt[, Tax := as.character(Tax)]

  dt <- dt[
    !grepl("TOTAL_READS_RETAINED", Tax) &
      !grepl("^Unassigned", Tax, ignore.case = TRUE) &
      !grepl("^unclassified$", Tax, ignore.case = TRUE) &
      !grepl("^root$", Tax, ignore.case = TRUE) &
      !grepl("^d__;p__;c__;o__;f__;g__;s__$", Tax)
  ]

  old_cols <- setdiff(names(dt), "Tax")
  new_cols <- normalize_sample_names(old_cols)

  dt[, (old_cols) := lapply(.SD, as.numeric), .SDcols = old_cols]

  if (any(duplicated(new_cols))) {

    message(
      "Collapsing duplicated sample columns after normalization: ",
      sum(duplicated(new_cols))
    )

    mat <- as.matrix(dt[, ..old_cols])
    colnames(mat) <- new_cols
    rownames(mat) <- dt$Tax
    mat[is.na(mat)] <- 0

    collapsed <- rowsum(
      t(mat),
      group = colnames(mat),
      reorder = FALSE
    )
    collapsed <- t(collapsed)

    out <- as.data.table(collapsed, keep.rownames = "Tax")
    return(out)
  }

  setnames(dt, old_cols, new_cols)
  dt
}

to_relative_abundance <- function(dt) {

  sample_cols <- setdiff(names(dt), "Tax")

  mat <- as.matrix(dt[, ..sample_cols])
  rownames(mat) <- dt$Tax
  mat[is.na(mat)] <- 0

  col_sums <- colSums(mat, na.rm = TRUE)

  message("Samples before abundance filtering: ", ncol(mat))
  message("Samples with sum > 0: ", sum(col_sums > 0))
  message("Samples with sum == 0: ", sum(col_sums == 0))

  if (sum(col_sums > 0) == 0) {
    stop("All samples have zero abundance after loading matrix.")
  }

  mat <- mat[, col_sums > 0, drop = FALSE]
  col_sums <- colSums(mat, na.rm = TRUE)

  if (max(col_sums, na.rm = TRUE) > 2) {
    mat <- sweep(mat, 2, col_sums, "/")
  }

  mat[is.na(mat)] <- 0

  as.data.table(mat, keep.rownames = "Tax")
}

collapse_to_rank <- function(dt, rank = "genus") {

  message("Collapsing to rank: ", rank)

  dt[, RankTax := extract_rank(Tax, rank = rank)]
  dt <- dt[!is.na(RankTax) & RankTax != "" & RankTax != "Unclassified"]

  # Remove host / contaminant taxa at collapsed rank level.
  dt <- dt[!RankTax %in% excluded_host_taxa]

  sample_cols <- setdiff(names(dt), c("Tax", "RankTax"))

  collapsed <- dt[
    ,
    lapply(.SD, sum, na.rm = TRUE),
    by = RankTax,
    .SDcols = sample_cols
  ]

  setnames(collapsed, "RankTax", "Tax")
  collapsed
}

# ============================================================
# KRAKEN2 MATRIX BUILDER
# ============================================================

read_single_kraken_file <- function(file, rank_level = "genus") {

  dt <- tryCatch(
    fread(file, sep = "\t", header = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )

  if (is.null(dt)) {
    warning("Skipping unreadable file: ", file)
    return(NULL)
  }

  if (!all(c("Run", "percent", "rank", "name") %in% names(dt))) {
    warning("Skipping malformed file: ", file)
    return(NULL)
  }

  if (rank_level == "genus") {
    dt <- dt[rank == "G"]
  } else if (rank_level == "species") {
    dt <- dt[rank == "S"]
  } else {
    stop("Unsupported Kraken rank_level: ", rank_level)
  }

  if (nrow(dt) == 0) return(NULL)

  dt <- dt[
    !grepl("^unclassified", name, ignore.case = TRUE) &
      !name %in% c("root", "cellular organisms")
  ]

  if (nrow(dt) == 0) return(NULL)

  run_id <- unique(dt$Run)[1]

  out <- dt[
    ,
    .(
      Abundance = sum(as.numeric(percent), na.rm = TRUE) / 100
    ),
    by = .(Tax = name)
  ]

  out[, Run := run_id]
  out[, .(Run, Tax, Abundance)]
}

build_kraken_matrix_from_long <- function(
    kraken_dir,
    out_path,
    rank_level = "genus",
    n_cores = 28
) {

  message("Building Kraken2 ", rank_level, " matrix from: ", kraken_dir)

  files <- list.files(
    kraken_dir,
    pattern = "\\.tsv$",
    full.names = TRUE
  )

  message("Kraken2 files found: ", length(files))

  if (length(files) == 0) {
    stop("No Kraken2 .tsv files found in: ", kraken_dir)
  }

  res <- mclapply(
    files,
    read_single_kraken_file,
    rank_level = rank_level,
    mc.cores = n_cores
  )

  long <- rbindlist(res, use.names = TRUE, fill = TRUE)

  if (nrow(long) == 0) {
    stop("No Kraken2 records retained for rank: ", rank_level)
  }

  message("Rows retained for Kraken2 ", rank_level, ": ", nrow(long))

  mat <- dcast(
    long,
    Tax ~ Run,
    value.var = "Abundance",
    fill = 0
  )

  setorder(mat, Tax)
  fwrite(mat, out_path, sep = "\t")

  message("Saved Kraken2 ", rank_level, " matrix: ", out_path)

  invisible(mat)
}

# ============================================================
# BIOPROJECT AGGREGATION AND ANALYSIS
# ============================================================

aggregate_to_bioproject <- function(dt, metadata, method_name) {

  sample_cols <- setdiff(names(dt), "Tax")
  common_samples <- intersect(sample_cols, metadata$SampleID)

  message(method_name, ": matrix samples = ", length(sample_cols))
  message(method_name, ": metadata samples = ", uniqueN(metadata$SampleID))
  message(method_name, ": common samples = ", length(common_samples))

  if (length(common_samples) < 5) {
    stop(method_name, ": too few common samples between matrix and metadata.")
  }

  dt <- dt[, c("Tax", common_samples), with = FALSE]

  long <- melt(
    dt,
    id.vars = "Tax",
    variable.name = "SampleID",
    value.name = "Abundance"
  )

  md <- metadata[
    SampleID %in% common_samples,
    .(
      SampleID,
      BioProject,
      CenterName,
      CenterName_standardized,
      Continent,
      Country,
      Platform
    )
  ]

  long <- merge(long, md, by = "SampleID", all.x = TRUE)

  # Aggregate strictly to one row per BioProject x Tax.
  # Platform is descriptive metadata only and must not define the
  # independent analytical unit. Multi-platform projects are therefore
  # collapsed to a single BioProject-level profile.
  bp <- long[
    !is.na(BioProject),
    .(
      mean_abundance = mean(Abundance, na.rm = TRUE),
      prevalence_sra = mean(Abundance > 0, na.rm = TRUE),
      n_sra = uniqueN(SampleID),
      CenterName_standardized = paste(sort(unique(na.omit(CenterName_standardized))), collapse = ";"),
      Continent = paste(sort(unique(na.omit(Continent))), collapse = ";"),
      Country = paste(sort(unique(na.omit(Country))), collapse = ";"),
      Platform = paste(sort(unique(na.omit(Platform))), collapse = ";")
    ),
    by = .(
      BioProject,
      Tax
    )
  ]

  bp_wide <- dcast(
    bp,
    BioProject + CenterName_standardized + Continent + Country + Platform ~ Tax,
    value.var = "mean_abundance",
    fill = 0
  )

  fwrite(bp, file.path(out_dir, paste0(method_name, "_BioProject_taxa_long.csv")))
  fwrite(bp_wide, file.path(out_dir, paste0(method_name, "_BioProject_taxa_wide.csv")))

  bp_wide
}

summarise_bioproject_prevalence <- function(bp_wide, method_name) {

  meta_cols <- c("BioProject", "CenterName_standardized", "Continent", "Country", "Platform")
  tax_cols <- setdiff(names(bp_wide), meta_cols)

  prev <- data.table(
    Tax = tax_cols,
    prevalence_bioproject = sapply(
      tax_cols,
      function(x) mean(bp_wide[[x]] > detection_threshold, na.rm = TRUE)
    ),
    mean_abundance_bioproject = sapply(
      tax_cols,
      function(x) mean(bp_wide[[x]], na.rm = TRUE)
    ),
    median_abundance_bioproject = sapply(
      tax_cols,
      function(x) median(bp_wide[[x]], na.rm = TRUE)
    ),
    n_bioproject_detected = sapply(
      tax_cols,
      function(x) sum(bp_wide[[x]] > detection_threshold, na.rm = TRUE)
    )
  )

  prev <- prev[order(-prevalence_bioproject, -mean_abundance_bioproject)]

  fwrite(
    prev,
    file.path(out_dir, paste0(method_name, "_BioProject_prevalence_summary.csv"))
  )

  audit_list <- lapply(
    names(prevalence_audit_thresholds),
    function(th_name){

      th <- prevalence_audit_thresholds[[th_name]]

      data.table(
        Tax = tax_cols,
        Threshold = th_name,
        Prevalence_percent = sapply(
          tax_cols,
          function(x) mean(bp_wide[[x]] > th, na.rm = TRUE) * 100
        ),
        Positive_BioProjects = sapply(
          tax_cols,
          function(x) sum(bp_wide[[x]] > th, na.rm = TRUE)
        )
      )
    }
  )

  audit_dt <- rbindlist(audit_list)

  fwrite(
    audit_dt,
    file.path(out_dir, paste0(method_name, "_prevalence_threshold_audit.csv"))
  )

  threshold_summary <- audit_dt[
    ,
    .(
      Mean_prevalence = mean(Prevalence_percent, na.rm = TRUE),
      Median_prevalence = median(Prevalence_percent, na.rm = TRUE),
      Genera_detected = sum(Positive_BioProjects > 0)
    ),
    by = Threshold
  ]

  fwrite(
    threshold_summary,
    file.path(out_dir, paste0(method_name, "_threshold_summary.csv"))
  )

  core <- prev[
    prevalence_bioproject >= min_bioprojects_for_core &
      mean_abundance_bioproject >= min_mean_abundance_core
  ]

  fwrite(
    core,
    file.path(out_dir, paste0(method_name, "_BioProject_core_taxa.csv"))
  )

  p <- prev[1:min(top_n_taxa, .N)] %>%
    mutate(Tax = reorder(Tax, prevalence_bioproject)) %>%
    ggplot(aes(x = prevalence_bioproject * 100, y = Tax)) +
    geom_col() +
    labs(
      title = paste0(method_name, ": dominant taxa prevalence across BioProjects"),
      subtitle = paste0("Detection threshold = ", detection_threshold * 100, "%"),
      x = "BioProject prevalence (%)",
      y = NULL
    ) +
    theme_minimal(base_size = 13)

  ggsave(file.path(out_dir, paste0(method_name, "_BioProject_top_taxa_prevalence.png")), p, width = 8, height = 7, dpi = 600)
  ggsave(file.path(out_dir, paste0(method_name, "_BioProject_top_taxa_prevalence.pdf")), p, width = 8, height = 7)

  prev
}

run_bioproject_pcoa <- function(bp_wide, method_name) {

  meta_cols <- c("BioProject", "CenterName_standardized", "Continent", "Country", "Platform")
  tax_cols <- setdiff(names(bp_wide), meta_cols)

  mat <- as.matrix(bp_wide[, ..tax_cols])
  rownames(mat) <- bp_wide$BioProject
  mat[is.na(mat)] <- 0

  mat <- mat[, colSums(mat, na.rm = TRUE) > 0, drop = FALSE]

  message(method_name, ": BioProjects before filtering = ", nrow(mat))
  message(method_name, ": taxa before filtering = ", ncol(mat))
  message(method_name, ": empty BioProjects = ", sum(rowSums(mat) == 0))
  message(method_name, ": empty taxa = ", sum(colSums(mat) == 0))

  keep_rows <- rowSums(mat, na.rm = TRUE) > 0

  if (sum(!keep_rows) > 0) {
    fwrite(
      data.table(BioProject = rownames(mat)[!keep_rows]),
      file.path(out_dir, paste0(method_name, "_empty_BioProjects_removed_before_PCoA.csv"))
    )
    message(method_name, ": removing empty BioProjects before Bray-Curtis PCoA: ", sum(!keep_rows))
  }

  mat <- mat[keep_rows, , drop = FALSE]

  meta <- bp_wide[BioProject %in% rownames(mat), ..meta_cols]
  meta <- meta[match(rownames(mat), BioProject)]

  if (nrow(mat) < 3) {
    warning(method_name, ": fewer than 3 non-empty BioProjects; skipping PCoA.")
    return(NULL)
  }

  dist_bc <- vegdist(mat, method = "bray")
  dist_mat <- as.matrix(dist_bc)

  if (anyNA(dist_mat)) {
    warning(method_name, ": NA values in Bray-Curtis distance; skipping PCoA.")
    return(NULL)
  }

  pcoa <- cmdscale(dist_bc, eig = TRUE, k = 2)

  coords <- data.table(
    BioProject = rownames(mat),
    PCoA1 = pcoa$points[, 1],
    PCoA2 = pcoa$points[, 2]
  )

  coords <- merge(coords, meta, by = "BioProject", all.x = TRUE)

  eig <- pcoa$eig
  var1 <- round(100 * eig[1] / sum(eig[eig > 0]), 2)
  var2 <- round(100 * eig[2] / sum(eig[eig > 0]), 2)

  fwrite(coords, file.path(out_dir, paste0(method_name, "_BioProject_PCoA_coordinates.csv")))

  formula_terms <- c()

  if (uniqueN(na.omit(meta$Continent)) >= 2) formula_terms <- c(formula_terms, "Continent")
  if (uniqueN(na.omit(meta$Platform)) >= 2) formula_terms <- c(formula_terms, "Platform")

  if (length(formula_terms) > 0) {

    formula_obj <- as.formula(paste("dist_bc ~", paste(formula_terms, collapse = " + ")))

    perm <- adonis2(
      formula_obj,
      data = as.data.frame(meta),
      permutations = 999,
      by = "margin",
      parallel = n_cores
    )

    perm_dt <- as.data.table(as.data.frame(perm), keep.rownames = "Term")
    fwrite(perm_dt, file.path(out_dir, paste0(method_name, "_BioProject_PERMANOVA_Continent_Platform.csv")))
  }

  p <- ggplot(
    coords,
    aes(x = PCoA1, y = PCoA2, color = Continent, shape = Platform)
  ) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(aes(label = BioProject), size = 2.2, max.overlaps = 20) +
    labs(
      title = paste0(method_name, ": BioProject-level Bray-Curtis PCoA"),
      x = paste0("PCoA1 (", var1, "%)"),
      y = paste0("PCoA2 (", var2, "%)")
    ) +
    theme_minimal(base_size = 13)

  ggsave(file.path(out_dir, paste0(method_name, "_BioProject_PCoA.png")), p, width = 9, height = 7, dpi = 600)
  ggsave(file.path(out_dir, paste0(method_name, "_BioProject_PCoA.pdf")), p, width = 9, height = 7)

  invisible(coords)
}

starter_nonstarter_analysis <- function(bp_wide, method_name) {

  meta_cols <- c("BioProject", "CenterName_standardized", "Continent", "Country", "Platform")
  tax_cols <- setdiff(names(bp_wide), meta_cols)

  selected <- intersect(c(starter_genera, nonstarter_genera), tax_cols)

  if (length(selected) == 0) {
    warning(method_name, ": no starter/non-starter genera present in matrix.")
    return(NULL)
  }

  dt <- bp_wide[, c(meta_cols, selected), with = FALSE]

  long <- melt(
    dt,
    id.vars = meta_cols,
    variable.name = "Genus",
    value.name = "Abundance"
  )

  long[, Group := fifelse(
    Genus %in% starter_genera,
    "Starter-associated",
    "Non-starter-associated"
  )]

  stat_group <- tryCatch(wilcox.test(Abundance ~ Group, data = long), error = function(e) NULL)
  stat_genus <- tryCatch(kruskal.test(Abundance ~ Genus, data = long), error = function(e) NULL)

  stats_out <- data.table(
    method = method_name,
    test = c("Wilcoxon starter_vs_nonstarter", "Kruskal-Wallis among_genera"),
    statistic = c(
      ifelse(is.null(stat_group), NA_real_, unname(stat_group$statistic)),
      ifelse(is.null(stat_genus), NA_real_, unname(stat_genus$statistic))
    ),
    p_value = c(
      ifelse(is.null(stat_group), NA_real_, stat_group$p.value),
      ifelse(is.null(stat_genus), NA_real_, stat_genus$p.value)
    )
  )

  fwrite(stats_out, file.path(out_dir, paste0(method_name, "_starter_nonstarter_tests.csv")))
  fwrite(long, file.path(out_dir, paste0(method_name, "_starter_nonstarter_long.csv")))

  p1 <- ggplot(
    long,
    aes(x = Genus, y = Abundance * 100, fill = Group)
  ) +
    geom_boxplot(outlier.alpha = 0.35) +
    geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
    scale_y_continuous(trans = "log1p") +
    labs(
      title = paste0(method_name, ": distribution of starter and non-starter genera"),
      x = NULL,
      y = "Mean relative abundance per BioProject (%)"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(out_dir, paste0(method_name, "_starter_nonstarter_distribution.png")), p1, width = 8, height = 6, dpi = 600)
  ggsave(file.path(out_dir, paste0(method_name, "_starter_nonstarter_distribution.pdf")), p1, width = 8, height = 6)

  invisible(stats_out)
}

wgs_species_richness_by_genus <- function(wgs_species_bp_wide, method_name) {

  meta_cols <- c("BioProject", "CenterName_standardized", "Continent", "Country", "Platform")
  species_cols <- setdiff(names(wgs_species_bp_wide), meta_cols)

  species_summary <- data.table(
    Species = species_cols,
    Genus = sub(" .*", "", species_cols),
    prevalence_bioproject = sapply(species_cols, function(x) mean(wgs_species_bp_wide[[x]] > 0)),
    n_bioproject_detected = sapply(species_cols, function(x) sum(wgs_species_bp_wide[[x]] > 0)),
    mean_abundance_bioproject = sapply(species_cols, function(x) mean(wgs_species_bp_wide[[x]]))
  )

  # Remove host / contaminant taxa from species-level richness summaries.
  species_summary <- species_summary[
    !Genus %in% excluded_host_taxa &
      !Species %in% excluded_host_taxa
  ]

  species_filtered <- species_summary[
    n_bioproject_detected >= species_detection_bioprojects &
      mean_abundance_bioproject >= species_detection_mean_abundance
  ]

  genus_richness <- species_filtered[
    ,
    .(
      n_species_detected = uniqueN(Species),
      mean_species_abundance = mean(mean_abundance_bioproject),
      max_species_prevalence = max(prevalence_bioproject),
      n_bioprojects_any_species = max(n_bioproject_detected)
    ),
    by = Genus
  ][order(-n_species_detected, -max_species_prevalence)]

  fwrite(species_summary, file.path(out_dir, paste0(method_name, "_species_detection_summary.csv")))
  fwrite(species_filtered, file.path(out_dir, paste0(method_name, "_species_detection_filtered.csv")))
  fwrite(genus_richness, file.path(out_dir, paste0(method_name, "_species_richness_by_genus.csv")))

  if (nrow(genus_richness) > 0) {
    p <- genus_richness[1:min(25, .N)] %>%
      mutate(Genus = reorder(Genus, n_species_detected)) %>%
      ggplot(aes(x = n_species_detected, y = Genus)) +
      geom_col() +
      labs(
        title = paste0(method_name, ": number of detected species within dominant genera"),
        x = "Number of species detected",
        y = NULL
      ) +
      theme_minimal(base_size = 13)

    ggsave(file.path(out_dir, paste0(method_name, "_species_richness_by_genus.png")), p, width = 8, height = 7, dpi = 600)
    ggsave(file.path(out_dir, paste0(method_name, "_species_richness_by_genus.pdf")), p, width = 8, height = 7)
  }

  genus_richness
}

method_concordance <- function(prev_a, prev_b, comparison_name) {

  # Raw table: all genera shared by both pipelines.
  cor_dt_raw <- merge(
    prev_a[, .(
      Tax,
      A_prevalence = prevalence_bioproject,
      A_mean_abundance = mean_abundance_bioproject
    )],
    prev_b[, .(
      Tax,
      B_prevalence = prevalence_bioproject,
      B_mean_abundance = mean_abundance_bioproject
    )],
    by = "Tax"
  )

  fwrite(
    cor_dt_raw,
    file.path(out_dir, paste0(comparison_name, "_genus_concordance_raw.csv"))
  )

  # Filter for biologically interpretable correlation/plot.
  # Keep genera detected in >= min_prevalence_for_correlation BioProjects in at least one pipeline
  # and with mean abundance >= min_mean_abundance_for_correlation in at least one pipeline.
  cor_dt <- cor_dt_raw[
    (
      A_prevalence >= min_prevalence_for_correlation |
        B_prevalence >= min_prevalence_for_correlation
    ) &
      (
        A_mean_abundance >= min_mean_abundance_for_correlation |
          B_mean_abundance >= min_mean_abundance_for_correlation
      )
  ]

  fwrite(
    cor_dt,
    file.path(out_dir, paste0(comparison_name, "_genus_concordance_filtered.csv"))
  )

  qc <- data.table(
    comparison = comparison_name,
    n_genera_before_filter = nrow(cor_dt_raw),
    n_genera_after_filter = nrow(cor_dt),
    prevalence_threshold = min_prevalence_for_correlation,
    mean_abundance_threshold = min_mean_abundance_for_correlation,
    n_A_prevalence_nonzero_raw = sum(cor_dt_raw$A_prevalence > 0),
    n_B_prevalence_nonzero_raw = sum(cor_dt_raw$B_prevalence > 0),
    n_A_prevalence_nonzero_filtered = sum(cor_dt$A_prevalence > 0),
    n_B_prevalence_nonzero_filtered = sum(cor_dt$B_prevalence > 0)
  )

  fwrite(
    qc,
    file.path(out_dir, paste0(comparison_name, "_concordance_qc.csv"))
  )

  if (nrow(cor_dt) < 3) {
    warning(comparison_name, ": fewer than 3 genera after filtering; skipping concordance stats.")
    return(NULL)
  }

  cor_prev <- cor.test(
    cor_dt$A_prevalence,
    cor_dt$B_prevalence,
    method = "spearman",
    exact = FALSE
  )

  cor_abund <- cor.test(
    cor_dt$A_mean_abundance,
    cor_dt$B_mean_abundance,
    method = "spearman",
    exact = FALSE
  )

  stats <- data.table(
    comparison = comparison_name,
    metric = c("BioProject prevalence", "Mean relative abundance"),
    method = "Spearman",
    rho = c(unname(cor_prev$estimate), unname(cor_abund$estimate)),
    p_value = c(cor_prev$p.value, cor_abund$p.value),
    p_value_scientific = c(
      formatC(cor_prev$p.value, format = "e", digits = 2),
      formatC(cor_abund$p.value, format = "e", digits = 2)
    ),
    n_common_genera_raw = nrow(cor_dt_raw),
    n_common_genera_filtered = nrow(cor_dt),
    prevalence_threshold = min_prevalence_for_correlation,
    mean_abundance_threshold = min_mean_abundance_for_correlation
  )

  fwrite(
    stats,
    file.path(out_dir, paste0(comparison_name, "_genus_concordance_stats.csv"))
  )

  # Prevalence plot based on filtered genera.
  p1 <- ggplot(
    cor_dt,
    aes(x = A_prevalence * 100, y = B_prevalence * 100, label = Tax)
  ) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_text_repel(size = 2.5, max.overlaps = 20) +
    labs(
      title = paste0(comparison_name, ": genus prevalence concordance"),
      subtitle = paste0(
        "Filtered genera: n = ", nrow(cor_dt),
        " / raw n = ", nrow(cor_dt_raw),
        "; Spearman rho = ", round(unname(cor_prev$estimate), 3),
        ", p = ", formatC(cor_prev$p.value, format = "e", digits = 2)
      ),
      x = "Method A BioProject prevalence (%)",
      y = "Method B BioProject prevalence (%)"
    ) +
    theme_minimal(base_size = 13)

  ggsave(
    file.path(out_dir, paste0(comparison_name, "_genus_prevalence_correlation_filtered.png")),
    p1,
    width = 7,
    height = 6,
    dpi = 600
  )

  ggsave(
    file.path(out_dir, paste0(comparison_name, "_genus_prevalence_correlation_filtered.pdf")),
    p1,
    width = 7,
    height = 6
  )

  # Abundance plot based on filtered genera.
  p2 <- ggplot(
    cor_dt,
    aes(x = A_mean_abundance * 100, y = B_mean_abundance * 100, label = Tax)
  ) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_text_repel(size = 2.5, max.overlaps = 20) +
    scale_x_continuous(trans = "log1p") +
    scale_y_continuous(trans = "log1p") +
    labs(
      title = paste0(comparison_name, ": genus abundance concordance"),
      subtitle = paste0(
        "Filtered genera: n = ", nrow(cor_dt),
        " / raw n = ", nrow(cor_dt_raw),
        "; Spearman rho = ", round(unname(cor_abund$estimate), 3),
        ", p = ", formatC(cor_abund$p.value, format = "e", digits = 2)
      ),
      x = "Method A mean relative abundance per BioProject (%)",
      y = "Method B mean relative abundance per BioProject (%)"
    ) +
    theme_minimal(base_size = 13)

  ggsave(
    file.path(out_dir, paste0(comparison_name, "_genus_abundance_correlation_filtered.png")),
    p2,
    width = 7,
    height = 6,
    dpi = 600
  )

  ggsave(
    file.path(out_dir, paste0(comparison_name, "_genus_abundance_correlation_filtered.pdf")),
    p2,
    width = 7,
    height = 6
  )

  invisible(stats)
}

# ============================================================
# LOAD FINAL STAGE5 METADATA WITH CURATED CENTER MAP
# ============================================================

message("Loading final stage5 metadata:")
message("- ", stage5_amplicon_metadata_path)
message("- ", stage5_wgs_metadata_path)

if (!file.exists(stage5_amplicon_metadata_path)) {
  stop("Missing stage5 amplicon metadata: ", stage5_amplicon_metadata_path)
}
if (!file.exists(stage5_wgs_metadata_path)) {
  stop("Missing stage5 WGS metadata: ", stage5_wgs_metadata_path)
}

metadata_amplicon <- fread(stage5_amplicon_metadata_path)
metadata_wgs <- fread(stage5_wgs_metadata_path)

metadata_amplicon[, sequencing_subset := "AMPLICON"]
metadata_wgs[, sequencing_subset := "WGS"]

metadata <- rbindlist(
  list(metadata_amplicon, metadata_wgs),
  use.names = TRUE,
  fill = TRUE
)

if ("Run" %in% names(metadata)) {
  metadata[, SampleID := as.character(Run)]
} else if ("SampleID" %in% names(metadata)) {
  metadata[, SampleID := as.character(SampleID)]
} else if ("Sample" %in% names(metadata)) {
  metadata[, SampleID := as.character(Sample)]
}

if (!"CenterName" %in% names(metadata)) metadata[, CenterName := NA_character_]
if (!"Country" %in% names(metadata)) metadata[, Country := NA_character_]
if (!"Continent" %in% names(metadata)) metadata[, Continent := NA_character_]
if (!"Platform" %in% names(metadata)) metadata[, Platform := NA_character_]
if (!"LibraryStrategy" %in% names(metadata)) metadata[, LibraryStrategy := sequencing_subset]

metadata[, CenterName := as.character(CenterName)]
metadata[, CenterName_norm := normalize_center_name(CenterName)]
metadata[, Country_original := as.character(Country)]
metadata[, Continent_original := as.character(Continent)]
metadata[, CenterName_standardized := CenterName_norm]

if (file.exists(center_country_map_path)) {

  message("Loading curated center-country map: ", center_country_map_path)

  center_map <- fread(center_country_map_path)

  if (!all(c("CenterName", "CenterName_standardized", "Country", "Continent") %in% names(center_map))) {
    stop("center_country_map_curated.csv must contain: CenterName, CenterName_standardized, Country, Continent")
  }

  center_map[, CenterName := as.character(CenterName)]
  center_map[, CenterName_norm := normalize_center_name(CenterName)]
  center_map[, CenterName_standardized := as.character(CenterName_standardized)]
  center_map[, Country_curated := as.character(Country)]
  center_map[, Continent_curated := as.character(Continent)]

  center_map_small <- unique(
    center_map[
      ,
      .(
        CenterName_norm,
        CenterName_standardized_curated = CenterName_standardized,
        Country_curated,
        Continent_curated
      )
    ],
    by = "CenterName_norm"
  )

  n_before <- nrow(metadata)
  bp_before <- uniqueN(metadata$BioProject)
  sra_before <- uniqueN(metadata$SampleID)

  metadata <- merge(
    metadata,
    center_map_small,
    by = "CenterName_norm",
    all.x = TRUE
  )

  if (nrow(metadata) != n_before || uniqueN(metadata$BioProject) != bp_before || uniqueN(metadata$SampleID) != sra_before) {
    stop("CenterName mapping changed the number of rows, BioProjects, or SRA records. Check duplicated CenterName_norm in center_country_map_curated.csv.")
  }

  metadata[
    !is.na(CenterName_standardized_curated) & CenterName_standardized_curated != "",
    CenterName_standardized := CenterName_standardized_curated
  ]

  metadata[
    !is.na(Country_curated) & Country_curated != "",
    Country := Country_curated
  ]

  metadata[
    !is.na(Continent_curated) & Continent_curated != "",
    Continent := Continent_curated
  ]

} else {
  message("Curated center-country map not found. Using metadata Country/Continent as provided.")
}

metadata[is.na(Country) | Country == "", Country := "Unknown"]
metadata[is.na(Continent) | Continent == "", Continent := "Unknown"]
metadata[is.na(CenterName_standardized) | CenterName_standardized == "", CenterName_standardized := "Unknown"]
metadata[is.na(Platform) | Platform == "", Platform := "Unknown"]

required_md <- c("SampleID", "BioProject")
if (!all(required_md %in% names(metadata))) {
  stop("Metadata must contain SampleID/Run and BioProject columns.")
}

metadata <- metadata[!is.na(SampleID) & !is.na(BioProject)]
metadata[, SampleID := normalize_sample_names(as.character(SampleID))]
metadata[, BioProject := as.character(BioProject)]
metadata[, Continent := as.character(Continent)]
metadata[, Country := as.character(Country)]
metadata[, Platform := as.character(Platform)]
metadata[, CenterName_standardized := as.character(CenterName_standardized)]
metadata[, LibraryStrategy := toupper(as.character(LibraryStrategy))]
metadata[is.na(LibraryStrategy) | LibraryStrategy == "", LibraryStrategy := sequencing_subset]
metadata[, sequencing_subset := toupper(as.character(sequencing_subset))]

# One row per SRA accession in the final analytical metadata.
metadata <- metadata[!duplicated(SampleID)]

fwrite(
  unique(
    metadata[
      ,
      .(
        CenterName,
        CenterName_standardized,
        Country,
        Continent,
        Country_original,
        Continent_original
      )
    ]
  ),
  file.path(out_dir, "CenterName_country_mapping_used.csv")
)

# BioProject-level metadata table.
bioproject_metadata_summary <- metadata[
  ,
  .(
    n_sra = .N,
    n_amplicon_sra = sum(sequencing_subset == "AMPLICON", na.rm = TRUE),
    n_wgs_sra = sum(sequencing_subset == "WGS", na.rm = TRUE),
    CenterName_standardized = unique(na.omit(CenterName_standardized))[1],
    Country = unique(na.omit(Country))[1],
    Continent = unique(na.omit(Continent))[1],
    Platform = unique(na.omit(Platform))[1]
  ),
  by = BioProject
][order(-n_sra)]

fwrite(
  bioproject_metadata_summary,
  file.path(out_dir, "BioProject_SRA_counts_metadata_summary.csv")
)

# ----------------------------------------------------------------
# v7/v8 geographic summaries: SRA and BioProject counts per Country and Continent
# ----------------------------------------------------------------

country_summary <- metadata[
  ,
  .(
    n_sra = .N,
    n_bioprojects = uniqueN(BioProject),
    n_amplicon_sra = sum(sequencing_subset == "AMPLICON", na.rm = TRUE),
    n_wgs_sra = sum(sequencing_subset == "WGS", na.rm = TRUE),
    Continents = paste(sort(unique(Continent)), collapse = "; ")
  ),
  by = Country
][order(-n_sra, -n_bioprojects, Country)]

continent_summary <- metadata[
  ,
  .(
    n_sra = .N,
    n_bioprojects = uniqueN(BioProject),
    n_countries = uniqueN(Country),
    n_amplicon_sra = sum(sequencing_subset == "AMPLICON", na.rm = TRUE),
    n_wgs_sra = sum(sequencing_subset == "WGS", na.rm = TRUE),
    Countries = paste(sort(unique(Country)), collapse = "; ")
  ),
  by = Continent
][order(-n_sra, -n_bioprojects, Continent)]

country_continent_summary <- metadata[
  ,
  .(
    n_sra = .N,
    n_bioprojects = uniqueN(BioProject),
    n_amplicon_sra = sum(sequencing_subset == "AMPLICON", na.rm = TRUE),
    n_wgs_sra = sum(sequencing_subset == "WGS", na.rm = TRUE)
  ),
  by = .(Country, Continent)
][order(-n_sra, -n_bioprojects, Country, Continent)]

fwrite(
  country_summary,
  file.path(out_dir, "Country_SRA_BioProject_summary.csv")
)

fwrite(
  continent_summary,
  file.path(out_dir, "Continent_SRA_BioProject_summary.csv")
)

fwrite(
  country_continent_summary,
  file.path(out_dir, "Country_Continent_SRA_BioProject_summary.csv")
)

message("Geographic summaries saved:")
message("- ", file.path(out_dir, "Country_SRA_BioProject_summary.csv"))
message("- ", file.path(out_dir, "Continent_SRA_BioProject_summary.csv"))
message("- ", file.path(out_dir, "Country_Continent_SRA_BioProject_summary.csv"))

message("Metadata samples after final stage5 curation: ", uniqueN(metadata$SampleID))
message("BioProjects after final stage5 curation: ", uniqueN(metadata$BioProject))
message("Countries after final stage5 curation: ", uniqueN(metadata$Country))
message("Continents after final stage5 curation: ", uniqueN(metadata$Continent))
message("Standardized institutions after final stage5 curation: ", uniqueN(metadata$CenterName_standardized))

# ============================================================
# BUILD KRAKEN2 MATRICES IF REQUESTED
# ============================================================

if (build_kraken_matrix) {

  build_kraken_matrix_from_long(
    kraken_dir = kraken_taxa_long_dir,
    out_path = kraken_genus_matrix_path,
    rank_level = "genus",
    n_cores = n_cores
  )

  build_kraken_matrix_from_long(
    kraken_dir = kraken_taxa_long_dir,
    out_path = kraken_species_matrix_path,
    rank_level = "species",
    n_cores = n_cores
  )

} else {
  message("Skipping Kraken2 matrix construction. Using existing matrices:")
  message(kraken_genus_matrix_path)
  message(kraken_species_matrix_path)
}

# ============================================================
# QIIME2 AMPLICON
# ============================================================

qiime2_raw <- read_tax_matrix(qiime2_amplicon_matrix_path)
qiime2_rel <- to_relative_abundance(qiime2_raw)
qiime2_genus <- collapse_to_rank(copy(qiime2_rel), rank = "genus")

# ============================================================
# GLOBAL AUDIT: stage5 metadata vs QIIME2 tax_matrix.tsv
# ============================================================

qiime2_samples <- setdiff(names(qiime2_genus), "Tax")

metadata_audit <- copy(metadata)

metadata_audit[
  ,
  in_qiime2_matrix := SampleID %in% qiime2_samples
]

sra_audit <- metadata_audit[
  ,
  .(
    n_sra_metadata = .N,
    n_sra_qiime2 = sum(in_qiime2_matrix),
    n_sra_missing = sum(!in_qiime2_matrix),
    pct_missing = round(
      100 * sum(!in_qiime2_matrix) / .N,
      2
    )
  )
]

fwrite(
  sra_audit,
  file.path(out_dir, "QIIME2_stage5_vs_tax_matrix_global_audit.csv")
)

bp_audit <- metadata_audit[
  ,
  .(
    n_sra_metadata = .N,
    n_sra_qiime2 = sum(in_qiime2_matrix),
    n_sra_missing = sum(!in_qiime2_matrix)
  ),
  by = BioProject
]

bp_audit[
  ,
  present_in_qiime2 := n_sra_qiime2 > 0
]

bp_audit <- merge(
  bp_audit,
  unique(
    metadata_audit[
      ,
      .(
        BioProject,
        Country,
        Continent,
        CenterName_standardized
      )
    ]
  ),
  by = "BioProject",
  all.x = TRUE
)

fwrite(
  bp_audit,
  file.path(out_dir, "QIIME2_stage5_vs_tax_matrix_BioProject_audit.csv")
)

message("")
message("===== GLOBAL QIIME2 AUDIT =====")
print(sra_audit)

message("")
message("BioProjects in metadata: ",
        uniqueN(metadata_audit$BioProject))

message("BioProjects represented in QIIME2 matrix: ",
        bp_audit[present_in_qiime2 == TRUE, uniqueN(BioProject)])

message("BioProjects completely missing from QIIME2 matrix: ",
        bp_audit[present_in_qiime2 == FALSE, uniqueN(BioProject)])

message("Missing SRA records: ",
        sum(bp_audit$n_sra_missing))

message("===============================")
message("")

# ============================================================
# DIAGNOSTIC: Turkish amplicon BioProjects missing from QIIME2 matrix
# ============================================================

qiime2_sample_cols <- setdiff(names(qiime2_genus), "Tax")

turkey_amplicon_metadata <- metadata[
  Country == "Turkey",
  .(
    BioProject,
    SampleID,
    CenterName,
    CenterName_standardized,
    Country,
    Continent,
    in_qiime2_matrix = SampleID %in% qiime2_sample_cols
  )
]

fwrite(
  turkey_amplicon_metadata,
  file.path(out_dir, "DIAGNOSTIC_Turkey_samples_in_QIIME2_matrix.csv")
)

turkey_amplicon_bioprojects <- turkey_amplicon_metadata[
  ,
  .(
    n_sra_metadata = uniqueN(SampleID),
    n_sra_in_qiime2_matrix = sum(in_qiime2_matrix),
    present_in_qiime2_matrix = any(in_qiime2_matrix),
    CenterName_standardized = unique(na.omit(CenterName_standardized))[1],
    Country = unique(na.omit(Country))[1],
    Continent = unique(na.omit(Continent))[1]
  ),
  by = BioProject
][order(BioProject)]

fwrite(
  turkey_amplicon_bioprojects,
  file.path(out_dir, "DIAGNOSTIC_Turkey_BioProjects_in_QIIME2_matrix.csv")
)

message("Turkey BioProjects in metadata:")
print(turkey_amplicon_bioprojects)

qiime2_bp <- aggregate_to_bioproject(qiime2_genus, metadata, "QIIME2_Amplicon_genus")
qiime2_prev <- summarise_bioproject_prevalence(qiime2_bp, "QIIME2_Amplicon_genus")
run_bioproject_pcoa(qiime2_bp, "QIIME2_Amplicon_genus")
starter_nonstarter_analysis(qiime2_bp, "QIIME2_Amplicon_genus")

# ============================================================
# KRAKEN2 WGS
# ============================================================

kraken_genus_raw <- read_tax_matrix(kraken_genus_matrix_path)
kraken_species_raw <- read_tax_matrix(kraken_species_matrix_path)

kraken_genus_rel <- to_relative_abundance(kraken_genus_raw)
kraken_species_rel <- to_relative_abundance(kraken_species_raw)

kraken_genus_bp <- aggregate_to_bioproject(kraken_genus_rel, metadata, "Kraken2_WGS_genus")
kraken_species_bp <- aggregate_to_bioproject(kraken_species_rel, metadata, "Kraken2_WGS_species")

kraken_prev <- summarise_bioproject_prevalence(kraken_genus_bp, "Kraken2_WGS_genus")
run_bioproject_pcoa(kraken_genus_bp, "Kraken2_WGS_genus")
starter_nonstarter_analysis(kraken_genus_bp, "Kraken2_WGS_genus")
wgs_species_richness_by_genus(kraken_species_bp, "Kraken2_WGS_species")

# ============================================================
# mOTUs WGS
# ============================================================

motus_raw <- read_tax_matrix(motus_matrix_path)
motus_rel <- to_relative_abundance(motus_raw)
motus_genus <- collapse_to_rank(copy(motus_rel), rank = "genus")
motus_species <- collapse_to_rank(copy(motus_rel), rank = "species")

motus_genus_bp <- aggregate_to_bioproject(motus_genus, metadata, "mOTUs_WGS_genus")
motus_species_bp <- aggregate_to_bioproject(motus_species, metadata, "mOTUs_WGS_species")

motus_prev <- summarise_bioproject_prevalence(motus_genus_bp, "mOTUs_WGS_genus")
run_bioproject_pcoa(motus_genus_bp, "mOTUs_WGS_genus")
starter_nonstarter_analysis(motus_genus_bp, "mOTUs_WGS_genus")
wgs_species_richness_by_genus(motus_species_bp, "mOTUs_WGS_species")

# ============================================================
# METHOD CONCORDANCE
# ============================================================

# QIIME2 vs WGS comparisons are aggregate cross-collection comparisons.
# The underlying BioProject collections are independent and therefore
# these comparisons describe concordance of aggregate taxonomic patterns,
# not paired method agreement.
method_concordance(qiime2_prev, kraken_prev, "QIIME2_vs_Kraken2")
method_concordance(qiime2_prev, motus_prev, "QIIME2_vs_mOTUs")

# Kraken2 and mOTUs were applied to the same WGS collection.
# Restrict both matrices to the identical set of BioProjects before
# calculating prevalence and abundance concordance.
common_wgs_bioprojects <- intersect(
  kraken_genus_bp$BioProject,
  motus_genus_bp$BioProject
)

kraken_genus_bp_paired <- kraken_genus_bp[
  BioProject %in% common_wgs_bioprojects
]

motus_genus_bp_paired <- motus_genus_bp[
  BioProject %in% common_wgs_bioprojects
]

message(
  "Paired Kraken2-mOTUs concordance: ",
  length(common_wgs_bioprojects),
  " common BioProjects"
)

message(
  "Excluded from paired Kraken2 set: ",
  paste(
    setdiff(kraken_genus_bp$BioProject, common_wgs_bioprojects),
    collapse = ", "
  )
)

kraken_prev_paired <- summarise_bioproject_prevalence(
  kraken_genus_bp_paired,
  "Kraken2_WGS_genus_paired19"
)

motus_prev_paired <- summarise_bioproject_prevalence(
  motus_genus_bp_paired,
  "mOTUs_WGS_genus_paired19"
)

method_concordance(
  kraken_prev_paired,
  motus_prev_paired,
  "Kraken2_vs_mOTUs"
)

# ============================================================
# DATASET AUDIT PLOTS
# ============================================================

bp_counts <- fread(file.path(out_dir, "BioProject_SRA_counts_metadata_summary.csv"))

p_counts <- ggplot(bp_counts, aes(x = n_sra)) +
  geom_histogram(bins = 40) +
  scale_x_continuous(trans = "log10") +
  labs(
    title = "Uneven distribution of SRA records across BioProjects",
    x = "Number of SRA records per BioProject, log10 scale",
    y = "Number of BioProjects"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(out_dir, "BioProject_SRA_count_distribution.png"), p_counts, width = 7, height = 5, dpi = 600)
ggsave(file.path(out_dir, "BioProject_SRA_count_distribution.pdf"), p_counts, width = 7, height = 5)

top_bp <- bp_counts[1:min(20, .N)]

p_top <- ggplot(top_bp, aes(x = reorder(BioProject, n_sra), y = n_sra)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top BioProjects by number of SRA records",
    x = NULL,
    y = "Number of SRA records"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(out_dir, "Top_BioProjects_by_SRA_count.png"), p_top, width = 8, height = 6, dpi = 600)
ggsave(file.path(out_dir, "Top_BioProjects_by_SRA_count.pdf"), p_top, width = 8, height = 6)


# ============================================================
# PUBLICATION TABLES AND FINAL ANALYTICAL DATASET SUMMARY
# ============================================================

fmt_int <- function(x) {
  ifelse(is.na(x), NA_character_, format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE))
}

write_table_all <- function(dt, basename) {
  csv_path <- file.path(out_dir, paste0(basename, ".csv"))
  tsv_path <- file.path(out_dir, paste0(basename, ".tsv"))
  fwrite(dt, csv_path)
  fwrite(dt, tsv_path, sep = "\t")

  xlsx_path <- file.path(out_dir, paste0(basename, ".xlsx"))
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    openxlsx::write.xlsx(as.data.frame(dt), xlsx_path, overwrite = TRUE)
  } else {
    message("Package openxlsx not installed; skipping XLSX: ", xlsx_path)
  }
}

bp_nonempty_count <- function(bp_wide) {
  meta_cols <- c("BioProject", "CenterName_standardized", "Continent", "Country", "Platform")
  tax_cols <- setdiff(names(bp_wide), meta_cols)
  if (length(tax_cols) == 0) return(0L)
  mat <- as.matrix(bp_wide[, ..tax_cols])
  mat[is.na(mat)] <- 0
  sum(rowSums(mat, na.rm = TRUE) > 0)
}

matrix_summary_from_objects <- function(dataset, raw_dt, rel_dt, bp_wide = NULL) {
  raw_cols <- setdiff(names(raw_dt), "Tax")
  rel_cols <- setdiff(names(rel_dt), "Tax")
  common <- intersect(rel_cols, metadata$SampleID)

  data.table(
    Dataset = dataset,
    Raw_matrix_columns_after_sample_normalization = length(raw_cols),
    Matrix_samples_after_zero_count_removal = length(rel_cols),
    Common_samples_with_stage5_metadata = length(common),
    BioProjects_represented_in_matrix = metadata[SampleID %in% common, uniqueN(BioProject)],
    BioProjects_after_BioProject_aggregation = if (!is.null(bp_wide)) nrow(bp_wide) else NA_integer_,
    Nonempty_BioProjects_used_for_PCoA = if (!is.null(bp_wide)) bp_nonempty_count(bp_wide) else NA_integer_
  )
}

# Static workflow counts from the metadata curation pipeline and dynamic final stage5 counts.
stage5_total_sra <- uniqueN(metadata$SampleID)
stage5_total_bp <- uniqueN(metadata$BioProject)
stage5_amp_sra <- metadata[sequencing_subset == "AMPLICON", uniqueN(SampleID)]
stage5_amp_bp <- metadata[sequencing_subset == "AMPLICON", uniqueN(BioProject)]
stage5_wgs_sra <- metadata[sequencing_subset == "WGS", uniqueN(SampleID)]
stage5_wgs_bp <- metadata[sequencing_subset == "WGS", uniqueN(BioProject)]

initial_count_check <- data.table(
  Source = c("Initial metadata file", "Final stage5 metadata", "Final stage5 amplicon", "Final stage5 WGS"),
  Path = c(initial_metadata_path, metadata_path, stage5_amplicon_metadata_path, stage5_wgs_metadata_path),
  SRA = c(
    if (file.exists(initial_metadata_path)) fread(initial_metadata_path, select = "Run")[, uniqueN(Run)] else NA_integer_,
    stage5_total_sra,
    stage5_amp_sra,
    stage5_wgs_sra
  ),
  BioProjects = c(
    if (file.exists(initial_metadata_path)) fread(initial_metadata_path, select = "BioProject")[, uniqueN(BioProject)] else NA_integer_,
    stage5_total_bp,
    stage5_amp_bp,
    stage5_wgs_bp
  )
)

write_table_all(initial_count_check, "Initial_and_final_metadata_count_check")

table1_numeric <- data.table(
  Metric = c(
    "Initial NCBI query",
    "Metadata-based filtering",
    "Automatic BioProject retention (>=80% cheese-associated SRA)",
    "Manual BioProject curation",
    "Final curated dataset",
    "Amplicon subset",
    "WGS subset"
  ),
  SRA = c(14203L, 10853L, 8961L, NA_integer_, stage5_total_sra, stage5_amp_sra, stage5_wgs_sra),
  BioProjects = c(233L, 163L, 92L, 45L, stage5_total_bp, stage5_amp_bp, stage5_wgs_bp),
  Notes = c(
    "Candidate records retrieved from the NCBI BioProject/SRA query.",
    "Records retained after run-level metadata whitelist/blacklist filtering.",
    "BioProjects automatically retained when >=80% of eligible SRA matched the whitelist.",
    "BioProjects requiring manual review; retained/rejected using BioProject descriptions, publications, and run annotations.",
    "Union of automatically retained BioProjects and manually accepted BioProjects; dynamic count from stage5 metadata.",
    "Dynamic count from stage5/runinfo_AMPLICON.csv.",
    "Dynamic count from stage5/runinfo_WGS.csv."
  )
)

# Minimal publication table with formatted numbers.
table1_publication_minimal <- copy(table1_numeric)
table1_publication_minimal[, SRA := fifelse(is.na(SRA), "--", fmt_int(SRA))]
table1_publication_minimal[, BioProjects := fmt_int(BioProjects)]
table1_publication_minimal[, Notes := NULL]

write_table_all(table1_numeric, "Table1_filtering_workflow_publication_with_notes")
write_table_all(table1_publication_minimal, "Table1_filtering_workflow_publication_minimal")

analysis_summary <- rbindlist(
  list(
    data.table(
      Dataset = "Final curated metadata",
      Raw_matrix_columns_after_sample_normalization = NA_integer_,
      Matrix_samples_after_zero_count_removal = stage5_total_sra,
      Common_samples_with_stage5_metadata = stage5_total_sra,
      BioProjects_represented_in_matrix = stage5_total_bp,
      BioProjects_after_BioProject_aggregation = NA_integer_,
      Nonempty_BioProjects_used_for_PCoA = NA_integer_
    ),
    data.table(
      Dataset = "Amplicon subset in curated metadata",
      Raw_matrix_columns_after_sample_normalization = NA_integer_,
      Matrix_samples_after_zero_count_removal = stage5_amp_sra,
      Common_samples_with_stage5_metadata = stage5_amp_sra,
      BioProjects_represented_in_matrix = stage5_amp_bp,
      BioProjects_after_BioProject_aggregation = NA_integer_,
      Nonempty_BioProjects_used_for_PCoA = NA_integer_
    ),
    data.table(
      Dataset = "WGS subset in curated metadata",
      Raw_matrix_columns_after_sample_normalization = NA_integer_,
      Matrix_samples_after_zero_count_removal = stage5_wgs_sra,
      Common_samples_with_stage5_metadata = stage5_wgs_sra,
      BioProjects_represented_in_matrix = stage5_wgs_bp,
      BioProjects_after_BioProject_aggregation = NA_integer_,
      Nonempty_BioProjects_used_for_PCoA = NA_integer_
    ),
    matrix_summary_from_objects("QIIME2 amplicon genus", qiime2_raw, qiime2_genus, qiime2_bp),
    matrix_summary_from_objects("Kraken2 WGS genus", kraken_genus_raw, kraken_genus_rel, kraken_genus_bp),
    matrix_summary_from_objects("Kraken2 WGS species", kraken_species_raw, kraken_species_rel, kraken_species_bp),
    matrix_summary_from_objects("mOTUs WGS genus", motus_raw, motus_genus, motus_genus_bp),
    matrix_summary_from_objects("mOTUs WGS species", motus_raw, motus_species, motus_species_bp)
  ),
  use.names = TRUE,
  fill = TRUE
)

write_table_all(analysis_summary, "Table2_final_analysis_dataset_summary_v7_metrics")

matrix_consistency <- rbindlist(
  list(
    data.table(
      Dataset = "QIIME2 amplicon genus matrix vs stage5 amplicon metadata",
      Expected_SRA_in_stage5 = stage5_amp_sra,
      Matrix_SRA_after_zero_count_removal = length(setdiff(names(qiime2_genus), "Tax")),
      Missing_in_matrix = length(setdiff(metadata[sequencing_subset == "AMPLICON", SampleID], setdiff(names(qiime2_genus), "Tax"))),
      Extra_in_matrix = length(setdiff(setdiff(names(qiime2_genus), "Tax"), metadata[sequencing_subset == "AMPLICON", SampleID]))
    ),
    data.table(
      Dataset = "Kraken2 WGS genus matrix vs stage5 WGS metadata",
      Expected_SRA_in_stage5 = stage5_wgs_sra,
      Matrix_SRA_after_zero_count_removal = length(setdiff(names(kraken_genus_rel), "Tax")),
      Missing_in_matrix = length(setdiff(metadata[sequencing_subset == "WGS", SampleID], setdiff(names(kraken_genus_rel), "Tax"))),
      Extra_in_matrix = length(setdiff(setdiff(names(kraken_genus_rel), "Tax"), metadata[sequencing_subset == "WGS", SampleID]))
    ),
    data.table(
      Dataset = "Kraken2 WGS species matrix vs stage5 WGS metadata",
      Expected_SRA_in_stage5 = stage5_wgs_sra,
      Matrix_SRA_after_zero_count_removal = length(setdiff(names(kraken_species_rel), "Tax")),
      Missing_in_matrix = length(setdiff(metadata[sequencing_subset == "WGS", SampleID], setdiff(names(kraken_species_rel), "Tax"))),
      Extra_in_matrix = length(setdiff(setdiff(names(kraken_species_rel), "Tax"), metadata[sequencing_subset == "WGS", SampleID]))
    ),
    data.table(
      Dataset = "mOTUs WGS matrix vs stage5 WGS metadata",
      Expected_SRA_in_stage5 = stage5_wgs_sra,
      Matrix_SRA_after_zero_count_removal = length(setdiff(names(motus_genus), "Tax")),
      Missing_in_matrix = length(setdiff(metadata[sequencing_subset == "WGS", SampleID], setdiff(names(motus_genus), "Tax"))),
      Extra_in_matrix = length(setdiff(setdiff(names(motus_genus), "Tax"), metadata[sequencing_subset == "WGS", SampleID]))
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

matrix_consistency[, Perfect_match := Missing_in_matrix == 0 & Extra_in_matrix == 0]
write_table_all(matrix_consistency, "Matrix_vs_stage5_consistency_checks_v7_metrics")

message("\nPublication-ready tables saved:")
message("- Table1_filtering_workflow_publication_with_notes.csv/tsv")
message("- Table1_filtering_workflow_publication_minimal.csv/tsv")
message("- Table2_final_analysis_dataset_summary_v7_metrics.csv/tsv")
message("- Matrix_vs_stage5_consistency_checks_v7_metrics.csv/tsv")

# ============================================================
# FINAL SUMMARY
# ============================================================

summary_lines <- c(
  "BioProject-level audit completed.",
  "Version: v8",
  paste0("Output directory: ", out_dir),
  paste0("Metadata: ", metadata_path),
  paste0("Initial metadata count check: ", initial_metadata_path),
  paste0("Curated center map: ", center_country_map_path),
  paste0("Curated center map exists: ", file.exists(center_country_map_path)),
  paste0("QIIME2 matrix: ", qiime2_amplicon_matrix_path),
  paste0("Kraken genus matrix: ", kraken_genus_matrix_path),
  paste0("Kraken species matrix: ", kraken_species_matrix_path),
  paste0("mOTUs matrix: ", motus_matrix_path),
  paste0("Threads used: ", n_cores),
  paste0("build_kraken_matrix: ", build_kraken_matrix),
  paste0("Detection threshold for prevalence: ", detection_threshold, " (", detection_threshold * 100, "%)"),
  paste0("Excluded host taxa: ", paste(excluded_host_taxa, collapse = ", ")),
  "",
  "Main output files:",
  "- README.md",
  "- RUN_SUMMARY.txt",
  "- CenterName_country_mapping_used.csv",
  "- BioProject_SRA_counts_metadata_summary.csv",
  "- Country_SRA_BioProject_summary.csv",
  "- Continent_SRA_BioProject_summary.csv",
  "- Country_Continent_SRA_BioProject_summary.csv",
  "- QIIME2_Amplicon_genus_BioProject_prevalence_summary.csv",
  "- Kraken2_WGS_genus_BioProject_prevalence_summary.csv",
  "- mOTUs_WGS_genus_BioProject_prevalence_summary.csv",
  "- QIIME2_vs_Kraken2_genus_concordance_raw.csv",
  "- QIIME2_vs_Kraken2_genus_concordance_filtered.csv",
  "- QIIME2_vs_Kraken2_concordance_qc.csv",
  "- QIIME2_vs_Kraken2_genus_concordance_stats.csv",
  "- QIIME2_vs_mOTUs_genus_concordance_raw.csv",
  "- QIIME2_vs_mOTUs_genus_concordance_filtered.csv",
  "- QIIME2_vs_mOTUs_concordance_qc.csv",
  "- QIIME2_vs_mOTUs_genus_concordance_stats.csv",
  "- Kraken2_vs_mOTUs_genus_concordance_raw.csv",
  "- Kraken2_vs_mOTUs_genus_concordance_filtered.csv",
  "- Kraken2_vs_mOTUs_concordance_qc.csv",
  "- Kraken2_vs_mOTUs_genus_concordance_stats.csv",
  "- Kraken2_WGS_species_species_richness_by_genus.csv",
  "- mOTUs_WGS_species_species_richness_by_genus.csv",
  "- Table1_filtering_workflow_publication_with_notes.csv",
  "- Table1_filtering_workflow_publication_minimal.csv",
  "- Table2_final_analysis_dataset_summary_v7_metrics.csv",
  "- Matrix_vs_stage5_consistency_checks_v7_metrics.csv"
)

writeLines(summary_lines, file.path(out_dir, "RUN_SUMMARY.txt"))
message(paste(summary_lines, collapse = "\n"))
