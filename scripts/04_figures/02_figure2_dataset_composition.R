#!/usr/bin/env Rscript

# ============================================================
# Figure 2. Composition and representativeness of final curated
# public cheese microbiome datasets
#
# Panels:
# A. SRA records per BioProject
# B. SRA records per year by library strategy
# C. BioProjects per year
# D. Sequencing platforms by SRA records
#
# Version: v6_final8910
# Dataset: final analyzed Stage 5 dataset only
#          runinfo_WGS.csv + runinfo_AMPLICON.csv
# Expected total: 8,910 SRA records
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(scales)
  library(stringi)
})

# ============================================================
# Paths
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
wgs_file <- file.path(
  project_dir,
  "metadata/curation/runinfo_WGS.csv"
)

amp_file <- file.path(
  project_dir,
  "metadata/curation/runinfo_AMPLICON.csv"
)

center_country_map_path <- file.path(
  project_dir,
  "metadata/curation/center_country_map_curated.csv"
)

out_dir <- file.path(
  project_dir,
  "results/figures/figure2"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper functions
# ============================================================

normalize_center <- function(x) {
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

standardize_unknown <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  dplyr::if_else(is.na(x) | x == "", "Unknown", x)
}

parse_release_date <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(lubridate::ymd_hms(x, quiet = TRUE))
  out_date <- as.Date(out)
  missing <- is.na(out_date)
  if (any(missing)) {
    out_date[missing] <- suppressWarnings(lubridate::ymd(x[missing], quiet = TRUE))
  }
  out_date
}

# ============================================================
# Load final curated runinfo: WGS + AMPLICON only
# ============================================================

if (!file.exists(wgs_file)) stop("Missing WGS file: ", wgs_file)
if (!file.exists(amp_file)) stop("Missing AMPLICON file: ", amp_file)

message("Loading final WGS runinfo: ", wgs_file)
message("Loading final AMPLICON runinfo: ", amp_file)

wgs_runinfo <- readr::read_csv(wgs_file, show_col_types = FALSE) %>%
  mutate(LibraryStrategy = "WGS")

amp_runinfo <- readr::read_csv(amp_file, show_col_types = FALSE) %>%
  mutate(LibraryStrategy = "Amplicon")

runinfo <- bind_rows(wgs_runinfo, amp_runinfo) %>%
  mutate(
    BioProject = stringr::str_squish(as.character(BioProject)),
    CenterName = if ("CenterName" %in% names(.)) {
      stringr::str_squish(as.character(CenterName))
    } else {
      NA_character_
    },
    Platform = if ("Platform" %in% names(.)) {
      stringr::str_squish(as.character(Platform))
    } else {
      NA_character_
    },
    LibraryStrategy = factor(
      LibraryStrategy,
      levels = c("Amplicon", "WGS")
    ),
    ReleaseDate_parsed = if ("ReleaseDate" %in% names(.)) {
      parse_release_date(ReleaseDate)
    } else {
      as.Date(NA)
    },
    Year = lubridate::year(ReleaseDate_parsed),
    Continent_original = if ("Continent" %in% names(.)) {
      stringr::str_squish(as.character(Continent))
    } else {
      NA_character_
    },
    Country_original = if ("Country" %in% names(.)) {
      stringr::str_squish(as.character(Country))
    } else {
      NA_character_
    },
    CenterName_norm = normalize_center(CenterName)
  ) %>%
  filter(
    !is.na(BioProject),
    BioProject != ""
  )

message("Final Figure 2 dataset SRA records: ", nrow(runinfo))
message("Final Figure 2 dataset BioProjects: ", n_distinct(runinfo$BioProject))
message("Final Figure 2 amplicon SRA: ", sum(runinfo$LibraryStrategy == "Amplicon", na.rm = TRUE))
message("Final Figure 2 WGS SRA: ", sum(runinfo$LibraryStrategy == "WGS", na.rm = TRUE))

# ============================================================
# Load curated CenterName map
# ============================================================

if (!file.exists(center_country_map_path)) {
  stop(
    "Curated center-country map not found: ", center_country_map_path,
    "\nCreate it first as center_country_map_curated.csv with columns: CenterName, CenterName_standardized, Country, Continent"
  )
}

message("Loading curated center map: ", center_country_map_path)

center_map <- readr::read_csv(
  center_country_map_path,
  show_col_types = FALSE
) %>%
  mutate(
    CenterName = stringr::str_squish(as.character(CenterName)),
    CenterName_norm = normalize_center(CenterName),
    CenterName_standardized = stringr::str_squish(as.character(CenterName_standardized)),
    Country_curated = stringr::str_squish(as.character(Country)),
    Continent_curated = stringr::str_squish(as.character(Continent))
  ) %>%
  distinct(
    CenterName_norm,
    .keep_all = TRUE
  ) %>%
  select(
    CenterName_norm,
    CenterName_standardized,
    Country_curated,
    Continent_curated
  )

# ============================================================
# Join curated metadata
# ============================================================

runinfo <- runinfo %>%
  left_join(
    center_map,
    by = "CenterName_norm"
  ) %>%
  mutate(
    Country = dplyr::coalesce(
      standardize_unknown(Country_curated),
      standardize_unknown(Country_original)
    ),
    Continent = dplyr::coalesce(
      standardize_unknown(Continent_curated),
      standardize_unknown(Continent_original)
    ),
    Country = standardize_unknown(Country),
    Continent = standardize_unknown(Continent),
    Platform = standardize_unknown(Platform),
    CenterName_standardized = if_else(
      is.na(CenterName_standardized) |
        CenterName_standardized == "",
      CenterName_norm,
      CenterName_standardized
    )
  )

# ============================================================
# Diagnostics
# ============================================================

mapping_diagnostics <- runinfo %>%
  distinct(
    CenterName,
    CenterName_norm,
    CenterName_standardized,
    Country,
    Continent
  ) %>%
  arrange(
    Country,
    CenterName_standardized,
    CenterName
  )

unmapped_centers <- mapping_diagnostics %>%
  filter(
    Country == "Unknown" | Continent == "Unknown"
  )

readr::write_csv(
  mapping_diagnostics,
  file.path(out_dir, "Figure2_CenterName_country_mapping_used.csv")
)

readr::write_csv(
  unmapped_centers,
  file.path(out_dir, "Figure2_unmapped_centers.csv")
)

message("Unique centers: ", nrow(mapping_diagnostics))
message("Unmapped centers: ", nrow(unmapped_centers))

# ============================================================
# Summary tables
# ============================================================

# A. SRA per BioProject
# Keep BioProject as the unit, but annotate each BioProject by its dominant
# library strategy so Panel A uses the same color language as Panel B.
bioproject_strategy <- runinfo %>%
  count(BioProject, LibraryStrategy, name = "n_strategy_sra") %>%
  group_by(BioProject) %>%
  slice_max(n_strategy_sra, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    BioProject,
    dominant_library_strategy = as.character(LibraryStrategy)
  )

bioproject_counts <- runinfo %>%
  count(BioProject, name = "n_sra") %>%
  left_join(bioproject_strategy, by = "BioProject") %>%
  arrange(desc(n_sra)) %>%
  mutate(
    rank = row_number(),
    cumulative_sra = cumsum(n_sra),
    cumulative_fraction = cumulative_sra / sum(n_sra),
    cumulative_percent = 100 * cumulative_fraction,
    dominant_library_strategy = factor(
      dominant_library_strategy,
      levels = c("Amplicon", "WGS")
    )
  )

# Metadata summary per BioProject, useful for supplementary tables
bioproject_metadata <- runinfo %>%
  group_by(BioProject) %>%
  summarise(
    n_sra = n(),
    Country = first(na.omit(Country)),
    Continent = first(na.omit(Continent)),
    Platform = paste(sort(unique(na.omit(Platform))), collapse = "; "),
    LibraryStrategy = paste(sort(unique(na.omit(as.character(LibraryStrategy)))), collapse = "; "),
    CenterName_standardized = first(na.omit(CenterName_standardized)),
    first_release = suppressWarnings(min(ReleaseDate_parsed, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sra))

# B. SRA per year + library strategy
sra_strategy_per_year <- runinfo %>%
  filter(!is.na(Year)) %>%
  count(Year, LibraryStrategy, name = "n_sra") %>%
  arrange(Year, LibraryStrategy)

sra_per_year <- runinfo %>%
  filter(!is.na(Year)) %>%
  count(Year, name = "n_sra") %>%
  arrange(Year)

# C. BioProjects per year
# First release year of each BioProject is used.
bioproject_year <- runinfo %>%
  filter(!is.na(Year)) %>%
  group_by(BioProject) %>%
  summarise(
    Year = min(Year, na.rm = TRUE),
    .groups = "drop"
  )

bioproject_per_year <- bioproject_year %>%
  count(Year, name = "n_bioprojects") %>%
  arrange(Year)

# D. SRA per sequencing platform
platform_sra <- runinfo %>%
  count(Platform, name = "n_sra") %>%
  arrange(desc(n_sra))

# Additional useful audit summaries, not directly plotted
country_sra_counts <- runinfo %>%
  count(Country, name = "n_sra") %>%
  arrange(desc(n_sra))

continent_sra_counts <- runinfo %>%
  count(Continent, name = "n_sra") %>%
  arrange(desc(n_sra))

strategy_sra_counts <- runinfo %>%
  count(LibraryStrategy, name = "n_sra") %>%
  arrange(desc(n_sra))

country_bioproject_counts <- runinfo %>%
  group_by(Country, Continent) %>%
  summarise(
    n_sra = n(),
    n_bioprojects = n_distinct(BioProject),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sra), desc(n_bioprojects))

continent_bioproject_counts <- runinfo %>%
  group_by(Continent) %>%
  summarise(
    n_sra = n(),
    n_bioprojects = n_distinct(BioProject),
    n_countries = n_distinct(Country),
    Countries = paste(sort(unique(Country)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sra), desc(n_bioprojects))

# ============================================================
# Save summary tables
# ============================================================

readr::write_csv(
  bioproject_counts,
  file.path(out_dir, "Figure2_A_BioProject_SRA_counts.csv")
)

readr::write_csv(
  bioproject_metadata,
  file.path(out_dir, "Figure2_BioProject_metadata_curated.csv")
)

readr::write_csv(
  sra_strategy_per_year,
  file.path(out_dir, "Figure2_B_SRA_per_year_by_library_strategy.csv")
)

readr::write_csv(
  sra_per_year,
  file.path(out_dir, "Figure2_B_SRA_per_year_total.csv")
)

readr::write_csv(
  bioproject_per_year,
  file.path(out_dir, "Figure2_C_BioProjects_per_year.csv")
)

readr::write_csv(
  platform_sra,
  file.path(out_dir, "Figure2_D_SRA_per_platform.csv")
)

readr::write_csv(
  country_sra_counts,
  file.path(out_dir, "Supplementary_SRA_per_country.csv")
)

readr::write_csv(
  continent_sra_counts,
  file.path(out_dir, "Supplementary_SRA_per_continent.csv")
)

readr::write_csv(
  strategy_sra_counts,
  file.path(out_dir, "Supplementary_SRA_per_library_strategy.csv")
)

readr::write_csv(
  country_bioproject_counts,
  file.path(out_dir, "Supplementary_Country_SRA_BioProject_summary.csv")
)

readr::write_csv(
  continent_bioproject_counts,
  file.path(out_dir, "Supplementary_Continent_SRA_BioProject_summary.csv")
)

# ============================================================
# Shared theme
# ============================================================

theme_fig2 <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    legend.position = "bottom"
  )

# Consistent color palette for Figure 2
# Keep colors stable across panels and exports.
col_amplicon <- "#1F78B4"
col_wgs      <- "#E31A1C"
col_bioproj  <- "#33A02C"
col_platform <- "#6A3D9A"
col_cumulative <- "#111111"

library_strategy_cols <- c(
  "Amplicon" = col_amplicon,
  "WGS" = col_wgs
)

platform_cols <- c(
  "ILLUMINA" = "#1F78B4",
  "OXFORD_NANOPORE" = "#E31A1C",
  "PACBIO_SMRT" = "#33A02C",
  "LS454" = "#FF7F00",
  "ION_TORRENT" = "#6A3D9A",
  "ABI_SOLID" = "#B15928",
  "ELEMENT" = "#A6CEE3",
  "Unknown" = "#999999"
)

# ============================================================
# Panel A: SRA records per BioProject
# ============================================================

max_sra_panel_a <- max(bioproject_counts$n_sra, na.rm = TRUE)

pA <- ggplot(
  bioproject_counts,
  aes(x = rank, y = n_sra)
) +
  geom_col(aes(fill = dominant_library_strategy), alpha = 0.95) +
  geom_line(
    aes(y = cumulative_fraction * max_sra_panel_a, group = 1),
    color = col_cumulative,
    linewidth = 0.8
  ) +
  geom_point(
    aes(y = cumulative_fraction * max_sra_panel_a),
    color = col_cumulative,
    size = 0.9,
    alpha = 0.8
  ) +
  scale_fill_manual(values = library_strategy_cols, drop = FALSE) +
  scale_y_continuous(
    labels = comma,
    sec.axis = sec_axis(
      ~ . / max_sra_panel_a * 100,
      name = "Cumulative SRA (%)",
      labels = function(x) paste0(round(x), "%")
    )
  ) +
  labs(
        x = "BioProjects ranked by number of SRA records",
    y = "Number of SRA records",
    fill = "Dominant library strategy"
  ) +
  theme_fig2 +
  theme(
    legend.position = "bottom"
  )

# ============================================================
# Panel B: SRA records per year + library strategy
# ============================================================

pB <- ggplot(
  sra_strategy_per_year,
  aes(
    x = Year,
    y = n_sra,
    fill = LibraryStrategy
  )
) +
  geom_col(position = "stack") +
  scale_fill_manual(values = library_strategy_cols, drop = FALSE) +
  scale_y_continuous(labels = comma) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 8)
  ) +
  labs(
        x = NULL,
    y = "Number of SRA records",
    fill = "Library strategy"
  ) +
  theme_fig2

# ============================================================
# Panel C: BioProjects per year
# ============================================================

pC <- ggplot(
  bioproject_per_year,
  aes(x = Year, y = n_bioprojects)
) +
  geom_col(fill = "#59A14F") +
  scale_y_continuous(
    labels = comma,
    breaks = scales::pretty_breaks(n = 5)
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 8)
  ) +
  labs(
        x = "Year",
    y = "Number of BioProjects"
  ) +
  theme_fig2 +
  theme(
    legend.position = "none"
  )

# ============================================================
# Panel D: Sequencing platforms
# ============================================================

platform_plot <- platform_sra %>%
  filter(Platform != "Unknown") %>%
  mutate(
    Platform = forcats::fct_reorder(Platform, n_sra)
  )

pD <- ggplot(
  platform_plot,
  aes(
    x = n_sra,
    y = Platform,
    fill = Platform
  )
) +
  geom_col() +
  scale_fill_manual(values = platform_cols, drop = FALSE) +
  scale_x_log10(
    labels = comma,
    breaks = scales::trans_breaks("log10", function(x) 10^x)
  ) +
  labs(
        x = "Number of SRA records, log10 scale",
    y = NULL
  ) +
  theme_fig2 +
  theme(
    legend.position = "none"
  )

# ============================================================
# Combine
# ============================================================

figure2 <- (pA | pB) / (pC | pD)

# ============================================================
# Save
# ============================================================

ggsave(
  filename = file.path(out_dir, "Figure2_dataset_composition_curated_v6_final8910.png"),
  plot = figure2,
  width = 14,
  height = 10,
  dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Figure2_dataset_composition_curated_v6_final8910.pdf"),
  plot = figure2,
  width = 14,
  height = 10
)

# ============================================================
# README / summary
# ============================================================

summary_lines <- c(
  "Figure 2 dataset composition completed.",
  "Version: v6_final8910",
  "Dataset: final analyzed Stage 5 WGS + AMPLICON only.",
  paste0("Output directory: ", out_dir),
  paste0("WGS runinfo: ", wgs_file),
  paste0("AMPLICON runinfo: ", amp_file),
  paste0("Curated center map: ", center_country_map_path),
  paste0("Total SRA records: ", nrow(runinfo)),
  paste0("Total BioProjects: ", n_distinct(runinfo$BioProject)),
  paste0("Amplicon SRA records: ", sum(runinfo$LibraryStrategy == "Amplicon", na.rm = TRUE)),
  paste0("WGS SRA records: ", sum(runinfo$LibraryStrategy == "WGS", na.rm = TRUE)),
  paste0("Amplicon BioProjects: ", n_distinct(runinfo$BioProject[runinfo$LibraryStrategy == "Amplicon"])),
  paste0("WGS BioProjects: ", n_distinct(runinfo$BioProject[runinfo$LibraryStrategy == "WGS"])),
  paste0("Countries represented: ", n_distinct(runinfo$Country[runinfo$Country != "Unknown"])),
  paste0("Continents represented: ", n_distinct(runinfo$Continent[runinfo$Continent != "Unknown"])),
  paste0("Library strategies represented: ", paste(sort(unique(as.character(runinfo$LibraryStrategy))), collapse = ", ")),
  paste0("Platforms represented: ", paste(sort(unique(runinfo$Platform)), collapse = ", ")),
  paste0("Unique centers: ", nrow(mapping_diagnostics)),
  paste0("Unmapped centers: ", nrow(unmapped_centers)),
  "",
  "Figure panels:",
  "A. SRA records per BioProject with cumulative contribution",
  "B. SRA records per year by library strategy",
  "C. BioProjects per year with temporal trend",
  "D. Sequencing platforms by SRA records, log10 scale",
  "",
  "Main output files:",
  "- Figure2_dataset_composition_curated_v6_final8910.png",
  "- Figure2_dataset_composition_curated_v6_final8910.pdf",
  "- Figure2_A_BioProject_SRA_counts.csv",
  "- Figure2_B_SRA_per_year_by_library_strategy.csv",
  "- Figure2_B_SRA_per_year_total.csv",
  "- Figure2_C_BioProjects_per_year.csv",
  "- Figure2_D_SRA_per_platform.csv",
  "- Figure2_BioProject_metadata_curated.csv",
  "- Figure2_CenterName_country_mapping_used.csv",
  "- Figure2_unmapped_centers.csv",
  "",
  "Supplementary output files:",
  "- Supplementary_SRA_per_country.csv",
  "- Supplementary_SRA_per_continent.csv",
  "- Supplementary_SRA_per_library_strategy.csv",
  "- Supplementary_Country_SRA_BioProject_summary.csv",
  "- Supplementary_Continent_SRA_BioProject_summary.csv"
)

writeLines(
  summary_lines,
  file.path(out_dir, "README_Figure2.txt")
)

message(paste(summary_lines, collapse = "\n"))

figure2
