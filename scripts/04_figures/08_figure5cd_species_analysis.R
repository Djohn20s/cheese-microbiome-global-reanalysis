#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(ggrepel)
  library(grid)
})

options(stringsAsFactors = FALSE)
options(scipen = 999)

# ============================================================
# PATHS
# ============================================================

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)

if (length(file_arg) == 1) {
  script_file <- normalizePath(sub("^--file=", "", file_arg))
  repo_dir <- normalizePath(file.path(dirname(script_file), "..", ".."))
} else {
  repo_dir <- normalizePath(getwd())
}

input_dir <- file.path(
  repo_dir,
  "processed_data/bioproject_matrices"
)

out_dir <- file.path(
  repo_dir,
  "results/figures/figure5"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

kraken_file <- file.path(
  input_dir,
  "Kraken2_WGS_species_BioProject_taxa_wide.csv"
)

motus_file <- file.path(
  input_dir,
  "mOTUs_WGS_species_BioProject_taxa_wide.csv"
)

for (f in c(kraken_file, motus_file)) {
  if (!file.exists(f)) {
    stop("Missing required input: ", f)
  }
}

# ============================================================
# SETTINGS
# ============================================================

species_detection_threshold <- 0.0001

# Figure 5C
species_detection_bioprojects <- 3
species_detection_mean_abundance <- 0.0001
top_n_richness_genera <- 20

# Figure 5D
min_species_concordance_prevalence <- 0.10
min_species_concordance_abundance <- 0.0001
top_n_density_outliers <- 14
pseudo_count <- 1e-6

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

forced_upper_left_outlier_species <- c(
  "Streptococcus pasteurianus",
  "Pseudoalteromonas prydzensis",
  "Marinomonas polaris",
  "Vreelandella profundi"
)

forced_upper_right_core_species <- c(
  "Streptococcus thermophilus",
  "Lactococcus lactis",
  "Lactococcus cremoris"
)

meta_cols <- c(
  "BioProject",
  "CenterName_standardized",
  "CenterName",
  "Continent",
  "Country",
  "Platform"
)

# ============================================================
# HELPERS
# ============================================================

clean_species_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

extract_genus_from_species <- function(x) {
  sub(" .*", "", clean_species_label(x))
}

collapse_bioproject_wide <- function(dt) {

  meta_present <- intersect(meta_cols, names(dt))

  if (!"BioProject" %in% meta_present) {
    stop("BioProject column not detected.")
  }

  tax_cols <- setdiff(names(dt), meta_present)

  dt <- as.data.table(dt)

  dt[
    ,
    lapply(.SD, mean, na.rm = TRUE),
    by = BioProject,
    .SDcols = tax_cols
  ]
}

summarise_species_wide <- function(wide_dt, species_cols) {

  data.table(
    Species = species_cols,

    prevalence = sapply(
      species_cols,
      function(x) {
        mean(
          wide_dt[[x]] > species_detection_threshold,
          na.rm = TRUE
        )
      }
    ),

    mean_abundance = sapply(
      species_cols,
      function(x) mean(wide_dt[[x]], na.rm = TRUE)
    ),

    n_bioproject_detected = sapply(
      species_cols,
      function(x) {
        sum(
          wide_dt[[x]] > species_detection_threshold,
          na.rm = TRUE
        )
      }
    )
  )
}

make_richness <- function(wide_dt, method_name) {

  species_cols <- setdiff(names(wide_dt), "BioProject")

  species_summary <- data.table(
    Species = species_cols,
    Genus = extract_genus_from_species(species_cols),

    prevalence_bioproject = sapply(
      species_cols,
      function(x) mean(wide_dt[[x]] > 0, na.rm = TRUE)
    ),

    n_bioproject_detected = sapply(
      species_cols,
      function(x) sum(wide_dt[[x]] > 0, na.rm = TRUE)
    ),

    mean_abundance_bioproject = sapply(
      species_cols,
      function(x) mean(wide_dt[[x]], na.rm = TRUE)
    )
  )

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

  genus_richness[, Method := method_name]

  list(
    species_summary = species_summary,
    species_filtered = species_filtered,
    genus_richness = genus_richness
  )
}

# ============================================================
# LOAD + COLLAPSE BIOPROJECTS
# ============================================================

kraken_raw <- fread(kraken_file)
motus_raw <- fread(motus_file)

kraken <- collapse_bioproject_wide(kraken_raw)
motus <- collapse_bioproject_wide(motus_raw)

message("Kraken2 BioProjects: ", nrow(kraken))
message("mOTUs BioProjects: ", nrow(motus))

if (nrow(kraken) != 20) {
  stop("Expected 20 unique Kraken2 BioProjects.")
}

if (nrow(motus) != 19) {
  stop("Expected 19 unique mOTUs BioProjects.")
}

# ============================================================
# PANEL C — SPECIES RICHNESS BY GENUS
# ============================================================

kraken_richness_obj <- make_richness(
  kraken,
  "Kraken2"
)

motus_richness_obj <- make_richness(
  motus,
  "mOTUs"
)

fwrite(
  kraken_richness_obj$genus_richness,
  file.path(out_dir, "Figure5C_Kraken2_genus_richness.csv")
)

fwrite(
  motus_richness_obj$genus_richness,
  file.path(out_dir, "Figure5C_mOTUs_genus_richness.csv")
)

richness <- bind_rows(
  kraken_richness_obj$genus_richness,
  motus_richness_obj$genus_richness
) %>%
  as_tibble()

top_genera <- richness %>%
  filter(Genus != "Unknown") %>%
  group_by(Genus) %>%
  summarise(
    max_species = max(
      n_species_detected,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(max_species)) %>%
  slice_head(n = top_n_richness_genera) %>%
  pull(Genus)

richness_plot <- richness %>%
  filter(Genus %in% top_genera) %>%
  mutate(
    Genus = factor(
      Genus,
      levels = rev(top_genera)
    )
  )

fwrite(
  as.data.table(richness_plot),
  file.path(
    out_dir,
    "Figure5C_species_richness_by_genus_plot_data.csv"
  )
)

pC <- ggplot(
  richness_plot,
  aes(
    x = n_species_detected,
    y = Genus,
    fill = Method
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "C",
    x = "Number of detected species",
    y = NULL,
    fill = NULL
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0
    ),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.position = "bottom",
    plot.margin = margin(5, 10, 5, 5)
  )

ggsave(
  file.path(
    out_dir,
    "Figure5C_species_richness_by_genus.png"
  ),
  pC,
  width = 8,
  height = 7,
  dpi = 600
)

ggsave(
  file.path(
    out_dir,
    "Figure5C_species_richness_by_genus.pdf"
  ),
  pC,
  width = 8,
  height = 7
)

# ============================================================
# PANEL D — KRAKEN2 VS MOTUS SPECIES CONCORDANCE
# ============================================================

kraken_species <- setdiff(
  names(kraken),
  "BioProject"
)

motus_species <- setdiff(
  names(motus),
  "BioProject"
)

common_species <- intersect(
  kraken_species,
  motus_species
)

message(
  "Common Kraken2+mOTUs species: ",
  length(common_species)
)

if (length(common_species) != 448) {
  stop(
    "Expected 448 common species; found ",
    length(common_species)
  )
}

kraken_summary <- summarise_species_wide(
  kraken,
  common_species
)

motus_summary <- summarise_species_wide(
  motus,
  common_species
)

species_summary <- merge(
  kraken_summary,
  motus_summary,
  by = "Species",
  suffixes = c("_Kraken2", "_mOTUs")
) %>%
  as_tibble() %>%
  mutate(
    Genus = extract_genus_from_species(Species),

    prevalence_mean =
      (prevalence_Kraken2 +
       prevalence_mOTUs) / 2,

    mean_abundance_mean =
      (mean_abundance_Kraken2 +
       mean_abundance_mOTUs) / 2
  )

species_concordance_plot <- species_summary %>%
  filter(
    prevalence_mean >=
      min_species_concordance_prevalence,

    mean_abundance_mean >=
      min_species_concordance_abundance
  ) %>%
  mutate(
    mean_abundance_Kraken2_safe =
      pmax(mean_abundance_Kraken2, 0) +
      pseudo_count,

    mean_abundance_mOTUs_safe =
      pmax(mean_abundance_mOTUs, 0) +
      pseudo_count,

    log10_Kraken2 =
      log10(mean_abundance_Kraken2_safe),

    log10_mOTUs =
      log10(mean_abundance_mOTUs_safe),

    log10_ratio_Kraken2_vs_mOTUs =
      log10_Kraken2 - log10_mOTUs,

    abs_log10_ratio =
      abs(log10_ratio_Kraken2_vs_mOTUs),

    fold_change_major_method =
      10^abs_log10_ratio,

    higher_method = case_when(
      log10_ratio_Kraken2_vs_mOTUs > 0 ~
        "Kraken2 higher",

      log10_ratio_Kraken2_vs_mOTUs < 0 ~
        "mOTUs higher",

      TRUE ~ "Equal"
    )
  ) %>%
  arrange(desc(mean_abundance_mean))

message(
  "Species retained for Figure 5D: ",
  nrow(species_concordance_plot)
)

if (nrow(species_concordance_plot) != 103) {
  stop(
    "Expected 103 species after filtering; found ",
    nrow(species_concordance_plot)
  )
}

species_concordance_test <- cor.test(
  species_concordance_plot$mean_abundance_Kraken2_safe,
  species_concordance_plot$mean_abundance_mOTUs_safe,
  method = "spearman",
  exact = FALSE
)

rho <- unname(
  species_concordance_test$estimate
)

p_value <- species_concordance_test$p.value

message(
  "Spearman rho: ",
  format(rho, digits = 8)
)

message(
  "P-value: ",
  format(p_value, scientific = TRUE)
)

if (abs(rho - 0.5338598) > 1e-6) {
  stop(
    "Unexpected Spearman rho: ",
    rho
  )
}

species_concordance_label <- paste0(
  "Spearman rho = ",
  sprintf("%.2f", rho),
  "\n",
  "p = ",
  formatC(
    p_value,
    format = "e",
    digits = 1
  ),
  "\n",
  "n species = ",
  nrow(species_concordance_plot)
)

annotation_x <- min(
  species_concordance_plot$log10_Kraken2,
  na.rm = TRUE
)

annotation_y <- max(
  species_concordance_plot$log10_mOTUs,
  na.rm = TRUE
)

# ============================================================
# OUTLIER LABELS — SAME RULES AS FINAL FIGURE
# ============================================================

density_outlier_auto <- species_concordance_plot %>%
  filter(
    !is.na(Species),
    Species != ""
  ) %>%
  arrange(
    desc(abs_log10_ratio),
    desc(prevalence_mean),
    desc(mean_abundance_mean)
  ) %>%
  slice_head(n = top_n_density_outliers)

density_outlier_forced_upper_right <-
  species_concordance_plot %>%
  filter(
    Species %in%
      forced_upper_right_core_species
  ) %>%
  mutate(
    forced_group =
      "Dominant concordant core"
  )

density_outlier_forced_upper_left <-
  species_concordance_plot %>%
  filter(
    Species %in%
      forced_upper_left_outlier_species
  ) %>%
  mutate(
    forced_group =
      "mOTUs-higher outlier"
  )

density_outlier_auto <-
  density_outlier_auto %>%
  mutate(
    forced_group =
      "Largest method-discordant outlier"
  )

density_outlier_species <- bind_rows(
  density_outlier_forced_upper_right,
  density_outlier_forced_upper_left,
  density_outlier_auto
) %>%
  distinct(Species, .keep_all = TRUE) %>%
  mutate(
    forced_priority = case_when(
      Species %in%
        forced_upper_right_core_species ~ 1L,

      Species %in%
        forced_upper_left_outlier_species ~ 2L,

      TRUE ~ 3L
    ),

    forced_order = case_when(
      Species ==
        "Streptococcus thermophilus" ~ 1L,

      Species ==
        "Lactococcus lactis" ~ 2L,

      Species ==
        "Lactococcus cremoris" ~ 3L,

      Species ==
        "Streptococcus pasteurianus" ~ 4L,

      Species ==
        "Pseudoalteromonas prydzensis" ~ 5L,

      Species ==
        "Marinomonas polaris" ~ 6L,

      TRUE ~ 100L
    )
  ) %>%
  arrange(
    forced_priority,
    forced_order,
    desc(abs_log10_ratio),
    desc(prevalence_mean),
    desc(mean_abundance_mean)
  ) %>%
  slice_head(
    n = top_n_density_outliers
  ) %>%
  mutate(
    Outlier_ID = seq_len(n()),
    Outlier_label =
      as.character(as.roman(Outlier_ID)),

    Species_chr =
      as.character(Species),

    Species = factor(
      Species_chr,
      levels = Species_chr
    ),

    legend_label = paste(
      sprintf("%-4s", Outlier_label),
      Species_chr
    )
  )

density_outlier_species$legend_label <-
  factor(
    density_outlier_species$legend_label,
    levels =
      density_outlier_species$legend_label
  )

# ============================================================
# EXPORT PANEL D SOURCE DATA
# ============================================================

fwrite(
  as.data.table(species_summary),
  file.path(
    out_dir,
    "Figure5D_common_species_summary.csv"
  )
)

fwrite(
  as.data.table(species_concordance_plot),
  file.path(
    out_dir,
    "Figure5D_species_Kraken2_mOTUs_concordance_plot_data.csv"
  )
)

fwrite(
  as.data.table(density_outlier_species),
  file.path(
    out_dir,
    "Figure5D_top14_density_outlier_species.csv"
  )
)

fwrite(
  data.table(
    rho = rho,
    p_value = p_value,
    n_species = nrow(species_concordance_plot),
    common_species_before_filtering =
      length(common_species),
    prevalence_threshold =
      min_species_concordance_prevalence,
    mean_abundance_threshold =
      min_species_concordance_abundance
  ),
  file.path(
    out_dir,
    "Figure5D_species_Kraken2_mOTUs_concordance_stats.csv"
  )
)

# ============================================================
# FIGURE 5D
# ============================================================

pD <- ggplot(
  species_concordance_plot,
  aes(
    x = log10_Kraken2,
    y = log10_mOTUs
  )
) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    alpha = 0.55,
    contour = TRUE,
    bins = 8
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_point(
    aes(size = prevalence_mean),
    alpha = 0.35,
    colour = "grey35"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.7,
    colour = "black"
  ) +
  geom_point(
    data = density_outlier_species,
    aes(colour = legend_label),
    size = 3.4,
    alpha = 0.98,
    show.legend = TRUE
  ) +
  ggrepel::geom_label_repel(
    data = density_outlier_species,
    aes(
      label = Outlier_label,
      colour = legend_label
    ),
    size = 3.2,
    fontface = "bold",
    fill = "white",
    label.size = 0.2,
    label.padding = unit(
      0.15,
      "lines"
    ),
    max.overlaps = Inf,
    box.padding = 0.45,
    point.padding = 0.25,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = annotation_x,
    y = annotation_y,
    label = species_concordance_label,
    hjust = 0,
    vjust = 1,
    size = 3.0
  ) +
  scale_size_continuous(
    labels = percent_format(
      accuracy = 1
    ),
    range = c(1.2, 4)
  ) +
  guides(
    colour = guide_legend(
      title = NULL,
      override.aes =
        list(size = 3.5, alpha = 1)
    ),
    fill =
      guide_colourbar(
        title = "Density"
      ),
    size =
      guide_legend(
        title = "Mean prevalence"
      )
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "D",
    x = expression(
      log[10](
        "mean relative abundance, Kraken2"
      )
    ),
    y = expression(
      log[10](
        "mean relative abundance, mOTUs"
      )
    ),
    fill = "Density",
    size = "Mean prevalence",
    colour = NULL
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0
    ),
    panel.grid.minor =
      element_blank(),
    legend.position = "right",
    legend.text =
      element_text(size = 7.2),
    legend.title =
      element_text(size = 8.5),
    legend.box.margin =
      margin(l = 18, r = 8),
    plot.margin =
      margin(
        t = 8,
        r = 140,
        b = 8,
        l = 8
      )
  )

ggsave(
  file.path(
    out_dir,
    "Figure5D_species_Kraken2_mOTUs_concordance.png"
  ),
  pD,
  width = 8.5,
  height = 7,
  dpi = 600
)

ggsave(
  file.path(
    out_dir,
    "Figure5D_species_Kraken2_mOTUs_concordance.pdf"
  ),
  pD,
  width = 8.5,
  height = 7
)

# ============================================================
# SUMMARY
# ============================================================

message("")
message("Figure 5C/D completed.")
message(
  "Figure 5C richness filter: >= ",
  species_detection_bioprojects,
  " BioProjects and mean abundance >= ",
  species_detection_mean_abundance
)

message(
  "Figure 5D: ",
  length(common_species),
  " common species -> ",
  nrow(species_concordance_plot),
  " retained species"
)

message(
  "Figure 5D Spearman rho = ",
  format(rho, digits = 8),
  "; p = ",
  format(p_value, scientific = TRUE)
)

message("Output directory: ", out_dir)
