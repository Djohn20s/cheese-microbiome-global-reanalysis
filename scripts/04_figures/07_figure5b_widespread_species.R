#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

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

source_file <- file.path(
  repo_dir,
  "processed_data/core_species/Figure5B_widespread_species_exact_counts.tsv"
)

out_dir <- file.path(
  repo_dir,
  "results/figures/figure5"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_file)) {
  stop("Missing Figure 5B source data: ", source_file)
}

# ============================================================
# LOAD FINAL AUDITED SPECIES TABLE
# ============================================================

widespread <- fread(source_file)

required_cols <- c(
  "Species",
  "Kraken2_detected_n",
  "Kraken2_denominator",
  "Kraken2_prevalence",
  "mOTUs_detected_n",
  "mOTUs_denominator",
  "mOTUs_prevalence",
  "Kraken2_mean_abundance",
  "mOTUs_mean_abundance",
  "prevalence_min",
  "widespread_class",
  "LOO_iterations_core",
  "LOO_iterations_total",
  "LOO_stability"
)

missing_cols <- setdiff(required_cols, names(widespread))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ============================================================
# VALIDATION
# ============================================================

if (nrow(widespread) != 12) {
  stop(
    "Expected 12 widespread species, found ",
    nrow(widespread)
  )
}

if (!all(widespread$Kraken2_denominator == 20)) {
  stop("Unexpected Kraken2 denominator.")
}

if (!all(widespread$mOTUs_denominator == 19)) {
  stop("Unexpected mOTUs denominator.")
}

if (!all(
  widespread$Kraken2_prevalence >= 0.50 &
  widespread$mOTUs_prevalence >= 0.50
)) {
  stop("At least one species does not satisfy the >=50% prevalence criterion.")
}

if (!all(
  widespread$Kraken2_mean_abundance >= 0.0001 &
  widespread$mOTUs_mean_abundance >= 0.0001
)) {
  stop("At least one species does not satisfy the >=0.01% mean-abundance criterion.")
}

setorder(
  widespread,
  -prevalence_min,
  Species
)

message("Widespread species: ", nrow(widespread))
message(
  "Highly widespread species: ",
  sum(widespread$widespread_class == "Highly widespread (≥75% in both)")
)
message(
  "Stable in all leave-one-BioProject-out iterations: ",
  sum(widespread$LOO_stability == 1)
)

# ============================================================
# LONG FORMAT
# ============================================================

plot_long <- bind_rows(
  widespread %>%
    transmute(
      Species,
      Pipeline = "Kraken2",
      prevalence = Kraken2_prevalence,
      detected_n = Kraken2_detected_n,
      denominator = Kraken2_denominator
    ),

  widespread %>%
    transmute(
      Species,
      Pipeline = "mOTUs",
      prevalence = mOTUs_prevalence,
      detected_n = mOTUs_detected_n,
      denominator = mOTUs_denominator
    )
)

species_order <- widespread$Species

plot_long$Species <- factor(
  plot_long$Species,
  levels = rev(species_order)
)

plot_long$Pipeline <- factor(
  plot_long$Pipeline,
  levels = c("Kraken2", "mOTUs")
)

plot_long <- plot_long %>%
  mutate(
    count_label = paste0(detected_n, "/", denominator)
  ) %>%
  group_by(Species) %>%
  mutate(
    label_side = case_when(
      prevalence == min(prevalence) &
        prevalence == max(prevalence) &
        Pipeline == "Kraken2" ~ "left",

      prevalence == min(prevalence) &
        prevalence == max(prevalence) &
        Pipeline == "mOTUs" ~ "right",

      prevalence == min(prevalence) ~ "left",

      TRUE ~ "right"
    ),

    label_x = if_else(
      label_side == "left",
      prevalence * 100 - 1.0,
      prevalence * 100 + 1.0
    ),

    label_hjust = if_else(
      label_side == "left",
      1,
      0
    )
  ) %>%
  ungroup()

fwrite(
  as.data.table(plot_long),
  file.path(out_dir, "Figure5B_plot_data.tsv"),
  sep = "\t"
)

# ============================================================
# FIGURE 5B
# ============================================================

p5b <- ggplot(
  plot_long,
  aes(
    x = prevalence * 100,
    y = Species,
    shape = Pipeline
  )
) +
  geom_vline(
    xintercept = 50,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_line(
    aes(group = Species),
    linewidth = 0.45
  ) +
  geom_point(
    size = 3.2,
    stroke = 0.8
  ) +
  geom_text(
    aes(
      x = label_x,
      label = count_label,
      hjust = label_hjust
    ),
    size = 3.0,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    limits = c(40, 110),
    breaks = c(50, 60, 70, 80, 90, 100),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "B",
    x = "BioProject prevalence",
    y = NULL,
    shape = NULL,
    caption = paste0(
      "Labels indicate detected BioProjects / total BioProjects. ",
      "Kraken2 denominator = 20; mOTUs denominator = 19."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),

    axis.text.y = element_text(
      face = "italic",
      size = 9
    ),

    plot.title = element_text(
      face = "bold",
      size = 12
    ),

    legend.position = "top",

    plot.margin = margin(
      8, 35, 8, 8
    )
  )

ggsave(
  file.path(
    out_dir,
    "Figure5B_REVISED_widespread_species.png"
  ),
  p5b,
  width = 8.3,
  height = 6.5,
  dpi = 600
)

ggsave(
  file.path(
    out_dir,
    "Figure5B_REVISED_widespread_species.pdf"
  ),
  p5b,
  width = 8.3,
  height = 6.5
)

message("Figure 5B completed.")
message("Input: ", source_file)
message("Output directory: ", out_dir)
