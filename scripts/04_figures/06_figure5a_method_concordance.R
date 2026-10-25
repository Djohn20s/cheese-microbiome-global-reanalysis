#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(scales)
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

source_file <- file.path(
  project_dir,
  "results/figures/figure3/Figure3_pipeline_concordance_statistics_all_genera.csv"
)

out_dir <- file.path(
  project_dir,
  "results/figures/figure5"
)

supp_dir <- file.path(
  project_dir,
  "supplementary"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_file)) {
  stop("Missing canonical Figure 3 concordance statistics: ", source_file)
}

# ============================================================
# LOAD CANONICAL CONCORDANCE RESULTS
# ============================================================

src <- fread(source_file)

required_cols <- c(
  "comparison",
  "n_common_genera",
  "n_filtered_genera",
  "prevalence_spearman_rho_filtered",
  "prevalence_p_value_filtered",
  "abundance_spearman_rho_filtered",
  "abundance_p_value_filtered"
)

missing_cols <- setdiff(required_cols, names(src))

if (length(missing_cols) > 0) {
  stop(
    "Missing required canonical columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Convert one row per comparison into two rows:
# prevalence and mean relative abundance.
dt <- rbindlist(
  list(
    src[, .(
      comparison,
      metric = "BioProject prevalence",
      rho = prevalence_spearman_rho_filtered,
      p_value = prevalence_p_value_filtered,
      n_common_genera_raw = n_common_genera,
      n_common_genera_filtered = n_filtered_genera
    )],

    src[, .(
      comparison,
      metric = "Mean relative abundance",
      rho = abundance_spearman_rho_filtered,
      p_value = abundance_p_value_filtered,
      n_common_genera_raw = n_common_genera,
      n_common_genera_filtered = n_filtered_genera
    )]
  ),
  use.names = TRUE
)

dt[, comparison := gsub("_", " ", as.character(comparison))]

# Explicit methodological scope.
dt[, comparison_scope := fifelse(
  comparison == "Kraken2 vs mOTUs",
  "Paired comparison: 19 common WGS BioProjects",
  "Aggregate cross-collection comparison"
)]

# S3 uses exactly the same canonical values as Figure 5A.
fwrite(
  dt,
  file.path(
    supp_dir,
    "Supplementary_Table_S3_pairwise_genus_concordance.csv"
  )
)

# ============================================================
# CLASSIFY COMPARISONS
# ============================================================

dt[, comparison_type := fifelse(
  comparison == "Kraken2 vs mOTUs",
  "Direct WGS cross-pipeline comparison",
  "Aggregate cross-collection comparison"
)]

# Order intentionally puts paired WGS comparison first
dt[, comparison := factor(
  comparison,
  levels = c(
    "Kraken2 vs mOTUs",
    "QIIME2 vs Kraken2",
    "QIIME2 vs mOTUs"
  )
)]

dt[, metric := factor(
  metric,
  levels = c(
    "BioProject prevalence",
    "Mean relative abundance"
  )
)]

dt[, comparison_type := factor(
  comparison_type,
  levels = c(
    "Direct WGS cross-pipeline comparison",
    "Aggregate cross-collection comparison"
  )
)]

# ============================================================
# LABELS
# ============================================================

dt[, comparison_label := paste0(
  as.character(comparison),
  "\\n",
  n_common_genera_raw,
  " shared genera; ",
  n_common_genera_filtered,
  " after filtering"
)]


dt[, value_label := sprintf("\u03c1 = %.2f", rho)]

dt[, n_label := paste0("n = ", n_common_genera_filtered)]

# Preserve ordering after construction of labels
label_levels <- dt[
  !duplicated(comparison),
  comparison_label
]

dt[, comparison_label := factor(
  comparison_label,
  levels = label_levels
)]

# ============================================================
# EXPORT EXACT DATA USED IN FIGURE
# ============================================================

fwrite(
  dt,
  file.path(out_dir, "Figure5A_plot_data.csv")
)

# ============================================================
# PLOT
# ============================================================

p <- ggplot(
  dt,
  aes(
    x = metric,
    y = rho,
    fill = metric
  )
) +
  geom_col(
    width = 0.62,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = value_label),
    vjust = -0.45,
    size = 3.7
  ) +
  geom_text(
    aes(label = n_label),
    vjust = 1.5,
    size = 3.0
  ) +
  facet_grid(
    comparison_type ~ comparison_label,
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) +
  scale_y_continuous(
    limits = c(0, 0.78),
    breaks = seq(0, 0.8, 0.2),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    title = "A",
    x = NULL,
    y = expression("Spearman " * rho)
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0
    ),
    strip.background = element_rect(
      fill = "grey95",
      colour = "grey50",
      linewidth = 0.4
    ),
    strip.text.x = element_text(
      face = "bold",
      size = 9.5,
      margin = margin(5, 5, 5, 5)
    ),
    strip.text.y.left = element_text(
      face = "bold",
      size = 9.5,
      angle = 90
    ),
    axis.text.x = element_text(
      size = 9,
      lineheight = 0.95
    ),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10),
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.8, "lines"),
    plot.margin = margin(8, 10, 8, 10)
  )

# ============================================================
# SAVE
# ============================================================

ggsave(
  file.path(out_dir, "Figure5A_REVISED_comment4.png"),
  p,
  width = 10.5,
  height = 6.0,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure5A_REVISED_comment4.pdf"),
  p,
  width = 10.5,
  height = 6.0,
  device = cairo_pdf
)

# ============================================================
# AUDIT SUMMARY
# ============================================================

summary_lines <- c(
  "Reviewer 3 Comment 4 - revised Figure 5A",
  "",
  "DIRECT WGS CROSS-PIPELINE COMPARISON:",
  sprintf(
    "Kraken2 vs mOTUs: prevalence rho = %.3f; abundance rho = %.3f; raw genera = %d; filtered genera = %d",
    dt[comparison == "Kraken2 vs mOTUs" & metric == "BioProject prevalence", rho],
    dt[comparison == "Kraken2 vs mOTUs" & metric == "Mean relative abundance", rho],
    dt[comparison == "Kraken2 vs mOTUs", unique(n_common_genera_raw)],
    dt[comparison == "Kraken2 vs mOTUs", unique(n_common_genera_filtered)]
  ),
  "",
  "AGGREGATE CROSS-COLLECTION COMPARISONS:",
  sprintf(
    "QIIME2 vs Kraken2: prevalence rho = %.3f; abundance rho = %.3f; raw genera = %d; filtered genera = %d",
    dt[comparison == "QIIME2 vs Kraken2" & metric == "BioProject prevalence", rho],
    dt[comparison == "QIIME2 vs Kraken2" & metric == "Mean relative abundance", rho],
    dt[comparison == "QIIME2 vs Kraken2", unique(n_common_genera_raw)],
    dt[comparison == "QIIME2 vs Kraken2", unique(n_common_genera_filtered)]
  ),
  sprintf(
    "QIIME2 vs mOTUs: prevalence rho = %.3f; abundance rho = %.3f; raw genera = %d; filtered genera = %d",
    dt[comparison == "QIIME2 vs mOTUs" & metric == "BioProject prevalence", rho],
    dt[comparison == "QIIME2 vs mOTUs" & metric == "Mean relative abundance", rho],
    dt[comparison == "QIIME2 vs mOTUs", unique(n_common_genera_raw)],
    dt[comparison == "QIIME2 vs mOTUs", unique(n_common_genera_filtered)]
  ),
  "",
  "IMPORTANT INTERPRETATION:",
  "- Kraken2 vs mOTUs is the direct cross-pipeline WGS comparison.",
  "- Comparisons involving QIIME2 are aggregate comparisons between independent study collections.",
  "- No mean correlation across the three pairwise comparisons should be reported."
)

writeLines(
  summary_lines,
  file.path(out_dir, "Figure5A_REVISED_audit_summary.txt")
)

cat(paste(summary_lines, collapse = "\n"))
cat("\n\nOutputs written to:\n", out_dir, "\n")