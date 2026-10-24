#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

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
rew <- file.path(
  project_dir,
  "processed_data/geographic_analysis"
)

out <- file.path(
  project_dir,
  "results/figures/figure4"
)

dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helpers
# ============================================================

read_tab <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

write_csv <- function(x, name) {
  path <- file.path(out, name)
  write.csv(x, path, row.names = FALSE, na = "")
  message("Written: ", path)
}

# ============================================================
# S5 — PERMANOVA
# ============================================================

primary <- read_tab(file.path(
  rew,
  "02_confirmed_origin",
  "PERMANOVA_Country_only_ALL_confirmed_origin.tsv"
))

min2 <- read_tab(file.path(
  rew,
  "03_country_min2",
  "PERMANOVA_min2_Country_only_ALL_eligible.tsv"
))

comp <- read_tab(file.path(
  rew,
  "02_confirmed_origin",
  "PERMANOVA_confirmed_origin_country_composition_type.tsv"
))

removed <- read_tab(file.path(
  rew,
  "06_starter_removed",
  "PERMANOVA_country_after_removing_5_starter_genera.tsv"
))

continuous <- read_tab(file.path(
  rew,
  "07_continuous_starter",
  "Country_effect_after_continuous_starter_proportions.tsv"
))

# Standardize arbitrary PERMANOVA column names
std_perm <- function(x) {
  nms <- names(x)

  getcol <- function(candidates) {
    hit <- candidates[candidates %in% nms]
    if (length(hit) == 0) return(rep(NA, nrow(x)))
    x[[hit[1]]]
  }

  data.frame(
    Term = as.character(getcol(c("Term", "term"))),
    Df = as.numeric(getcol(c("Df", "df"))),
    SumOfSqs = as.numeric(getcol(c("SumOfSqs", "SumOfSqs.", "SumOfSquares"))),
    R2 = as.numeric(getcol(c("R2", "R²"))),
    pseudo_F = as.numeric(getcol(c("F", "pseudo_F"))),
    p_value = as.numeric(getcol(c("p", "p_value", "Pr(>F)"))),
    stringsAsFactors = FALSE
  )
}

p1 <- std_perm(primary)
p2 <- std_perm(min2)
pr <- std_perm(removed)
pc <- std_perm(continuous)

# Composition table already contains explicit Model
comp_out <- data.frame(
  Analysis = "Composition-type analysis",
  Subset = "Starter-positive, single-origin BioProjects",
  N_BioProjects = 68,
  N_Countries = 20,
  Model = comp$Model,
  Term = comp$Term,
  Df = comp$Df,
  SumOfSqs = comp$SumOfSqs,
  R2 = comp$R2,
  pseudo_F = comp$F,
  p_value = comp$p,
  Permutations = 999,
  Transformation = "Hellinger",
  Distance = "Bray-Curtis",
  Notes = ifelse(
    grepl("\\+", comp$Model),
    "Sequential PERMANOVA (adonis2, by='terms'); term R2 depends on model order.",
    "Single-factor PERMANOVA."
  ),
  stringsAsFactors = FALSE
)

make_row <- function(
  z,
  analysis,
  subset,
  n,
  countries,
  model,
  notes = ""
) {
  # keep Country/model term rather than residual/total if present
  idx <- which(z$Term %in% c("Country", "country"))
  if (length(idx) == 0) idx <- which(!is.na(z$R2))[1]
  idx <- idx[1]

  data.frame(
    Analysis = analysis,
    Subset = subset,
    N_BioProjects = n,
    N_Countries = countries,
    Model = model,
    Term = z$Term[idx],
    Df = z$Df[idx],
    SumOfSqs = z$SumOfSqs[idx],
    R2 = z$R2[idx],
    pseudo_F = z$pseudo_F[idx],
    p_value = z$p_value[idx],
    Permutations = 999,
    Transformation = "Hellinger",
    Distance = "Bray-Curtis",
    Notes = notes,
    stringsAsFactors = FALSE
  )
}

S5 <- rbind(
  make_row(
    p1,
    "Primary geographic association",
    "All single-origin BioProjects",
    73, 20,
    "Country",
    "Country defined from BioSample-reported sample-origin metadata."
  ),

  make_row(
    p2,
    "Replicated-country sensitivity",
    "Countries represented by >=2 BioProjects",
    63, 10,
    "Country",
    "Sensitivity analysis excluding countries represented by a single BioProject."
  ),

  comp_out,

  make_row(
    pr,
    "Five composition-associated genera removed",
    "All single-origin BioProjects",
    73, 20,
    "Country",
    paste0(
      "Community matrix excludes Lactococcus, Streptococcus, ",
      "Lacticaseibacillus, Lactobacillus and Leuconostoc."
    )
  ),

  make_row(
    pc,
    "Continuous composition sensitivity",
    "Starter-positive, single-origin BioProjects",
    68, 20,
    paste0(
      "Lactococcus + Streptococcus + Lactobacillus + ",
      "Leuconostoc + Country"
    ),
    paste0(
      "Country entered after continuous proportions of four genera; ",
      "Lacticaseibacillus omitted as the reference because the five ",
      "composition proportions sum to 1. Sequential PERMANOVA."
    )
  )
)

write_csv(
  S5,
  "Supplementary_Table_S5_Geographic_PERMANOVA_FINAL.csv"
)

# ============================================================
# S5b — PERMDISP
# ============================================================

disp_primary <- read_tab(file.path(
  rew,
  "02_confirmed_origin",
  "PERMDISP_Country_ALL_confirmed_origin.tsv"
))

disp_min2 <- read_tab(file.path(
  rew,
  "03_country_min2",
  "PERMDISP_min2_Country_ALL_eligible.tsv"
))

disp_comp <- read_tab(file.path(
  rew,
  "02_confirmed_origin",
  "PERMDISP_confirmed_origin.tsv"
))

std_disp <- function(x) {

  nms <- names(x)

  getcol <- function(candidates) {
    hit <- candidates[candidates %in% nms]
    if (length(hit) == 0) return(rep(NA, nrow(x)))
    x[[hit[1]]]
  }

  data.frame(
    Term = as.character(getcol(c("Term", "term", "Grouping"))),
    PERMDISP_F = as.numeric(getcol(c("F", "PERMDISP_F"))),
    PERMDISP_p = as.numeric(getcol(c("p", "p_value", "PERMDISP_p", "Pr(>F)"))),
    stringsAsFactors = FALSE
  )
}

dp <- std_disp(disp_primary)
dm <- std_disp(disp_min2)
dc <- std_disp(disp_comp)

S5b <- data.frame(
  Analysis = c(
    "Primary geographic association",
    "Replicated-country sensitivity",
    "Composition-type subset",
    "Composition-type subset"
  ),
  Subset = c(
    "All single-origin BioProjects",
    "Countries represented by >=2 BioProjects",
    "Starter-positive, single-origin BioProjects",
    "Starter-positive, single-origin BioProjects"
  ),
  N_BioProjects = c(73, 63, 68, 68),
  N_Countries = c(20, 10, 20, 20),
  Grouping_variable = c(
    "Country",
    "Country",
    "Country",
    "Composition type"
  ),
  PERMDISP_F = c(
    dp$PERMDISP_F[1],
    dm$PERMDISP_F[1],
    dc$PERMDISP_F[1],
    dc$PERMDISP_F[2]
  ),
  PERMDISP_p = c(
    dp$PERMDISP_p[1],
    dm$PERMDISP_p[1],
    dc$PERMDISP_p[1],
    dc$PERMDISP_p[2]
  ),
  Permutations = 999,
  Distance = "Bray-Curtis on Hellinger-transformed profiles",
  Interpretation = c(
    "Significant heterogeneity of multivariate dispersion; primary PERMANOVA must not be interpreted solely as centroid separation.",
    "PERMDISP did not reach statistical significance in the replicated-country sensitivity subset.",
    "Significant heterogeneity of dispersion among countries in the composition-type subset.",
    "Dispersion differs among descriptive composition types."
  ),
  stringsAsFactors = FALSE
)

write_csv(
  S5b,
  "Supplementary_Table_S5b_Geographic_PERMDISP_FINAL.csv"
)

# ============================================================
# S7 — balanced geographic sensitivity
# ============================================================

S7 <- read_tab(file.path(
  rew,
  "04_balanced_subsampling",
  "balanced_subsampling_summary.tsv"
))

S7$Interpretation <- ""

idx <- grep("2.5 percentile", S7$Metric, ignore.case = TRUE)
if (length(idx))
  S7$Interpretation[idx] <-
    "Lower bound of empirical percentile interval across subsampling iterations; not a confidence interval."

idx <- grep("97.5 percentile", S7$Metric, ignore.case = TRUE)
if (length(idx))
  S7$Interpretation[idx] <-
    "Upper bound of empirical percentile interval across subsampling iterations; not a confidence interval."

idx <- grep("Fraction PERMANOVA p<=0.05", S7$Metric, fixed = TRUE)
if (length(idx))
  S7$Interpretation[idx] <-
    "Fraction of balanced subsampling iterations with PERMANOVA p <= 0.05."

idx <- grep("Fraction PERMDISP p<=0.05", S7$Metric, fixed = TRUE)
if (length(idx))
  S7$Interpretation[idx] <-
    "Fraction of balanced subsampling iterations with PERMDISP p <= 0.05."

write_csv(
  S7,
  "Supplementary_Table_S7_Balanced_subsampling_summary_FINAL.csv"
)

# ============================================================
# Variance partition — useful source/supplementary table
# ============================================================

vp <- read_tab(file.path(
  rew,
  "02_confirmed_origin",
  "Variance_partition_confirmed_origin.tsv"
))

vp$Percent <- 100 * vp$R2
vp$Definition <- c(
  "Country-associated fraction remaining after accounting for composition type.",
  "Fraction jointly attributable to country and composition type in the sequential partition.",
  "Composition-type-associated fraction remaining after accounting for country.",
  "Unexplained fraction of the community dissimilarity model."
)

write_csv(
  vp,
  "Supplementary_Table_S5c_Variance_partition_FINAL.csv"
)

# ============================================================
# Final index
# ============================================================

index <- data.frame(
  File = c(
    "Supplementary_Table_S4_Geographic_origin.csv",
    "Supplementary_Table_S4b_Country_coverage.csv",
    "Supplementary_Table_S5_Geographic_PERMANOVA_FINAL.csv",
    "Supplementary_Table_S5b_Geographic_PERMDISP_FINAL.csv",
    "Supplementary_Table_S5c_Variance_partition_FINAL.csv",
    "Supplementary_Table_S7_Balanced_subsampling_summary_FINAL.csv",
    "Source_Data_S7_Balanced_subsampling_all_iterations.csv",
    "Source_Data_Figure4_cluster_membership.csv"
  ),
  Description = c(
    "BioProject-level geographic provenance based on BioSample-reported sample origin.",
    "Country-level coverage of BioProjects used in geographic analyses.",
    "PERMANOVA models underlying revised geographic and composition analyses.",
    "PERMDISP diagnostics corresponding to geographic and composition analyses.",
    "Partition of variation between country, composition type and their shared fraction.",
    "Summary of 1000 balanced geographic subsampling iterations.",
    "Complete iteration-level results from balanced geographic subsampling.",
    "BioProject membership of the descriptive composition types used in Figure 4."
  ),
  stringsAsFactors = FALSE
)

write_csv(
  index,
  "Figure4_Supplementary_FINAL_INDEX.csv"
)

message("")
message("======================================================")
message("FINAL FIGURE 4 SUPPLEMENTARY TABLES CREATED")
message("Output: ", out)
message("======================================================")
