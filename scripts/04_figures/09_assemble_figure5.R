#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(magick)
})

# ============================================================
# REPOSITORY ROOT
# ============================================================

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)

if (length(file_arg) == 1) {
  script_file <- normalizePath(sub("^--file=", "", file_arg))
  repo_dir <- normalizePath(file.path(dirname(script_file), "..", ".."))
} else {
  repo_dir <- normalizePath(getwd())
}

figure_dir <- file.path(
  repo_dir,
  "results/figures/figure5"
)

# ============================================================
# INPUT PANELS
# ============================================================

panel_A <- file.path(
  figure_dir,
  "Figure5A_REVISED_comment4.png"
)

panel_B <- file.path(
  figure_dir,
  "Figure5B_REVISED_widespread_species.png"
)

panel_C <- file.path(
  figure_dir,
  "Figure5C_species_richness_by_genus.png"
)

panel_D <- file.path(
  figure_dir,
  "Figure5D_species_Kraken2_mOTUs_concordance.png"
)

input_files <- c(
  A = panel_A,
  B = panel_B,
  C = panel_C,
  D = panel_D
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing Figure 5 panels:\n",
    paste(
      names(missing_files),
      missing_files,
      sep = ": ",
      collapse = "\n"
    )
  )
}

message("Figure 5 input panels:")
for (nm in names(input_files)) {
  message("  ", nm, ": ", input_files[[nm]])
}

# ============================================================
# READ PANELS
# ============================================================

A <- image_read(panel_A)
B <- image_read(panel_B)
C <- image_read(panel_C)
D <- image_read(panel_D)

# Flatten in case any PNG contains transparency.
A <- image_flatten(A)
B <- image_flatten(B)
C <- image_flatten(C)
D <- image_flatten(D)

# ============================================================
# NORMALIZE PANEL WIDTHS
#
# A and B form the upper row.
# C and D form the lower row.
#
# Each panel is resized to the same width while preserving its
# aspect ratio. Rows are then padded to equal height.
# ============================================================

target_width <- 3000

resize_panel <- function(img, width = target_width) {
  image_resize(
    img,
    paste0(width, "x")
  )
}

A <- resize_panel(A)
B <- resize_panel(B)
C <- resize_panel(C)
D <- resize_panel(D)

pad_to_height <- function(img, target_height) {

  info <- image_info(img)

  if (info$height == target_height) {
    return(img)
  }

  image_extent(
    img,
    geometry = paste0(
      info$width,
      "x",
      target_height
    ),
    gravity = "center",
    color = "white"
  )
}

# ============================================================
# TOP ROW: A + B
# ============================================================

top_height <- max(
  image_info(A)$height,
  image_info(B)$height
)

A <- pad_to_height(A, top_height)
B <- pad_to_height(B, top_height)

top_row <- image_append(
  c(A, B),
  stack = FALSE
)

# ============================================================
# BOTTOM ROW: C + D
# ============================================================

bottom_height <- max(
  image_info(C)$height,
  image_info(D)$height
)

C <- pad_to_height(C, bottom_height)
D <- pad_to_height(D, bottom_height)

bottom_row <- image_append(
  c(C, D),
  stack = FALSE
)

# ============================================================
# NORMALIZE ROW WIDTHS
# ============================================================

row_width <- max(
  image_info(top_row)$width,
  image_info(bottom_row)$width
)

pad_row <- function(img, target_width) {

  info <- image_info(img)

  if (info$width == target_width) {
    return(img)
  }

  image_extent(
    img,
    geometry = paste0(
      target_width,
      "x",
      info$height
    ),
    gravity = "center",
    color = "white"
  )
}

top_row <- pad_row(
  top_row,
  row_width
)

bottom_row <- pad_row(
  bottom_row,
  row_width
)

# ============================================================
# FINAL FIGURE
# ============================================================

figure5 <- image_append(
  c(top_row, bottom_row),
  stack = TRUE
)

# Small white border around complete figure.
figure5 <- image_border(
  figure5,
  color = "white",
  geometry = "40x40"
)

# ============================================================
# SAVE
# ============================================================

png_out <- file.path(
  figure_dir,
  "Figure5_FINAL_REVISED.png"
)

pdf_out <- file.path(
  figure_dir,
  "Figure5_FINAL_REVISED.pdf"
)

image_write(
  figure5,
  path = png_out,
  format = "png"
)

image_write(
  figure5,
  path = pdf_out,
  format = "pdf"
)

# ============================================================
# AUDIT
# ============================================================

final_info <- image_info(figure5)

audit <- data.frame(
  output = c(
    "Figure5_FINAL_REVISED.png",
    "Figure5_FINAL_REVISED.pdf"
  ),
  width_px = final_info$width,
  height_px = final_info$height,
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  file.path(
    figure_dir,
    "Figure5_FINAL_REVISED_audit.csv"
  ),
  row.names = FALSE
)

message("")
message("Figure 5 assembly completed.")
message(
  "Final dimensions: ",
  final_info$width,
  " x ",
  final_info$height,
  " px"
)
message("PNG: ", png_out)
message("PDF: ", pdf_out)
