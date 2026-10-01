# ============================================================
# FIGURE 3 v2
# DMD-ASSOCIATED MOLECULAR SIGNATURE OF SATELLITE STATE 4
# Publication-quality version
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(scales)
  library(grid)
})

# ------------------------------------------------------------
# DIRECTORIES
# ------------------------------------------------------------

base_dir <- "/home/ubuntu/scRNAseq/data/GSE288958"

out_dir <- file.path(
  base_dir,
  "Publication",
  "Figure3_v2"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# SOURCE FILES
# ------------------------------------------------------------

deg_file <- file.path(
  base_dir,
  "Satellite_cluster7_DEG_DMD_vs_Normal_edgeR_annotated.csv"
)

topdeg_file <- file.path(
  base_dir,
  "Satellite_cluster7_DEG_Top30_DMD_vs_Normal.csv"
)

signature_file <- file.path(
  base_dir,
  "Satellite_cluster7_DMD_19gene_signature_scores.csv"
)

loso_file <- file.path(
  base_dir,
  "Satellite_cluster7_DMD_19gene_LOSO_robustness.csv"
)

heatmap_file <- file.path(
  base_dir,
  "Satellite_cluster7_Significant_DEG_Sample_Expression_scaled.csv"
)

# ------------------------------------------------------------
# CHECK FILES
# ------------------------------------------------------------

source_files <- c(
  deg_file,
  topdeg_file,
  signature_file,
  loso_file,
  heatmap_file
)

cat("\n===== SOURCE FILE CHECK =====\n")

for (f in source_files) {
  cat(
    basename(f),
    " : ",
    file.exists(f),
    "\n",
    sep = ""
  )
}

if (!all(file.exists(source_files))) {
  stop("One or more required source files are missing.")
}

# ------------------------------------------------------------
# READ DATA
# ------------------------------------------------------------

deg <- read.csv(
  deg_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

topdeg <- read.csv(
  topdeg_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

signature <- read.csv(
  signature_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

loso <- read.csv(
  loso_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

heat <- read.csv(
  heatmap_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# CLEAN GENE COLUMN
# ------------------------------------------------------------

if ("" %in% colnames(deg)) {
  colnames(deg)[colnames(deg) == ""] <- "gene"
}

if ("" %in% colnames(topdeg)) {
  colnames(topdeg)[colnames(topdeg) == ""] <- "gene"
}

if ("" %in% colnames(heat)) {
  colnames(heat)[colnames(heat) == ""] <- "gene"
}

# ------------------------------------------------------------
# GENERAL THEME
# ------------------------------------------------------------

pub_theme <- theme_classic(base_size = 15) +

  theme(

    axis.title = element_text(
      size = 17,
      face = "bold",
      colour = "black"
    ),

    axis.text = element_text(
      size = 13,
      face = "bold",
      colour = "black"
    ),

    plot.title = element_text(
      size = 17,
      face = "bold",
      hjust = 0,
      colour = "black"
    ),

    plot.subtitle = element_text(
      size = 12,
      colour = "grey20"
    ),

    legend.title = element_text(
      size = 13,
      face = "bold"
    ),

    legend.text = element_text(
      size = 12,
      face = "bold"
    ),

    axis.line = element_line(
      linewidth = 0.9,
      colour = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.8,
      colour = "black"
    ),

    axis.ticks.length = unit(
      0.18,
      "cm"
    ),

    panel.border = element_blank(),

    plot.margin = margin(
      12,
      12,
      12,
      12
    )
  )

# ============================================================
# PANEL A
# VOLCANO PLOT
# ============================================================

cat("\nCreating Panel A...\n")

deg_plot <- deg %>%
  mutate(
    gene = as.character(gene),
    FDR = as.numeric(FDR),
    logFC = as.numeric(logFC),
    neglog10FDR = -log10(pmax(FDR, 1e-300))
  )

deg_plot <- deg_plot %>%
  mutate(
    status = case_when(

      FDR < 0.05 &
        logFC > 0 ~ "DMD-up",

      FDR < 0.05 &
        logFC < 0 ~ "Normal-up",

      TRUE ~ "Not significant"
    )
  )

label_genes <- c(
  "BTNL8",
  "SORCS1",
  "COL27A1",
  "CCND1",
  "ARHGAP22",
  "SEMA3A",
  "CCBE1",
  "MAPK4",
  "SYT1",
  "TACR3"
)

label_df <- deg_plot %>%
  filter(gene %in% label_genes)

pA <- ggplot(
  deg_plot,
  aes(
    x = logFC,
    y = neglog10FDR
  )
) +

  geom_point(
    data = subset(
      deg_plot,
      status == "Not significant"
    ),
    size = 2.0,
    alpha = 0.55,
    colour = "grey70"
  ) +

  geom_point(
    data = subset(
      deg_plot,
      status == "DMD-up"
    ),
    size = 2.5,
    alpha = 0.75,
    colour = "#C62828"
  ) +

  geom_point(
    data = subset(
      deg_plot,
      status == "Normal-up"
    ),
    size = 2.5,
    alpha = 0.75,
    colour = "#1565C0"
  ) +

  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey35"
  ) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey35"
  ) +

  geom_text(
    data = label_df,
    aes(label = gene),
    size = 4.0,
    fontface = "bold",
    nudge_y = 0.18,
    check_overlap = TRUE
  ) +

  labs(
    title = "DMD-associated differential expression",
    x = "log2 fold change",
    y = expression(-log[10]~FDR)
  ) +

  scale_x_continuous(
    expand = expansion(
      mult = c(0.03, 0.08)
    )
  ) +

  pub_theme

# ============================================================
# PANEL B
# SAMPLE-LEVEL EXPRESSION HEATMAP
# ============================================================

cat("\nCreating Panel B...\n")

heat_long <- heat %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expression"
  )

sample_order <- c(
  "Normal_1",
  "Normal_2",
  "Normal_3",
  "Normal_4",
  "Normal_5",
  "BMD_1",
  "BMD_2",
  "BMD_3",
  "DMD_1",
  "DMD_2",
  "DMD_3"
)

heat_long$sample <- factor(
  heat_long$sample,
  levels = sample_order
)

heat_long$gene <- factor(
  heat_long$gene,
  levels = rev(unique(heat_long$gene))
)

pB <- ggplot(
  heat_long,
  aes(
    x = sample,
    y = gene,
    fill = expression
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.35
  ) +

  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Scaled\nexpression"
  ) +

  labs(
    title = "Sample-level expression of DMD-associated genes",
    x = NULL,
    y = NULL
  ) +

  theme_minimal(base_size = 14) +

  theme(

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11,
      face = "bold",
      colour = "black"
    ),

    axis.text.y = element_text(
      size = 11,
      face = "bold",
      colour = "black"
    ),

    axis.title = element_text(
      size = 15,
      face = "bold"
    ),

    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0
    ),

    legend.title = element_text(
      size = 12,
      face = "bold"
    ),

    legend.text = element_text(
      size = 10
    ),

    panel.grid = element_blank(),

    plot.margin = margin(
      12,
      15,
      12,
      12
    )
  )

# ============================================================
# PANEL C
# 19-GENE SIGNATURE SCORE
# ============================================================

cat("\nCreating Panel C...\n")

signature <- signature %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "Normal",
        "BMD",
        "DMD"
      )
    )
  )

pC <- ggplot(
  signature,
  aes(
    x = group,
    y = DMD_signature_score
  )
) +

  geom_violin(
    aes(fill = group),
    alpha = 0.20,
    colour = "black",
    linewidth = 0.8,
    trim = FALSE
  ) +

  geom_jitter(
    width = 0.08,
    size = 4.2,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 1.0
  ) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey40"
  ) +

  scale_fill_manual(
    values = c(
      "Normal" = "#4D79A7",
      "BMD" = "#59A14F",
      "DMD" = "#E15759"
    )
  ) +

  labs(
    title = "19-gene DMD molecular signature",
    x = NULL,
    y = "Signature score"
  ) +

  pub_theme +

  theme(
    legend.position = "none"
  )

# ============================================================
# PANEL D
# LOSO ROBUSTNESS
# ============================================================

cat("\nCreating Panel D...\n")

loso <- loso %>%
  mutate(
    held_out_sample = as.character(
      held_out_sample
    ),

    group = factor(
      group,
      levels = c(
        "Normal",
        "DMD"
      )
    )
  )

loso_long <- bind_rows(

  loso %>%
    transmute(
      sample = held_out_sample,
      group = group,
      metric = "Held-out score",
      value = held_out_score
    ),

  loso %>%
    transmute(
      sample = held_out_sample,
      group = group,
      metric = "Training-group reference",
      value = ifelse(
        group == "Normal",
        training_Normal_mean,
        training_DMD_mean
      )
    )
)

pD <- ggplot(
  loso_long,
  aes(
    x = sample,
    y = value,
    group = metric
  )
) +

  geom_line(
    aes(
      colour = metric
    ),
    linewidth = 1.0
  ) +

  geom_point(
    aes(
      colour = metric
    ),
    size = 4
  ) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey40"
  ) +

  scale_colour_manual(
    values = c(
      "Held-out score" = "#C62828",
      "Training-group reference" = "#424242"
    )
  ) +

  labs(
    title = "Leave-one-sample-out robustness",
    x = NULL,
    y = "Signature score",
    colour = NULL
  ) +

  pub_theme +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11
    ),

    legend.position = "bottom"
  )

# ============================================================
# PANEL E
# MOLECULAR ARCHITECTURE
# ============================================================

cat("\nCreating Panel E...\n")

architecture <- data.frame(

  category = c(
    "Immune / transcriptional",
    "ECM / development",
    "Cytoskeleton / Rho",
    "Neural guidance",
    "Cell cycle / proliferation"
  ),

  genes = c(
    "BTNL8",
    "COL27A1, CCBE1",
    "ARHGAP22, DIAPH3",
    "SEMA3A",
    "CCND1, RRM2, CENPF, UHRF1, NUSAP1, MKI67, TOP2A"
  ),

  stringsAsFactors = FALSE
)

architecture$n_genes <- sapply(
  strsplit(
    architecture$genes,
    ", "
  ),
  length
)

architecture$category <- factor(
  architecture$category,
  levels = rev(
    architecture$category
  )
)

pE <- ggplot(
  architecture,
  aes(
    x = n_genes,
    y = category
  )
) +

  geom_segment(
    aes(
      x = 0,
      xend = n_genes,
      y = category,
      yend = category
    ),
    linewidth = 1.8,
    colour = "grey65"
  ) +

  geom_point(
    aes(
      size = n_genes
    ),
    shape = 21,
    fill = "#C62828",
    colour = "black",
    stroke = 1.0
  ) +

  geom_text(
    aes(
      label = genes
    ),
    hjust = 0,
    nudge_x = 0.20,
    size = 3.6,
    fontface = "bold"
  ) +

  scale_size_continuous(
    range = c(
      6,
      13
    )
  ) +

  labs(
    title = "Functional architecture of the DMD-associated program",
    x = "Genes represented",
    y = NULL
  ) +

  pub_theme +

  theme(
    legend.position = "none",

    axis.text.y = element_text(
      size = 12
    ),

    plot.margin = margin(
      12,
      120,
      12,
      12
    )
  ) +

  coord_cartesian(
    clip = "off"
  )

# ============================================================
# PANEL F
# INTEGRATED SIGNATURE VALIDATION
# ============================================================

cat("\nCreating Panel F...\n")

signature_summary <- signature %>%
  group_by(group) %>%
  summarise(
    mean_score = mean(
      DMD_signature_score,
      na.rm = TRUE
    ),

    sd_score = sd(
      DMD_signature_score,
      na.rm = TRUE
    ),

    n = n(),

    .groups = "drop"
  )

pF <- ggplot(
  signature,
  aes(
    x = DMD_signature_score,
    y = group
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.8,
    colour = "grey40"
  ) +

  geom_point(
    size = 4.5,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 1.0,
    position = position_jitter(
      height = 0.10
    )
  ) +

  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 6,
    fill = "#C62828",
    colour = "black"
  ) +

  labs(
    title = "Sample-level validation of the DMD signature",
    x = "19-gene signature score",
    y = NULL
  ) +

  pub_theme

# ============================================================
# SAVE INDIVIDUAL PANELS
# ============================================================

cat("\n===== SAVING INDIVIDUAL PANELS =====\n")

panel_list <- list(
  A = pA,
  B = pB,
  C = pC,
  D = pD,
  E = pE,
  F = pF
)

for (nm in names(panel_list)) {

  ggsave(
    file.path(
      out_dir,
      paste0(
        "Figure3",
        nm,
        "_v2.pdf"
      )
    ),
    panel_list[[nm]],
    width = 7.2,
    height = 6.0,
    units = "in",
    device = cairo_pdf
  )

  ggsave(
    file.path(
      out_dir,
      paste0(
        "Figure3",
        nm,
        "_v2.png"
      )
    ),
    panel_list[[nm]],
    width = 7.2,
    height = 6.0,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

# ============================================================
# COMBINE FIGURE
# ============================================================

cat("\n===== COMBINING FIGURE 3 =====\n")

figure3 <- (

  pA | pB

) / (

  pC | pD

) / (

  pE | pF

) +

  plot_annotation(
    tag_levels = "A"
  ) &

  theme(
    plot.tag = element_text(
      size = 20,
      face = "bold"
    )
  )

# ------------------------------------------------------------
# SAVE COMBINED
# ------------------------------------------------------------

pdf_file <- file.path(
  out_dir,
  "Figure3_DMD_Molecular_Signature_PUBLICATION_v2.pdf"
)

png_file <- file.path(
  out_dir,
  "Figure3_DMD_Molecular_Signature_PUBLICATION_v2.png"
)

tiff_file <- file.path(
  out_dir,
  "Figure3_DMD_Molecular_Signature_PUBLICATION_v2.tiff"
)

ggsave(
  pdf_file,
  figure3,
  width = 14,
  height = 17,
  units = "in",
  device = cairo_pdf
)

ggsave(
  png_file,
  figure3,
  width = 14,
  height = 17,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  tiff_file,
  figure3,
  width = 14,
  height = 17,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# ============================================================
# SAVE RDS
# ============================================================

saveRDS(
  figure3,
  file.path(
    out_dir,
    "Figure3_DMD_Molecular_Signature_PUBLICATION_v2.rds"
  )
)

# ============================================================
# SAVE SOURCE SUMMARY
# ============================================================

summary_file <- file.path(
  out_dir,
  "Figure3_source_summary.txt"
)

sink(summary_file)

cat("FIGURE 3 SOURCE SUMMARY\n")
cat("========================\n\n")

cat("DEG source:\n")
cat(deg_file, "\n\n")

cat("19-gene signature source:\n")
cat(signature_file, "\n\n")

cat("LOSO source:\n")
cat(loso_file, "\n\n")

cat("Heatmap source:\n")
cat(heatmap_file, "\n\n")

cat("Number of DEG rows:", nrow(deg), "\n")
cat("Number of signature samples:", nrow(signature), "\n")
cat("Number of LOSO samples:", nrow(loso), "\n")
cat("Number of heatmap genes:", nrow(heat), "\n")

sink()

# ============================================================
# FINAL VERIFICATION
# ============================================================

cat("\n============================================\n")
cat("FIGURE 3 v2 COMPLETE\n")
cat("============================================\n\n")

files_created <- list.files(
  out_dir,
  full.names = TRUE
)

print(
  basename(files_created)
)

cat("\n===== FILE SIZES =====\n")

for (f in files_created) {

  cat(
    basename(f),
    " : ",
    file.info(f)$size,
    " bytes\n",
    sep = ""
  )
}

cat("\n===== FILE EXISTENCE =====\n")

print(
  file.exists(files_created)
)

cat("\n============================================\n")
cat("DONE\n")
cat("============================================\n")

