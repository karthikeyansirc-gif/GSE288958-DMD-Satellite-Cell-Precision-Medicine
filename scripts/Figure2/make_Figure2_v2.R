
# ============================================================
# FIGURE 2 v2
# DMD-associated satellite-cell State 4
# Publication-style redesign
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(grid)
  library(gridExtra)
})

# ------------------------------------------------------------
# DIRECTORIES
# ------------------------------------------------------------

base_dir <- "/home/ubuntu/scRNAseq/data/GSE288958"
out_dir  <- file.path(base_dir, "Publication", "Figure2_v2")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# FILES
# ------------------------------------------------------------

seurat_file <- file.path(
  base_dir,
  "Satellite_cluster7_reclustered.rds"
)

program_file <- file.path(
  base_dir,
  "STEP10_CORRECTED_State4_14gene_sample_expression.csv"
)

functional_file <- file.path(
  base_dir,
  "STEP14_State4_functional_category_summary.csv"
)

abundance_file <- file.path(
  base_dir,
  "STEP18B_State4_abundance_sample_table.csv"
)

signature_file <- file.path(
  base_dir,
  "STEP19B_State4_DMD_signature_sample_validation.csv"
)

network_file <- file.path(
  base_dir,
  "State4_DoRothEA_AB_regulatory_network.csv"
)

# ------------------------------------------------------------
# PUBLICATION THEME
# ------------------------------------------------------------

pub_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 10,
      color = "grey30"
    ),
    axis.title = element_text(
      size = 11,
      face = "bold"
    ),
    axis.text = element_text(
      size = 9,
      color = "black"
    ),
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(
      size = 9
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    plot.margin = margin(8, 8, 8, 8)
  )

# ------------------------------------------------------------
# COLORS
# ------------------------------------------------------------

group_cols <- c(
  "BMD" = "#4C78A8",
  "DMD" = "#D62728",
  "Normal" = "#777777"
)

category_cols <- c(
  "Cell_cycle_proliferation" = "#D62728",
  "Cytoskeleton_Rho_signaling" = "#1F77B4",
  "Extracellular_matrix_development" = "#2CA02C",
  "Immune_transcriptional_regulation" = "#9467BD",
  "Neural_guidance_signaling" = "#FF7F0E",
  "Unknown_or_context_dependent" = "#7F7F7F"
)

gene_cols <- c(
  "DIAPH3",
  "BTNL8",
  "RRM2",
  "COL27A1",
  "CENPF",
  "UHRF1",
  "ARHGAP22",
  "RFX8",
  "NUSAP1",
  "SEMA3A",
  "CCND1",
  "CCBE1",
  "MKI67",
  "TOP2A"
)

# ============================================================
# PANEL A
# STATE 4 LOCALIZATION UMAP
# ============================================================

message("Creating Panel A...")

obj <- readRDS(seurat_file)

# Find cluster column
cluster_col <- NULL

if ("seurat_clusters" %in% colnames(obj@meta.data)) {
  cluster_col <- "seurat_clusters"
} else {
  possible <- grep(
    "cluster|subcluster",
    colnames(obj@meta.data),
    ignore.case = TRUE,
    value = TRUE
  )
  if (length(possible) > 0) {
    cluster_col <- possible[1]
  }
}

if (is.null(cluster_col)) {
  stop("Could not identify cluster column in Seurat metadata.")
}

clusters <- as.character(obj@meta.data[[cluster_col]])

# State 4 is cluster 4 based on the validated State4 object
state4_flag <- clusters == "4"

# Prefer umap_satellite
reduction_use <- if ("umap_satellite" %in% names(obj@reductions)) {
  "umap_satellite"
} else {
  "umap"
}

emb <- Embeddings(obj, reduction = reduction_use)

umap_df <- data.frame(
  UMAP1 = emb[, 1],
  UMAP2 = emb[, 2],
  State = ifelse(state4_flag, "State 4", "Other satellite cells")
)

pA <- ggplot() +

  geom_point(
    data = umap_df %>% filter(State == "Other satellite cells"),
    aes(UMAP1, UMAP2),
    color = "#D0D0D0",
    size = 0.45,
    alpha = 0.65
  ) +

  geom_point(
    data = umap_df %>% filter(State == "State 4"),
    aes(UMAP1, UMAP2),
    color = "#C62828",
    size = 1.05,
    alpha = 0.95
  ) +

  labs(
    title = "A. State 4 localization",
    subtitle = "Satellite-cell UMAP",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +

  pub_theme +

  theme(
    legend.position = "none",
    aspect.ratio = 1
  )

ggsave(
  file.path(out_dir, "Figure2A_State4_Localization_v2.png"),
  pA,
  width = 5.2,
  height = 4.7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2A_State4_Localization_v2.pdf"),
  pA,
  width = 5.2,
  height = 4.7,
  device = cairo_pdf
)

# ============================================================
# PANEL B
# 14-GENE STATE 4 PROGRAM HEATMAP
# ============================================================

message("Creating Panel B...")

program <- read_csv(
  program_file,
  show_col_types = FALSE
)

genes_present <- intersect(
  gene_cols,
  colnames(program)
)

if (length(genes_present) < 5) {
  stop("Too few State4 genes found in program expression file.")
}

heat_df <- program %>%
  select(
    sample,
    group,
    all_of(genes_present)
  ) %>%
  distinct(sample, .keep_all = TRUE)

# Matrix
heat_mat <- as.matrix(
  heat_df[, genes_present]
)

rownames(heat_mat) <- heat_df$sample

# log transform if needed
heat_mat <- log1p(heat_mat)

# z-score genes
heat_scaled <- t(
  scale(t(heat_mat))
)

heat_long <- as.data.frame(heat_scaled) %>%
  mutate(sample = rownames(.)) %>%
  pivot_longer(
    cols = -sample,
    names_to = "gene",
    values_to = "z"
  ) %>%
  left_join(
    heat_df %>% select(sample, group),
    by = "sample"
  )

# Order samples
group_order <- c("Normal", "BMD", "DMD")

heat_long$group <- factor(
  heat_long$group,
  levels = group_order
)

sample_order <- heat_long %>%
  distinct(sample, group) %>%
  arrange(group, sample) %>%
  pull(sample)

heat_long$sample <- factor(
  heat_long$sample,
  levels = sample_order
)

gene_order <- rev(genes_present)

heat_long$gene <- factor(
  heat_long$gene,
  levels = gene_order
)

pB <- ggplot(
  heat_long,
  aes(
    x = sample,
    y = gene,
    fill = z
  )
) +

  geom_tile(
    color = "white",
    linewidth = 0.35
  ) +

  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Expression\nZ-score"
  ) +

  labs(
    title = "B. State 4 molecular program",
    subtitle = "14-gene expression across biological samples",
    x = NULL,
    y = NULL
  ) +

  pub_theme +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    ),
    axis.text.y = element_text(
      size = 9,
      face = "italic"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
  )

ggsave(
  file.path(out_dir, "Figure2B_State4_14gene_Heatmap_v2.png"),
  pB,
  width = 6.5,
  height = 5.0,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2B_State4_14gene_Heatmap_v2.pdf"),
  pB,
  width = 6.5,
  height = 5.0,
  device = cairo_pdf
)

# ============================================================
# PANEL C
# FUNCTIONAL ARCHITECTURE - LOLLIPOP/DOT PLOT
# ============================================================

message("Creating Panel C...")

functional <- read_csv(
  functional_file,
  show_col_types = FALSE
)

functional <- functional %>%
  mutate(
    category = factor(
      category,
      levels = rev(category)
    )
  )

pC <- ggplot(
  functional,
  aes(
    x = mean_log2FC,
    y = category
  )
) +

  geom_segment(
    aes(
      x = 0,
      xend = mean_log2FC,
      y = category,
      yend = category,
      color = category
    ),
    linewidth = 1
  ) +

  geom_point(
    aes(
      size = n_present,
      fill = category
    ),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +

  scale_color_manual(
    values = category_cols,
    guide = "none"
  ) +

  scale_fill_manual(
    values = category_cols,
    guide = "none"
  ) +

  scale_size_continuous(
    range = c(3, 8),
    name = "Genes present"
  ) +

  labs(
    title = "C. Functional architecture",
    subtitle = "State 4-associated biological programs",
    x = "Mean log2 fold-change",
    y = NULL
  ) +

  pub_theme +

  theme(
    axis.text.y = element_text(
      size = 9
    )
  )

ggsave(
  file.path(out_dir, "Figure2C_State4_Functional_Architecture_v2.png"),
  pC,
  width = 6.4,
  height = 4.7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2C_State4_Functional_Architecture_v2.pdf"),
  pC,
  width = 6.4,
  height = 4.7,
  device = cairo_pdf
)

# ============================================================
# PANEL D
# STATE 4 ABUNDANCE
# ============================================================

message("Creating Panel D...")

abundance <- read_csv(
  abundance_file,
  show_col_types = FALSE
)

abundance <- abundance %>%
  mutate(
    Group = factor(
      Group,
      levels = c("BMD", "DMD", "Normal")
    )
  )

pD <- ggplot(
  abundance,
  aes(
    x = Group,
    y = State4_percent,
    fill = Group
  )
) +

  geom_violin(
    alpha = 0.18,
    color = "black",
    linewidth = 0.6,
    trim = FALSE
  ) +

  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    color = "black",
    fill = "white",
    linewidth = 0.7
  ) +

  geom_jitter(
    width = 0.07,
    size = 3,
    shape = 21,
    color = "black"
  ) +

  scale_fill_manual(
    values = group_cols
  ) +

  labs(
    title = "D. State 4 abundance",
    subtitle = "Biological sample-level validation",
    x = NULL,
    y = "State 4 cells (%)"
  ) +

  pub_theme +

  theme(
    legend.position = "none"
  )

ggsave(
  file.path(out_dir, "Figure2D_State4_Abundance_v2.png"),
  pD,
  width = 4.8,
  height = 4.7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2D_State4_Abundance_v2.pdf"),
  pD,
  width = 4.8,
  height = 4.7,
  device = cairo_pdf
)

# ============================================================
# PANEL E
# STATE 4 PROGRAM SCORE
# ============================================================

message("Creating Panel E...")

signature <- read_csv(
  signature_file,
  show_col_types = FALSE
)

# Use DMD vs Normal only
score_df <- signature %>%
  filter(group %in% c("DMD", "Normal")) %>%
  mutate(
    group = factor(
      group,
      levels = c("DMD", "Normal")
    )
  )

pE <- ggplot(
  score_df,
  aes(
    x = group,
    y = DMD_signature_score,
    fill = group
  )
) +

  geom_violin(
    alpha = 0.20,
    color = "black",
    linewidth = 0.6,
    trim = FALSE
  ) +

  geom_boxplot(
    width = 0.20,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.7
  ) +

  geom_jitter(
    width = 0.07,
    size = 3,
    shape = 21,
    color = "black"
  ) +

  scale_fill_manual(
    values = c(
      "DMD" = "#D62728",
      "Normal" = "#777777"
    )
  ) +

  labs(
    title = "E. State 4 program score",
    subtitle = "DMD versus Normal",
    x = NULL,
    y = "State 4 program score"
  ) +

  pub_theme +

  theme(
    legend.position = "none"
  )

ggsave(
  file.path(out_dir, "Figure2E_State4_Program_Score_v2.png"),
  pE,
  width = 4.8,
  height = 4.7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2E_State4_Program_Score_v2.pdf"),
  pE,
  width = 4.8,
  height = 4.7,
  device = cairo_pdf
)

# ============================================================
# PANEL F
# REGULATORY NETWORK
# ============================================================

message("Creating Panel F...")

network <- read_csv(
  network_file,
  show_col_types = FALSE
)

# Create node positions
tf_nodes <- network %>%
  distinct(tf) %>%
  arrange(tf) %>%
  mutate(
    x = 0,
    y = seq(
      from = n() + 1,
      to = 2,
      length.out = n()
    )
  )

target_nodes <- network %>%
  distinct(target) %>%
  arrange(target) %>%
  mutate(
    x = 1,
    y = seq(
      from = n() + 1,
      to = 2,
      length.out = n()
    )
  )

edges <- network %>%
  left_join(
    tf_nodes,
    by = "tf"
  ) %>%
  rename(
    x_start = x,
    y_start = y
  ) %>%
  left_join(
    target_nodes,
    by = "target"
  ) %>%
  rename(
    x_end = x,
    y_end = y
  )

pF <- ggplot() +

  geom_curve(
    data = edges,
    aes(
      x = x_start,
      y = y_start,
      xend = x_end,
      yend = y_end,
      linetype = factor(mor)
    ),
    curvature = 0.12,
    color = "#555555",
    linewidth = 0.8,
    arrow = arrow(
      length = unit(0.12, "cm"),
      type = "closed"
    )
  ) +

  geom_point(
    data = tf_nodes,
    aes(x, y),
    size = 7,
    shape = 21,
    fill = "#D62728",
    color = "black",
    stroke = 0.8
  ) +

  geom_point(
    data = target_nodes,
    aes(x, y),
    size = 6,
    shape = 21,
    fill = "#4C78A8",
    color = "black",
    stroke = 0.8
  ) +

  geom_text(
    data = tf_nodes,
    aes(
      x = x - 0.06,
      y = y,
      label = tf
    ),
    hjust = 1,
    size = 4,
    fontface = "bold"
  ) +

  geom_text(
    data = target_nodes,
    aes(
      x = x + 0.06,
      y = y,
      label = target
    ),
    hjust = 0,
    size = 3.7
  ) +

  scale_linetype_manual(
    values = c(
      "1" = "solid",
      "-1" = "dashed"
    ),
    labels = c(
      "1" = "Activation",
      "-1" = "Repression"
    ),
    name = "Regulatory mode"
  ) +

  xlim(-0.55, 1.55) +

  labs(
    title = "F. State 4 regulatory architecture",
    subtitle = "DoRothEA-supported TF-target relationships",
    x = NULL,
    y = NULL
  ) +

  pub_theme +

  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  file.path(out_dir, "Figure2F_State4_TF_Network_v2.png"),
  pF,
  width = 6.4,
  height = 4.7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2F_State4_TF_Network_v2.pdf"),
  pF,
  width = 6.4,
  height = 4.7,
  device = cairo_pdf
)

# ============================================================
# COMBINE FIGURE
# ============================================================

message("Combining Figure 2...")

# Resize individual panels for balanced layout
pA2 <- pA + theme(
  plot.title = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 8)
)

pB2 <- pB + theme(
  plot.title = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 8)
)

pC2 <- pC + theme(
  plot.title = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 8)
)

pD2 <- pD + theme(
  plot.title = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 8)
)

pE2 <- pE + theme(
  plot.title = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 8)
)

pF2 <- pF + theme(
  plot.title = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 8)
)

figure2 <- (
  (pA2 | pB2) /
  (pC2 | pD2) /
  (pE2 | pF2)
) +

  plot_annotation(
    title = "Figure 2. DMD-associated satellite-cell State 4",
    theme = theme(
      plot.title = element_text(
        size = 19,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 12
        )
      )
    )
  )

# ------------------------------------------------------------
# SAVE MASTER FIGURES
# ------------------------------------------------------------

ggsave(
  file.path(
    out_dir,
    "Figure2_DMD_Satellite_State4_PUBLICATION_v2.pdf"
  ),
  figure2,
  width = 13,
  height = 16,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Figure2_DMD_Satellite_State4_PUBLICATION_v2.png"
  ),
  figure2,
  width = 13,
  height = 16,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Figure2_DMD_Satellite_State4_PUBLICATION_v2.tiff"
  ),
  figure2,
  width = 13,
  height = 16,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# ============================================================
# SAVE R OBJECT
# ============================================================

saveRDS(
  figure2,
  file.path(
    out_dir,
    "Figure2_DMD_Satellite_State4_PUBLICATION_v2.rds"
  )
)

# ============================================================
# SUMMARY
# ============================================================

cat("\n============================================\n")
cat("FIGURE 2 v2 COMPLETE\n")
cat("============================================\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("Files:\n")
print(list.files(
  out_dir,
  full.names = FALSE
))

cat("\n============================================\n")
cat("DONE\n")
cat("============================================\n")

EOF\
cd /home/ubuntu/scRNAseq/data/GSE288958

cat > make_Figure2_v2.R <<'EOF'

# ============================================================
# FIGURE 2 v2
# DMD-associated satellite-cell State 4
# Publication-style redesign
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(grid)
  library(gridExtra)
})

# ------------------------------------------------------------
# DIRECTORIES
# ------------------------------------------------------------

base_dir <- "/home/ubuntu/scRNAseq/data/GSE288958"
out_dir  <- file.path(base_dir, "Publication", "Figure2_v2")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# FILES
# ------------------------------------------------------------

seurat_file <- file.path(
  base_dir,
  "Satellite_cluster7_reclustered.rds"
)

program_file <- file.path(
  base_dir,
  "STEP10_CORRECTED_State4_14gene_sample_expression.csv"
)

functional_file <- file.path(
  base_dir,
  "STEP14_State4_functional_category_summary.csv"
)

abundance_file <- file.path(
  base_dir,
  "STEP18B_State4_abundance_sample_table.csv"
)

signature_file <- file.path(
  base_dir,
  "STEP19B_State4_DMD_signature_sample_validation.csv"
)

network_file <- file.path(
  base_dir,
  "State4_DoRothEA_AB_regulatory_network.csv"
)

# ------------------------------------------------------------
# PUBLICATION THEME
# ------------------------------------------------------------

pub_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 10,
      color = "grey30"
    ),
    axis.title = element_text(
      size = 11,
      face = "bold"
    ),
    axis.text = element_text(
      size = 9,
      color = "black"
    ),
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(
      size = 9
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    ),
    plot.margin = margin(8, 8, 8, 8)
  )

# ------------------------------------------------------------
# COLORS
# ------------------------------------------------------------

group_cols <- c(
  "BMD" = "#4C78A8",
  "DMD" = "#D62728",
  "Normal" = "#777777"
)

category_cols <- c(
  "Cell_cycle_proliferation" = "#D62728",
  "Cytoskeleton_Rho_signaling" = "#1F77B4",
  "Extracellular_matrix_development" = "#2CA02C",
  "Immune_transcriptional_regulation" = "#9467BD",
  "Neural_guidance_signaling" = "#FF7F0E",
  "Unknown_or_context_dependent" = "#7F7F7F"
)

gene_cols <- c(
  "DIAPH3",
  "BTNL8",
  "RRM2",
  "COL27A1",
  "CENPF",
  "UHRF1",
  "ARHGAP22",
  "RFX8",
  "NUSAP1",
  "SEMA3A",
  "CCND1",
  "CCBE1",
  "MKI67",
  "TOP2A"
)

# ============================================================
# PANEL A
# STATE 4 LOCALIZATION UMAP
# ============================================================

message("Creating Panel A...")

obj <- readRDS(seurat_file)

# Find cluster column
cluster_col <- NULL

if ("seurat_clusters" %in% colnames(obj@meta.data)) {
  cluster_col <- "seurat_clusters"
} else {
  possible <- grep(
    "cluster|subcluster",
    colnames(obj@meta.data),
    ignore.case = TRUE,
    value = TRUE
  )
  if (length(possible) > 0) {
    cluster_col <- possible[1]
  }
}

if (is.null(cluster_col)) {
  stop("Could not identify cluster column in Seurat metadata.")
}

clusters <- as.character(obj@meta.data[[cluster_col]])

# State 4 is cluster 4 based on the validated State4 object
state4_flag <- clusters == "4"

# Prefer umap_satellite
reduction_use <- if ("umap_satellite" %in% names(obj@reductions)) {
  "umap_satellite"
} else {
  "umap"
}

emb <- Embeddings(obj, reduction = reduction_use)

umap_df <- data.frame(
  UMAP1 = emb[, 1],
  UMAP2 = emb[, 2],
  State = ifelse(state4_flag, "State 4", "Other satellite cells")
)

pA <- ggplot() +

  geom_point(
    data = umap_df %>% filter(State == "Other satellite cells"),
    aes(UMAP1, UMAP2),
    color = "#D0D0D0",
    size = 0.45,
    alpha = 0.65
  ) +

  geom_point(
    data = umap_df %>% filter(State == "State 4"),
    aes(UMAP1, UMAP2),
    color = "#C62828",
    size = 1.05,
    alpha = 0.95
  ) +

  labs(
    title = "A. State 4 localization",
    subtitle = "Satellite-cell UMAP",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +

  pub_theme +

  theme(
    legend.position = "none",
    aspect.ratio = 1
  )

ggsave(
  file.path(out_dir, "Figure2A_State4_Localization_v2.png"),
  pA,
  width = 5.2,
  height = 4.7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2A_State4_Localization_v2.pdf"),
  pA,
  width = 5.2,
  height = 4.7,
  device = cairo_pdf
)

# ============================================================
# PANEL B
# 14-GENE STATE 4 PROGRAM HEATMAP
# ============================================================

message("Creating Panel B...")

program <- read_csv(
  program_file,
  show_col_types = FALSE
)

genes_present <- intersect(
  gene_cols,
  colnames(program)
)

if (length(genes_present) < 5) {
  stop("Too few State4 genes found in program expression file.")
}

heat_df <- program %>%
  select(
    sample,
    group,
    all_of(genes_present)
  ) %>%
  distinct(sample, .keep_all = TRUE)

# Matrix
heat_mat <- as.matrix(
  heat_df[, genes_present]
)

rownames(heat_mat) <- heat_df$sample

# log transform if needed
heat_mat <- log1p(heat_mat)

# z-score genes
heat_scaled <- t(
  scale(t(heat_mat))
)

heat_long <- as.data.frame(heat_scaled) %>%
  mutate(sample = rownames(.)) %>%
  pivot_longer(
    cols = -sample,
    names_to = "gene",
    values_to = "z"
  ) %>%
  left_join(
    heat_df %>% select(sample, group),
    by = "sample"
  )

# Order samples
group_order <- c("Normal", "BMD", "DMD")

heat_long$group <- factor(
  heat_long$group,
  levels = group_order
)

sample_order <- heat_long %>%
  distinct(sample, group) %>%
  arrange(group, sample) %>%
  pull(sample)

heat_long$sample <- factor(
  heat_long$sample,
  levels = sample_order
)

gene_order <- rev(genes_present)

heat_long$gene <- factor(
  heat_long$gene,
  levels = gene_order
)

pB <- ggplot(
  heat_long,
  aes(
    x = sample,
    y = gene,
    fill = z
  )
) +

  geom_tile(
    color = "white",
    linewidth = 0.35
  ) +

  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Expression\nZ-score"
  ) +

  labs(
    title = "B. State 4 molecular program",
    subtitle = "14-gene expression across biological samples",
    x = NULL,
    y = NULL
  ) +

  pub_theme +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 9
    ),
    axis.text.y = element_text(
      size = 9,
      face = "italic"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
  )

ggsave(
  file.path(out_dir, "Figure2B_State4_14gene_Heatmap_v2.png"),
  pB,
  width = 6.5,
  height = 5.0,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2B_State4_14gene_Heatmap_v2.pdf"),
  pB,
  width = 6.5,
  height = 5.0,
  device = cairo_pdf
)

# ============================================================
# PANEL C
# FUNCTIONAL ARCHITECTURE - LOLLIPOP/DOT PLOT
# ============================================================

message("Creating Panel C...")

functional <- read_csv(
  functional_file,
  show_col_types = FALSE
)

functional <- functional %>%
  mutate(
    category = factor(
      category,
      levels = rev(category)
    )
  )

pC <- ggplot(
  functional,
  aes(
    x = mean_log2FC,
    y = category
  )
) +

  geom_segment(
    aes(
      x = 0,
      xend = mean_log2FC,
      y = category,
      yend = category,
      color = category
    ),
    linewidth = 1
  ) +

  geom_point(
    aes(
      size = n_present,
      fill = category
    ),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +

  scale_color_manual(
    values = category_cols,
    guide = "none"
  ) +

  scale_fill_manual(
    values = category_cols,
    guide = "none"
  ) +

  scale_size_continuous(
    range = c(3, 8),
    name = "Genes present"
  ) +

  labs(
    title = "C. Functional architecture",
    subtitle = "State 4-associated biological programs",
    x = "Mean log2 fold-change",
    y = NULL
  ) +

  pub_theme +

  theme(
    axis.text.y = element_text(
      size = 9
    )
  )

ggsave(
  file.path(out_dir, "Figure2C_State4_Functional_Architecture_v2.png"),
  pC,
  width = 6.4,
  height = 4.7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2C_State4_Functional_Architecture_v2.pdf"),
  pC,
  width = 6.4,
  height = 4.7,
  device = cairo_pdf
)

# ============================================================
# PANEL D
# STATE 4 ABUNDANCE
# ============================================================

message("Creating Panel D...")

abundance <- read_csv(
  abundance_file,
  show_col_types = FALSE
)

abundance <- abundance %>%
  mutate(
    Group = factor(
      Group,
      levels = c("BMD", "DMD", "Normal")
    )
  )

pD <- ggplot(
  abundance,
  aes(
    x = Group,
    y = State4_percent,
    fill = Group
  )
) +

  geom_violin(
    alpha = 0.18,
    color = "black",
    linewidth = 0.6,
    trim = FALSE
  ) +

  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    color = "black",
    fill = "white",
    linewidth = 0.7
  ) +

  geom_jitter(
    width = 0.07,
    size = 3,
    shape = 21,
    color = "black"
  ) +

  scale_fill_manual(
    values = group_cols
  ) +

  labs(
    title = "D. State 4 abundance",
    subtitle = "Biological sample-level validation",
    x = NULL,
    y = "State 4 cells (%)"
  ) +

  pub_theme +

  theme(
    legend.position = "none"
  )

ggsave(
  file.path(out_dir, "Figure2D_State4_Abundance_v2.png"),
  pD,
  width = 4.8,
  height = 4.7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2D_State4_Abundance_v2.pdf"),
  pD,
  width = 4.8,
  height = 4.7,
  device = cairo_pdf
)

# ============================================================
# PANEL E
# STATE 4 PROGRAM SCORE
# ============================================================

message("Creating Panel E...")

signature <- read_csv(
  signature_file,
  show_col_types = FALSE
)

# Use DMD vs Normal only
score_df <- signature %>%
  filter(group %in% c("DMD", "Normal")) %>%
  mutate(
    group = factor(
      group,
      levels = c("DMD", "Normal")
    )
  )

pE <- ggplot(
  score_df,
  aes(
    x = group,
    y = DMD_signature_score,
    fill = group
  )
) +

  geom_violin(
    alpha = 0.20,
    color = "black",
    linewidth = 0.6,
    trim = FALSE
  ) +

  geom_boxplot(
    width = 0.20,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.7
  ) +

  geom_jitter(
    width = 0.07,
    size = 3,
    shape = 21,
    color = "black"
  ) +

  scale_fill_manual(
    values = c(
      "DMD" = "#D62728",
      "Normal" = "#777777"
    )
  ) +

  labs(
    title = "E. State 4 program score",
    subtitle = "DMD versus Normal",
    x = NULL,
    y = "State 4 program score"
  ) +

  pub_theme +

  theme(
    legend.position = "none"
  )

ggsave(
  file.path(out_dir, "Figure2E_State4_Program_Score_v2.png"),
  pE,
  width = 4.8,
  height = 4.7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2E_State4_Program_Score_v2.pdf"),
  pE,
  width = 4.8,
  height = 4.7,
  device = cairo_pdf
)

# ============================================================
# PANEL F
# REGULATORY NETWORK
# ============================================================

message("Creating Panel F...")

network <- read_csv(
  network_file,
  show_col_types = FALSE
)

# Create node positions
tf_nodes <- network %>%
  distinct(tf) %>%
  arrange(tf) %>%
  mutate(
    x = 0,
    y = seq(
      from = n() + 1,
      to = 2,
      length.out = n()
    )
  )

target_nodes <- network %>%
  distinct(target) %>%
  arrange(target) %>%
  mutate(
    x = 1,
    y = seq(
      from = n() + 1,
      to = 2,
      length.out = n()
    )
  )

edges <- network %>%
  left_join(
    tf_nodes,
    by = "tf"
  ) %>%
  rename(
    x_start = x,
    y_start = y
  ) %>%
  left_join(
    target_nodes,
    by = "target"
  ) %>%
  rename(
    x_end = x,
    y_end = y
  )

pF <- ggplot() +

  geom_curve(
    data = edges,
    aes(
      x = x_start,
      y = y_start,
      xend = x_end,
      yend = y_end,
      linetype = factor(mor)
    ),
    curvature = 0.12,
    color = "#555555",
    linewidth = 0.8,
    arrow = arrow(
      length = unit(0.12, "cm"),
      type = "closed"
    )
  ) +

  geom_point(
    data = tf_nodes,
    aes(x, y),
    size = 7,
    shape = 21,
    fill = "#D62728",
    color = "black",
    stroke = 0.8
  ) +

  geom_point(
    data = target_nodes,
    aes(x, y),
    size = 6,
    shape = 21,
    fill = "#4C78A8",
    color = "black",
    stroke = 0.8
  ) +

  geom_text(
    data = tf_nodes,
    aes(
      x = x - 0.06,
      y = y,
      label = tf
    ),
    hjust = 1,
    size = 4,
    fontface = "bold"
  ) +

  geom_text(
    data = target_nodes,
    aes(
      x = x + 0.06,
      y = y,
      label = target
    ),
    hjust = 0,
    size = 3.7
  ) +

  scale_linetype_manual(
    values = c(
      "1" = "solid",
      "-1" = "dashed"
    ),
    labels = c(
      "1" = "Activation",
      "-1" = "Repression"
    ),
    name = "Regulatory mode"
  ) +

  xlim(-0.55, 1.55) +

  labs(
    title = "F. State 4 regulatory architecture",
    subtitle = "DoRothEA-supported TF-target relationships",
    x = NULL,
    y = NULL
  ) +

  pub_theme +

  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

ggsave(
  file.path(out_dir, "Figure2F_State4_TF_Network_v2.png"),
  pF,
  width = 6.4,
  height = 4.7,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, "Figure2F_State4_TF_Network_v2.pdf"),
  pF,
  width = 6.4,
  height = 4.7,
  device = cairo_pdf
)

# ============================================================
# COMBINE FIGURE
# ============================================================

message("Combining Figure 2...")

# Resize individual panels for balanced layout
pA2 <- pA + theme(
  plot.title = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 8)
)

pB2 <- pB + theme(
  plot.title = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 8)
)

pC2 <- pC + theme(
  plot.title = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 8)
)

pD2 <- pD + theme(
  plot.title = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 8)
)

pE2 <- pE + theme(
  plot.title = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 8)
)

pF2 <- pF + theme(
  plot.title = element_text(size = 13, face = "bold"),
  plot.subtitle = element_text(size = 8)
)

figure2 <- (
  (pA2 | pB2) /
  (pC2 | pD2) /
  (pE2 | pF2)
) +

  plot_annotation(
    title = "Figure 2. DMD-associated satellite-cell State 4",
    theme = theme(
      plot.title = element_text(
        size = 19,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 12
        )
      )
    )
  )

# ------------------------------------------------------------
# SAVE MASTER FIGURES
# ------------------------------------------------------------

ggsave(
  file.path(
    out_dir,
    "Figure2_DMD_Satellite_State4_PUBLICATION_v2.pdf"
  ),
  figure2,
  width = 13,
  height = 16,
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Figure2_DMD_Satellite_State4_PUBLICATION_v2.png"
  ),
  figure2,
  width = 13,
  height = 16,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "Figure2_DMD_Satellite_State4_PUBLICATION_v2.tiff"
  ),
  figure2,
  width = 13,
  height = 16,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# ============================================================
# SAVE R OBJECT
# ============================================================

saveRDS(
  figure2,
  file.path(
    out_dir,
    "Figure2_DMD_Satellite_State4_PUBLICATION_v2.rds"
  )
)

# ============================================================
# SUMMARY
# ============================================================

cat("\n============================================\n")
cat("FIGURE 2 v2 COMPLETE\n")
cat("============================================\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("Files:\n")
print(list.files(
  out_dir,
  full.names = FALSE
))

cat("\n============================================\n")
cat("DONE\n")
cat("============================================\n")

