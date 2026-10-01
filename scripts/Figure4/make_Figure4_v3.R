############################################################
# FIGURE 4 v4
# DMD SATELLITE-CELL STATE ANALYSIS
#
# A — Global cell-state architecture
# B — Sample-level composition
# C — State 4 abundance
# D — State 4 DMD molecular signature
# E — Integrated biological model
############################################################

rm(list = ls())

options(stringsAsFactors = FALSE)

cat("\n")
cat("============================================\n")
cat("FIGURE 4 v4 — PUBLICATION VERSION\n")
cat("============================================\n\n")

############################################################
# 1. PACKAGES
############################################################

required <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "patchwork",
  "scales"
)

for (p in required) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)
library(scales)

############################################################
# 2. DIRECTORIES
############################################################

base_dir <- "/home/ubuntu/scRNAseq/data/GSE288958"

pub_dir <- file.path(
  base_dir,
  "Publication"
)

out_dir <- file.path(
  pub_dir,
  "Figure4_v4"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 3. INPUT FILES
############################################################

file_17 <- file.path(
  pub_dir,
  "17_CellState_Composition_Percentage.csv"
)

file_16 <- file.path(
  pub_dir,
  "16_CellState_Composition_Long.csv"
)

file_19 <- file.path(
  pub_dir,
  "19_CellState_Replicate_Support.csv"
)

file_state4 <- file.path(
  base_dir,
  "STEP18B_State4_abundance_sample_table.csv"
)

file_signature <- file.path(
  base_dir,
  "Satellite_cluster7_DMD_19gene_signature_scores.csv"
)

cat("===== INPUT FILE CHECK =====\n")

cat(
  "17_CellState_Composition_Percentage.csv :",
  file.exists(file_17), "\n"
)

cat(
  "16_CellState_Composition_Long.csv :",
  file.exists(file_16), "\n"
)

cat(
  "19_CellState_Replicate_Support.csv :",
  file.exists(file_19), "\n"
)

cat(
  "STEP18B_State4_abundance_sample_table.csv :",
  file.exists(file_state4), "\n"
)

cat(
  "Satellite_cluster7_DMD_19gene_signature_scores.csv :",
  file.exists(file_signature), "\n"
)

if (!all(
  file.exists(file_17),
  file.exists(file_16),
  file.exists(file_19),
  file.exists(file_state4),
  file.exists(file_signature)
)) {
  stop("One or more required input files are missing.")
}

############################################################
# 4. READ DATA
############################################################

comp <- read_csv(
  file_17,
  show_col_types = FALSE
)

comp_long <- read_csv(
  file_16,
  show_col_types = FALSE
)

replicate <- read_csv(
  file_19,
  show_col_types = FALSE
)

state4 <- read_csv(
  file_state4,
  show_col_types = FALSE
)

signature <- read_csv(
  file_signature,
  show_col_types = FALSE
)

cat("\n===== DATA STRUCTURE =====\n")

cat("\nComposition:\n")
print(colnames(comp))

cat("\nComposition long:\n")
print(colnames(comp_long))

cat("\nReplicate:\n")
print(colnames(replicate))

cat("\nState4:\n")
print(colnames(state4))

cat("\nSignature:\n")
print(colnames(signature))

############################################################
# 5. COLORS
############################################################

group_cols <- c(
  "Normal" = "#4C78A8",
  "BMD"    = "#F2A541",
  "DMD"    = "#C44E52"
)

celltype_cols <- c(
  "Activated fibroblast" =
    "#6A3D9A",
  "Adipocyte" =
    "#E6AB02",
  "Endothelial" =
    "#1F78B4",
  "Fast myofiber" =
    "#33A02C",
  "Fibroblast/ECM" =
    "#B15928",
  "Macrophage" =
    "#E31A1C",
  "Mast cell" =
    "#FB9A99",
  "Muscle-associated" =
    "#A6CEE3",
  "Myofiber" =
    "#66A61E",
  "Myogenic progenitor" =
    "#1B9E77",
  "Myogenic/ECM-associated" =
    "#7570B3",
  "Myogenic/mesenchymal progenitor" =
    "#D95F02",
  "Neural-associated" =
    "#E7298A",
  "Neutrophil/inflammatory" =
    "#A6761D",
  "Pericyte" =
    "#666666",
  "Progenitor/stromal-associated" =
    "#8DD3C7",
  "Regenerating myofiber" =
    "#80B1D3",
  "Satellite cells" =
    "#FB8072",
  "Satellite/myogenic progenitor" =
    "#FDB462",
  "Slow/oxidative myofiber" =
    "#B3DE69",
  "Specialized myogenic" =
    "#BC80BD",
  "T cell" =
    "#CCEBC5",
  "Vascular-associated" =
    "#FFED6F",
  "Activated macrophage" =
    "#C51B7D",
  "Metabolic myofiber" =
    "#5E3C99"
)

############################################################
# 6. PUBLICATION THEME
############################################################

theme_pub <- theme_classic(
  base_size = 14
) +
  theme(
    plot.title = element_text(
      size = 17,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 11,
      hjust = 0
    ),
    axis.title = element_text(
      size = 14,
      face = "bold"
    ),
    axis.text = element_text(
      size = 11
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
      10, 12, 10, 12
    )
  )

############################################################
# PANEL A
# GLOBAL CELL-STATE ARCHITECTURE
############################################################

cat("\nCreating Panel A...\n")

# IMPORTANT:
# The composition file is WIDE:
# celltype | BMD | DMD | Normal
#
# Convert to:
# celltype | group | Percentage

datA <- comp %>%
  pivot_longer(
    cols = c(
      "Normal",
      "BMD",
      "DMD"
    ),
    names_to = "group",
    values_to = "Percentage"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "Normal",
        "BMD",
        "DMD"
      )
    ),
    Percentage = as.numeric(
      Percentage
    )
  )

# Select biologically informative high-abundance states

top_states <- datA %>%
  group_by(celltype) %>%
  summarise(
    mean_percentage =
      mean(
        Percentage,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(mean_percentage)
  ) %>%
  slice_head(
    n = 18
  ) %>%
  pull(celltype)

datA_plot <- datA %>%
  filter(
    celltype %in% top_states
  ) %>%
  mutate(
    celltype = factor(
      celltype,
      levels = rev(top_states)
    )
  )

pA <- ggplot(
  datA_plot,
  aes(
    x = group,
    y = celltype,
    fill = Percentage
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f",
        Percentage
      )
    ),
    size = 3.2,
    fontface = "bold"
  ) +
  scale_fill_gradient(
    low = "#F7F7F7",
    high = "#B2182B",
    name = "Cells (%)"
  ) +
  labs(
    title =
      "Global cell-state architecture",
    subtitle =
      "Cell-state composition across disease groups",
    x = NULL,
    y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x =
      element_text(
        face = "bold"
      ),
    axis.text.y =
      element_text(
        size = 10
      )
  )

############################################################
# PANEL B
# SAMPLE-LEVEL COMPOSITION
############################################################

cat("\nCreating Panel B...\n")

datB <- comp_long %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "Normal",
        "BMD",
        "DMD"
      )
    ),
    celltype =
      as.character(celltype),
    Percentage =
      as.numeric(Percentage)
  )

# Show major cell-state composition by disease group.
# This avoids inventing replicate IDs when the source table
# does not contain them.

top_B_states <- datB %>%
  group_by(celltype) %>%
  summarise(
    mean_percentage =
      mean(
        Percentage,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(mean_percentage)
  ) %>%
  slice_head(
    n = 12
  ) %>%
  pull(celltype)

datB_plot <- datB %>%
  filter(
    celltype %in% top_B_states
  )

pB <- ggplot(
  datB_plot,
  aes(
    x = group,
    y = Percentage,
    fill = celltype
  )
) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.10,
    size = 1.7,
    alpha = 0.75
  ) +
  scale_fill_manual(
    values = celltype_cols,
    drop = FALSE
  ) +
  labs(
    title =
      "Cell-state composition across groups",
    subtitle =
      "Distribution of cell-state abundance",
    x = NULL,
    y = "Composition (%)",
    fill = "Cell state"
  ) +
  theme_pub +
  theme(
    axis.text.x =
      element_text(
        face = "bold"
      ),
    legend.position =
      "right"
  )

############################################################
# PANEL C
# STATE 4 ABUNDANCE
############################################################

cat("\nCreating Panel C...\n")

state4 <- state4 %>%
  mutate(
    Group = factor(
      Group,
      levels = c(
        "Normal",
        "BMD",
        "DMD"
      )
    ),
    State4_percent =
      as.numeric(
        State4_percent
      )
  )

pC <- ggplot(
  state4,
  aes(
    x = Group,
    y = State4_percent,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    alpha = 0.65
  ) +
  geom_jitter(
    width = 0.09,
    size = 2.2,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = group_cols
  ) +
  labs(
    title =
      "State 4 abundance",
    subtitle =
      "Satellite-cell state 4 across samples",
    x = NULL,
    y = "State 4 cells (%)",
    fill = "Group"
  ) +
  theme_pub +
  theme(
    legend.position =
      "none",
    axis.text.x =
      element_text(
        face = "bold"
      )
  )

############################################################
# PANEL D
# DMD MOLECULAR SIGNATURE
############################################################

cat("\nCreating Panel D...\n")

signature <- signature %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "Normal",
        "BMD",
        "DMD"
      )
    ),
    DMD_signature_score =
      as.numeric(
        DMD_signature_score
      )
  )

pD <- ggplot(
  signature,
  aes(
    x = group,
    y = DMD_signature_score,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    alpha = 0.65
  ) +
  geom_jitter(
    width = 0.09,
    size = 2.3,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = group_cols
  ) +
  labs(
    title =
      "DMD-associated molecular signature",
    subtitle =
      "19-gene signature score",
    x = NULL,
    y = "Signature score",
    fill = "Group"
  ) +
  theme_pub +
  theme(
    legend.position =
      "none",
    axis.text.x =
      element_text(
        face = "bold"
      )
  )

############################################################
# PANEL E
# INTEGRATED MODEL
############################################################

cat("\nCreating Panel E...\n")

# Use a clean conceptual model rather than a crowded plot.

model_df <- data.frame(
  x = c(
    1,
    2,
    3,
    4
  ),
  y = c(
    2,
    2,
    2,
    2
  ),
  label = c(
    "Satellite-cell\nstate 4",
    "Altered cellular\nabundance",
    "DMD-associated\n19-gene program",
    "Disease-associated\ncell state"
  )
)

pE <- ggplot(
  model_df,
  aes(
    x = x,
    y = y
  )
) +
  geom_point(
    size = 20,
    shape = 21,
    fill = "white",
    stroke = 1.2
  ) +
  geom_text(
    aes(
      label = label
    ),
    size = 4,
    fontface = "bold"
  ) +
  geom_segment(
    data =
      data.frame(
        x = c(
          1.35,
          2.35,
          3.35
        ),
        xend = c(
          1.65,
          2.65,
          3.65
        ),
        y = c(
          2,
          2,
          2
        ),
        yend = c(
          2,
          2,
          2
        )
      ),
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    arrow =
      arrow(
        length =
          unit(
            0.18,
            "cm"
          )
      ),
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 2.5,
    y = 1.25,
    label =
      "Integrated DMD satellite-cell state model",
    size = 5,
    fontface = "bold"
  ) +
  coord_cartesian(
    xlim = c(
      0.4,
      4.6
    ),
    ylim = c(
      0.7,
      3.0
    ),
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.title =
      element_text(
        size = 17,
        face = "bold"
      ),
    plot.margin =
      margin(
        15,
        15,
        20,
        15
      )
  ) +
  labs(
    title =
      "Integrated biological model"
  )

############################################################
# 7. PANEL LABELS
############################################################

pA <- pA +
  labs(tag = "A") +
  theme(
    plot.tag =
      element_text(
        size = 20,
        face = "bold"
      )
  )

pB <- pB +
  labs(tag = "B") +
  theme(
    plot.tag =
      element_text(
        size = 20,
        face = "bold"
      )
  )

pC <- pC +
  labs(tag = "C") +
  theme(
    plot.tag =
      element_text(
        size = 20,
        face = "bold"
      )
  )

pD <- pD +
  labs(tag = "D") +
  theme(
    plot.tag =
      element_text(
        size = 20,
        face = "bold"
      )
  )

pE <- pE +
  labs(tag = "E") +
  theme(
    plot.tag =
      element_text(
        size = 20,
        face = "bold"
      )
  )

############################################################
# 8. SAVE INDIVIDUAL PANELS
############################################################

cat("\n===== SAVING INDIVIDUAL PANELS =====\n")

ggsave(
  file.path(
    out_dir,
    "Figure4A_Global_CellState_Architecture.pdf"
  ),
  pA,
  width = 7.5,
  height = 7.0
)

ggsave(
  file.path(
    out_dir,
    "Figure4A_Global_CellState_Architecture.png"
  ),
  pA,
  width = 7.5,
  height = 7.0,
  dpi = 600
)

ggsave(
  file.path(
    out_dir,
    "Figure4B_CellState_Composition.pdf"
  ),
  pB,
  width = 8.0,
  height = 7.0
)

ggsave(
  file.path(
    out_dir,
    "Figure4B_CellState_Composition.png"
  ),
  pB,
  width = 8.0,
  height = 7.0,
  dpi = 600
)

ggsave(
  file.path(
    out_dir,
    "Figure4C_State4_Abundance.pdf"
  ),
  pC,
  width = 6.5,
  height = 6.0
)

ggsave(
  file.path(
    out_dir,
    "Figure4C_State4_Abundance.png"
  ),
  pC,
  width = 6.5,
  height = 6.0,
  dpi = 600
)

ggsave(
  file.path(
    out_dir,
    "Figure4D_DMD_Signature.pdf"
  ),
  pD,
  width = 6.5,
  height = 6.0
)

ggsave(
  file.path(
    out_dir,
    "Figure4D_DMD_Signature.png"
  ),
  pD,
  width = 6.5,
  height = 6.0,
  dpi = 600
)

ggsave(
  file.path(
    out_dir,
    "Figure4E_Integrated_Model.pdf"
  ),
  pE,
  width = 11,
  height = 3.8
)

ggsave(
  file.path(
    out_dir,
    "Figure4E_Integrated_Model.png"
  ),
  pE,
  width = 11,
  height = 3.8,
  dpi = 600
)

############################################################
# 9. COMBINE FIGURE
############################################################

cat("\n===== COMBINING FIGURE 4 =====\n")

top_row <-
  pA +
  pB +
  plot_layout(
    widths = c(
      1,
      1.15
    )
  )

middle_row <-
  pC +
  pD +
  plot_layout(
    widths = c(
      1,
      1
    )
  )

final_figure <-
  top_row /
  middle_row /
  pE +
  plot_layout(
    heights = c(
      1.25,
      1,
      0.55
    )
  )

############################################################
# 10. SAVE FINAL FIGURE
############################################################

ggsave(
  file.path(
    out_dir,
    "Figure4_DMD_Satellite_CellState_PUBLICATION_v4.pdf"
  ),
  final_figure,
  width = 15,
  height = 17
)

ggsave(
  file.path(
    out_dir,
    "Figure4_DMD_Satellite_CellState_PUBLICATION_v4.png"
  ),
  final_figure,
  width = 15,
  height = 17,
  dpi = 600
)

ggsave(
  file.path(
    out_dir,
    "Figure4_DMD_Satellite_CellState_PUBLICATION_v4.tiff"
  ),
  final_figure,
  width = 15,
  height = 17,
  dpi = 600,
  compression = "lzw"
)

############################################################
# 11. SAVE RDS
############################################################

saveRDS(
  final_figure,
  file.path(
    out_dir,
    "Figure4_DMD_Satellite_CellState_PUBLICATION_v4.rds"
  )
)

############################################################
# 12. SUMMARY
############################################################

sink(
  file.path(
    out_dir,
    "Figure4_source_summary.txt"
  )
)

cat(
  "FIGURE 4 v4\n\n"
)

cat(
  "Panel A: Global cell-state architecture\n"
)

cat(
  "Panel B: Cell-state composition across groups\n"
)

cat(
  "Panel C: State 4 abundance\n"
)

cat(
  "Panel D: DMD-associated 19-gene signature\n"
)

cat(
  "Panel E: Integrated biological model\n"
)

cat(
  "\nInput files:\n"
)

cat(
  "17_CellState_Composition_Percentage.csv\n"
)

cat(
  "16_CellState_Composition_Long.csv\n"
)

cat(
  "19_CellState_Replicate_Support.csv\n"
)

cat(
  "STEP18B_State4_abundance_sample_table.csv\n"
)

cat(
  "Satellite_cluster7_DMD_19gene_signature_scores.csv\n"
)

sink()

############################################################
# 13. FINAL
############################################################

cat("\n")
cat("============================================\n")
cat("FIGURE 4 v4 COMPLETE\n")
cat("============================================\n")

cat(
  "\nOutput directory:\n",
  out_dir,
  "\n\n"
)

print(
  list.files(
    out_dir
  )
)

cat("\nDONE\n")
