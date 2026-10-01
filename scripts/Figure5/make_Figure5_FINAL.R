############################################################
# FIGURE 5 FINAL
# DMD PRECISION MEDICINE — 4-PANEL PUBLICATION FIGURE
#
# A. State-4 therapeutic candidate landscape
# B. Transcriptomic reversal
# C. State-4 target evidence
# D. Precision-medicine therapeutic prioritization
#
# Uses existing STEP27 / STEP50 outputs.
# NO upstream analysis is rerun.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(patchwork)
  library(grid)
})

cat("\n============================================\n")
cat("FIGURE 5 FINAL — DMD PRECISION MEDICINE\n")
cat("4-PANEL PUBLICATION FIGURE\n")
cat("============================================\n\n")

############################################################
# 1. DIRECTORIES
############################################################

base_dir <- "/home/ubuntu/scRNAseq/data/GSE288958"

out_dir <- file.path(
  base_dir,
  "Publication",
  "Figure5_FINAL"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# 2. INPUT FILES
############################################################

therapeutic_file <- file.path(
  base_dir,
  "STEP50_DMD_precision_medicine_synthesis",
  "STEP50_DMD_therapeutic_candidates.csv"
)

reversal_file <- file.path(
  base_dir,
  "STEP27_DrugReversal",
  "STEP29_perturbagen_reversal_consistency.csv"
)

############################################################
# 3. CHECK INPUTS
############################################################

cat("===== INPUT CHECK =====\n")

if (!file.exists(therapeutic_file)) {
  stop("Therapeutic candidate file not found:\n", therapeutic_file)
}

if (!file.exists(reversal_file)) {
  stop("Reversal file not found:\n", reversal_file)
}

cat("Therapeutic candidates : TRUE\n")
cat("Reversal data          : TRUE\n\n")

############################################################
# 4. READ DATA
############################################################

therapeutic <- read_csv(
  therapeutic_file,
  show_col_types = FALSE
)

reversal <- read_csv(
  reversal_file,
  show_col_types = FALSE
)

############################################################
# 5. CLEAN THERAPEUTIC DATA
############################################################

therapeutic <- therapeutic %>%
  filter(
    !is.na(compound_name),
    compound_name != "",
    compound_name != "BRD-K28366633"
  ) %>%
  filter(
    !is.na(State4_percent_cells),
    is.finite(State4_percent_cells)
  )

############################################################
# 6. KEEP RELEVANT CANDIDATES
############################################################

preferred_order <- c(
  "PD-184352",
  "manumycin A",
  "QL-XII-47"
)

therapeutic <- therapeutic %>%
  filter(compound_name %in% preferred_order) %>%
  mutate(
    compound_name = factor(
      compound_name,
      levels = preferred_order
    )
  )

############################################################
# 7. REVERSAL DATA
############################################################

reversal_plot <- reversal %>%
  filter(
    compound_name %in% preferred_order
  ) %>%
  group_by(compound_name) %>%
  summarise(
    mean_connectivity = mean(
      mean_connectivity,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    compound_name = factor(
      compound_name,
      levels = preferred_order
    )
  )

############################################################
# 8. PANEL A
# BUBBLE / LOLLIPOP THERAPEUTIC LANDSCAPE
############################################################

cat("Creating Panel A...\n")

panel_A <- ggplot(
  therapeutic,
  aes(
    x = State4_percent_cells,
    y = compound_name
  )
) +

  geom_segment(
    aes(
      x = 0,
      xend = State4_percent_cells,
      y = compound_name,
      yend = compound_name
    ),
    linewidth = 1.2,
    alpha = 0.30
  ) +

  geom_point(
    aes(
      size = State4_percent_cells
    ),
    alpha = 0.85
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.1f%%",
        State4_percent_cells
      )
    ),
    nudge_x = 1.0,
    fontface = "bold",
    size = 4
  ) +

  scale_size_continuous(
    range = c(5, 14),
    guide = "none"
  ) +

  scale_x_continuous(
    limits = c(0, 34),
    expand = expansion(
      mult = c(0, 0.05)
    )
  ) +

  labs(
    title = "A  State-4 therapeutic candidate landscape",
    subtitle =
      "Candidate-target expression within the disease-associated State-4 population",
    x = "State-4 cells expressing candidate target (%)",
    y = NULL
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 11
    ),
    axis.title.x = element_text(
      face = "bold"
    ),
    panel.grid.major.x = element_line(
      linewidth = 0.25,
      colour = "grey85"
    )
  )

############################################################
# 9. PANEL B
# TRANSCRIPTOMIC REVERSAL
############################################################

cat("Creating Panel B...\n")

panel_B <- ggplot(
  reversal_plot,
  aes(
    x = mean_connectivity,
    y = compound_name
  )
) +

  geom_vline(
    xintercept = 0,
    linewidth = 0.8,
    alpha = 0.7
  ) +

  geom_segment(
    aes(
      x = 0,
      xend = mean_connectivity,
      y = compound_name,
      yend = compound_name
    ),
    linewidth = 1.1,
    alpha = 0.30
  ) +

  geom_point(
    size = 5
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        mean_connectivity
      )
    ),
    nudge_x = -0.012,
    fontface = "bold",
    size = 3.8
  ) +

  labs(
    title = "B  Transcriptomic reversal",
    subtitle =
      "Negative connectivity indicates reversal of the DMD disease signature",
    x = NULL,
    y = NULL
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      size = 10
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(
      face = "bold",
      size = 11
    ),
    axis.line.x = element_blank(),
    panel.grid = element_blank()
  )

############################################################
# 10. PANEL C
# CIRCULAR STATE-4 TARGET EVIDENCE
############################################################

cat("Creating Panel C...\n")

circle_data <- therapeutic %>%
  mutate(
    angle = seq(
      0,
      360 - 360 / n(),
      length.out = n()
    ),
    fraction = State4_percent_cells / 30
  )

############################################################
# Build circular plot using ggplot polar coordinates
############################################################

panel_C_data <- circle_data %>%
  mutate(
    x = seq_len(n()),
    ymin = 0.35,
    ymax = 0.35 + fraction * 0.65
  )

panel_C <- ggplot(
  panel_C_data
) +

  geom_col(
    aes(
      x = x,
      y = ymax,
      fill = compound_name
    ),
    width = 0.75,
    ymin = 0.35,
    colour = "black",
    linewidth = 0.4,
    show.legend = FALSE
  ) +

  geom_text(
    aes(
      x = x,
      y = ymax + 0.08,
      label = paste0(
        compound_name,
        "\n",
        sprintf(
          "%.1f%%",
          State4_percent_cells
        )
      )
    ),
    fontface = "bold",
    size = 3.5
  ) +

  coord_polar(
    start = 0
  ) +

  ylim(
    0,
    1.25
  ) +

  labs(
    title =
      "C  State-4 target evidence",
    subtitle =
      "Radial extent represents the fraction of State-4 cells expressing candidate targets"
  ) +

  theme_void(
    base_size = 12
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 10,
      hjust = 0
    )
  )

############################################################
# 11. PANEL D
# PRECISION MEDICINE FUNNEL / CONE
############################################################

cat("Creating Panel D...\n")

############################################################
# Use grid graphics for the funnel
############################################################

make_panel_D <- function() {

  g <- grid::grobTree(

    grid::textGrob(
      "D  Precision-medicine therapeutic prioritization",
      x = 0.01,
      y = 0.97,
      just = c("left", "top"),
      gp = grid::gpar(
        fontsize = 15,
        fontface = "bold"
      )
    ),

    grid::textGrob(
      "Convergence of disease-state biology, transcriptomic reversal and target evidence",
      x = 0.01,
      y = 0.90,
      just = c("left", "top"),
      gp = grid::gpar(
        fontsize = 10
      )
    ),

    ########################################################
    # Funnel body
    ########################################################

    grid::polygonGrob(
      x = c(
        0.10, 0.90,
        0.76, 0.24
      ),
      y = c(
        0.80, 0.80,
        0.23, 0.23
      ),
      gp = grid::gpar(
        fill = "grey92",
        col = "black",
        lwd = 1.5
      )
    ),

    ########################################################
    # Horizontal levels
    ########################################################

    grid::linesGrob(
      x = c(
        0.10, 0.90,
        NA,
        0.17, 0.83,
        NA,
        0.24, 0.76,
        NA,
        0.31, 0.69
      ),
      y = c(
        0.80, 0.80,
        NA,
        0.64, 0.64,
        NA,
        0.48, 0.48,
        NA,
        0.34, 0.34
      ),
      gp = grid::gpar(
        lwd = 1
      )
    ),

    ########################################################
    # Labels
    ########################################################

    grid::textGrob(
      "DMD disease signature",
      x = 0.50,
      y = 0.74,
      gp = grid::gpar(
        fontsize = 11,
        fontface = "bold"
      )
    ),

    grid::textGrob(
      "Disease-associated State-4 program",
      x = 0.50,
      y = 0.58,
      gp = grid::gpar(
        fontsize = 10,
        fontface = "bold"
      )
    ),

    grid::textGrob(
      "Transcriptomic reversal",
      x = 0.50,
      y = 0.43,
      gp = grid::gpar(
        fontsize = 10,
        fontface = "bold"
      )
    ),

    grid::textGrob(
      "State-4 target evidence",
      x = 0.50,
      y = 0.30,
      gp = grid::gpar(
        fontsize = 10,
        fontface = "bold"
      )
    ),

    ########################################################
    # Candidate hypotheses
    ########################################################

    grid::textGrob(
      "Candidate therapeutic hypotheses",
      x = 0.50,
      y = 0.16,
      gp = grid::gpar(
        fontsize = 11,
        fontface = "bold"
      )
    ),

    ########################################################
    # Candidate nodes
    ########################################################

    grid::pointsGrob(
      x = c(
        0.28,
        0.50,
        0.72
      ),
      y = c(
        0.07,
        0.07,
        0.07
      ),
      pch = 21,
      size = grid::unit(
        0.12,
        "npc"
      ),
      gp = grid::gpar(
        fill = "white",
        col = "black",
        lwd = 1.5
      )
    ),

    grid::textGrob(
      "PD-184352",
      x = 0.28,
      y = 0.01,
      gp = grid::gpar(
        fontsize = 9,
        fontface = "bold"
      )
    ),

    grid::textGrob(
      "QL-XII-47",
      x = 0.50,
      y = 0.01,
      gp = grid::gpar(
        fontsize = 9,
        fontface = "bold"
      )
    ),

    grid::textGrob(
      "manumycin A",
      x = 0.72,
      y = 0.01,
      gp = grid::gpar(
        fontsize = 9,
        fontface = "bold"
      )
    )
  )

  return(g)
}

panel_D <- make_panel_D()

############################################################
# 12. INDIVIDUAL PANEL EXPORT
############################################################

cat("Saving individual panels...\n")

ggsave(
  file.path(
    out_dir,
    "Figure5A_State4_Therapeutic_Candidates.pdf"
  ),
  panel_A,
  width = 8,
  height = 5,
  units = "in"
)

ggsave(
  file.path(
    out_dir,
    "Figure5A_State4_Therapeutic_Candidates.png"
  ),
  panel_A,
  width = 8,
  height = 5,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(
    out_dir,
    "Figure5A_State4_Therapeutic_Candidates.tiff"
  ),
  panel_A,
  width = 8,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

############################################################

ggsave(
  file.path(
    out_dir,
    "Figure5B_Transcriptomic_Reversal.pdf"
  ),
  panel_B,
  width = 8,
  height = 5,
  units = "in"
)

ggsave(
  file.path(
    out_dir,
    "Figure5B_Transcriptomic_Reversal.png"
  ),
  panel_B,
  width = 8,
  height = 5,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(
    out_dir,
    "Figure5B_Transcriptomic_Reversal.tiff"
  ),
  panel_B,
  width = 8,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

############################################################

ggsave(
  file.path(
    out_dir,
    "Figure5C_State4_Target_Evidence.pdf"
  ),
  panel_C,
  width = 7,
  height = 6,
  units = "in"
)

ggsave(
  file.path(
    out_dir,
    "Figure5C_State4_Target_Evidence.png"
  ),
  panel_C,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(
    out_dir,
    "Figure5C_State4_Target_Evidence.tiff"
  ),
  panel_C,
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

############################################################
# Panel D export
############################################################

pdf(
  file.path(
    out_dir,
    "Figure5D_Precision_Medicine_Framework.pdf"
  ),
  width = 8.5,
  height = 5.8
)

grid::grid.newpage()
grid::grid.draw(panel_D)

dev.off()

############################################################

png(
  file.path(
    out_dir,
    "Figure5D_Precision_Medicine_Framework.png"
  ),
  width = 5100,
  height = 3480,
  res = 600
)

grid::grid.newpage()
grid::grid.draw(panel_D)

dev.off()

############################################################

tiff(
  file.path(
    out_dir,
    "Figure5D_Precision_Medicine_Framework.tiff"
  ),
  width = 5100,
  height = 3480,
  res = 600,
  compression = "lzw"
)

grid::grid.newpage()
grid::grid.draw(panel_D)

dev.off()

############################################################
# 13. COMBINED FIGURE
############################################################

cat("Creating combined Figure 5...\n")

combined_plot <-
  (panel_A | panel_B) /
  (panel_C | patchwork::wrap_elements(panel_D)) +

  plot_annotation(
    title =
      "DMD precision medicine: disease-state-informed therapeutic prioritization",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 20,
        hjust = 0.5
      )
    )
  )

############################################################
# 14. SAVE COMBINED FIGURE
############################################################

combined_pdf <- file.path(
  out_dir,
  "Figure5_DMD_PrecisionMedicine_FINAL.pdf"
)

combined_png <- file.path(
  out_dir,
  "Figure5_DMD_PrecisionMedicine_FINAL.png"
)

combined_tiff <- file.path(
  out_dir,
  "Figure5_DMD_PrecisionMedicine_FINAL.tiff"
)

ggsave(
  combined_pdf,
  combined_plot,
  width = 14,
  height = 10,
  units = "in"
)

ggsave(
  combined_png,
  combined_plot,
  width = 14,
  height = 10,
  units = "in",
  dpi = 600
)

ggsave(
  combined_tiff,
  combined_plot,
  width = 14,
  height = 10,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

############################################################
# 15. SOURCE SUMMARY
############################################################

summary_file <- file.path(
  out_dir,
  "Figure5_source_summary.txt"
)

sink(summary_file)

cat("FIGURE 5 FINAL — SOURCE SUMMARY\n")
cat("============================================\n\n")

cat("Source files:\n")
cat(therapeutic_file, "\n")
cat(reversal_file, "\n\n")

cat("Therapeutic candidates represented:\n")

print(
  therapeutic %>%
    select(
      compound_name,
      State4_percent_cells
    )
)

cat("\nMean transcriptomic connectivity:\n")

print(
  reversal_plot
)

cat("\nInterpretation used in the figure:\n")
cat(
  "A: State-4 candidate-target expression landscape.\n"
)

cat(
  "B: Negative transcriptomic connectivity represents reversal of the DMD disease signature.\n"
)

cat(
  "C: Radial extent represents State-4 target-expression prevalence.\n"
)

cat(
  "D: Funnel representation of disease signature -> State-4 program -> "
)

cat(
  "transcriptomic reversal -> target evidence -> candidate therapeutic hypotheses.\n"
)

sink()

############################################################
# 16. FINAL CHECK
############################################################

cat("\n============================================\n")
cat("FIGURE 5 FINAL COMPLETE\n")
cat("============================================\n\n")

cat("Output directory:\n")
cat(out_dir, "\n\n")

print(
  list.files(
    out_dir
  )
)

cat("\nDONE\n")

