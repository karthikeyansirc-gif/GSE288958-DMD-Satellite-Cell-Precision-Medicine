# ============================================================
# Figure 2D — State 4 proportion among satellite cells
# Publication-ready corrected version
# ============================================================

library(ggplot2)
library(dplyr)
library(readr)

# ------------------------------------------------------------
# 1. Input
# ------------------------------------------------------------

input_file <- "STEP18B_State4_abundance_sample_table.csv"

dat <- read_csv(input_file, show_col_types = FALSE)

cat("\n===== COLUMN NAMES =====\n")
print(names(dat))

cat("\n===== FIRST ROWS =====\n")
print(head(dat))

# ------------------------------------------------------------
# 2. Identify sample/group/abundance columns
# ------------------------------------------------------------

sample_col <- names(dat)[
  grepl("sample", names(dat), ignore.case = TRUE)
][1]

group_col <- names(dat)[
  grepl("group|condition|diagnosis|disease",
        names(dat), ignore.case = TRUE)
][1]

value_candidates <- names(dat)[
  grepl("state.?4|abundance|proportion|percent|percentage",
        names(dat), ignore.case = TRUE)
]

cat("\nDetected sample column:", sample_col, "\n")
cat("Detected group column:", group_col, "\n")
cat("Possible State-4 columns:\n")
print(value_candidates)

# ------------------------------------------------------------
# 3. IMPORTANT:
# Select the actual State-4 abundance column
# ------------------------------------------------------------

value_col <- "State4_percent"
cat("\nUsing verified abundance column:", value_col, "\n")


# ------------------------------------------------------------
# 4. Prepare data
# ------------------------------------------------------------

plot_dat <- dat %>%
  select(
    Sample = all_of(sample_col),
    Group = all_of(group_col),
    State4 = all_of(value_col)
  ) %>%
  filter(
    !is.na(Sample),
    !is.na(Group),
    !is.na(State4)
  )

# ------------------------------------------------------------
# 5. Check the values BEFORE plotting
# ------------------------------------------------------------

cat("\n===== RAW STATE 4 VALUES =====\n")
print(plot_dat)

cat("\n===== RANGE =====\n")
print(range(plot_dat$State4, na.rm = TRUE))

# ------------------------------------------------------------
# 6. Convert percentage to proportion if necessary
# ------------------------------------------------------------


plot_dat <- plot_dat %>%
  mutate(State4 = State4 / 100)

cat("\n===== STATE 4 PROPORTION =====\n")
print(plot_dat)

cat("\n===== PROPORTION RANGE =====\n")
print(range(plot_dat$State4, na.rm = TRUE))

# ------------------------------------------------------------
# 7. Keep biologically valid range
# ------------------------------------------------------------

if (any(plot_dat$State4 < 0 | plot_dat$State4 > 1)) {

  warning(
    "Some State 4 values are outside the valid 0–1 range. ",
    "Check the input table before publication."
  )

}

# ------------------------------------------------------------
# 8. Order groups
# ------------------------------------------------------------

plot_dat$Group <- factor(
  plot_dat$Group,
  levels = c("Normal", "BMD", "DMD")
)

# ------------------------------------------------------------
# 9. Publication plot
# ------------------------------------------------------------

p <- ggplot(
  plot_dat,
  aes(
    x = Group,
    y = State4
  )
) +

  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ) +

  geom_jitter(
    width = 0.08,
    size = 3,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.7
  ) +

  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = function(x) sprintf("%.1f", x)
  ) +

  labs(
    x = NULL,
    y = "State 4 proportion among satellite cells"
  ) +

  theme_classic(base_size = 12) +

  theme(
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 10),
    axis.line = element_line(linewidth = 0.5),
    plot.margin = margin(8, 8, 8, 8)
  )

# ------------------------------------------------------------
# 10. Output directory
# ------------------------------------------------------------

outdir <- "Publication/Figure2_v2"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 11. Save PDF
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    outdir,
    "Figure2D_State4_proportion.pdf"
  ),
  plot = p,
  width = 4.2,
  height = 4.2,
  units = "in"
)

# ------------------------------------------------------------
# 12. Save PNG
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    outdir,
    "Figure2D_State4_proportion.png"
  ),
  plot = p,
  width = 4.2,
  height = 4.2,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------
# 13. Save the plotting data
# ------------------------------------------------------------

write.csv(
  plot_dat,
  file.path(
    outdir,
    "Figure2D_State4_proportion_plot_data.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("Figure 2D generated successfully.\n")
cat("========================================\n")
