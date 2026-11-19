rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(pairwiseAdonis)
library(ggplot2)
library(ggpubr)
library(cowplot)

# Run from the repository root.
input_dir <- file.path("data", "inputdata")
output_dir <- file.path("data", "output")
metadata_sheet <- "baseline"
profile_sheet <- "profile"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("Harbin_Normal", "Harbin_CRC", "SYSUCC_CRC")
group_colors <- c(Harbin_Normal = "#1B9E77", Harbin_CRC = "#3A86C8", SYSUCC_CRC = "#E16A86")
group_comparisons <- combn(group_levels, 2, simplify = FALSE)

# Read inputs and align samples.
metadata <- read_excel(file.path(input_dir, "MMR_study_metadata.xlsx"), sheet = metadata_sheet)
profile <- read_excel(file.path(input_dir, "metaphlan4_profile.xlsx"), sheet = profile_sheet)

community_matrix <- profile %>%
  pivot_longer(cols = -Species, names_to = "sampleID", values_to = "Abundance") %>%
  pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("sampleID")

metadata_clean <- metadata %>%
  filter(
    sampleID %in% rownames(community_matrix),
    type == "stool",
    cohort == "Harbin" | (cohort == "SYSUCC" & MMR_study == "Yes_primary")
  ) %>%
  mutate(cohort = factor(cohort), disease = factor(disease), Sex = factor(Sex), Age = as.numeric(Age)) %>%
  drop_na(cohort, disease, Age, Sex) %>%
  mutate(group = factor(paste(cohort, disease, sep = "_"), levels = group_levels))

stopifnot(!anyNA(metadata_clean$group), all(table(metadata_clean$group) > 0))
community_matrix <- community_matrix[metadata_clean$sampleID, , drop = FALSE]
stopifnot(identical(rownames(community_matrix), metadata_clean$sampleID))

# Bray-Curtis PERMANOVA.
bc_dist <- vegdist(community_matrix, method = "bray")

run_permanova <- function(model_formula) {
  set.seed(123)
  adonis2(model_formula, data = metadata_clean, permutations = 999, by = "margin")
}

# Test disease and cohort separately, adjusting for age and sex.
permanova_disease <- run_permanova(bc_dist ~ Age + Sex + disease)
permanova_cohort <- run_permanova(bc_dist ~ Age + Sex + cohort)

# Test disease and cohort in a joint model.
permanova_result <- run_permanova(bc_dist ~ Age + Sex + cohort + disease)

# Use the adjusted three-group test for panel A.
permanova_result_2 <- run_permanova(bc_dist ~ Age + Sex + group)

cat("\nCRC vs Normal (adjusted for age and sex)\n")
print(permanova_disease)
cat("\nSYSUCC vs Harbin (adjusted for age and sex)\n")
print(permanova_cohort)
cat("\nDisease and cohort (joint model)\n")
print(permanova_result)
cat("\nOverall three-group PERMANOVA (panel A)\n")
print(permanova_result_2)

# Pairwise PERMANOVA.
set.seed(123)
pairwise_res <- pairwise.adonis2(bc_dist ~ group, data = metadata_clean, permutations = 999)
print(pairwise_res)

# PCoA and negative eigenvalues.
pcoa_result <- cmdscale(bc_dist, k = 2, eig = TRUE)
eig_raw <- pcoa_result$eig
eig_tol <- sqrt(.Machine$double.eps) * max(abs(eig_raw))
negative_eig <- eig_raw[eig_raw < -eig_tol]
negative_fraction <- sum(abs(negative_eig)) / sum(abs(eig_raw))
message(sprintf(
  "Bray-Curtis PCoA: %d negative eigenvalues (%.2f%% of absolute inertia).",
  length(negative_eig), 100 * negative_fraction
))

# Express axis variance as a fraction of positive inertia.
pcoa_eig <- pcoa_result$eig
pcoa_axis_pct <- 100 * pcoa_eig[1:2] / sum(pcoa_eig[pcoa_eig > 0])
pcoa_scores <- as.data.frame(pcoa_result$points)
names(pcoa_scores) <- c("PCoA1", "PCoA2")
stopifnot(identical(rownames(pcoa_scores), metadata_clean$sampleID))
pcoa_data <- bind_cols(metadata_clean, pcoa_scores)

# Panel A: PCoA with marginal boxplots.
group_r2 <- permanova_result_2["group", "R2"]
group_p <- permanova_result_2["group", "Pr(>F)"]
label_text <- sprintf("R² = %.3f\nP = %.3f", group_r2, group_p)

p_main <- ggplot(pcoa_data, aes(PCoA1, PCoA2, color = group, shape = disease)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = c(CRC = 16, Normal = 17)) +
  labs(
    x = sprintf("PCoA1 (%.2f%%)", pcoa_axis_pct[1]),
    y = sprintf("PCoA2 (%.2f%%)", pcoa_axis_pct[2]),
    color = "Group", shape = "Disease"
  ) +
  annotate(
    "text", x = Inf, y = Inf, label = label_text,
    hjust = 1.1, vjust = 1.1, size = 4, fontface = "italic"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "right")

p_top <- ggplot(pcoa_data, aes(group, PCoA1, fill = group)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.2, size = 1) +
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = group_comparisons) +
  theme_void() +
  theme(legend.position = "none")

p_right <- ggplot(pcoa_data, aes(group, PCoA2, fill = group)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.2, size = 1) +
  geom_boxplot(alpha = 0.8, width = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", comparisons = group_comparisons) +
  theme_void() +
  theme(legend.position = "none")

# Align panel ranges, including significance brackets.
pcoa1_limits <- ggplot_build(p_top)$layout$panel_params[[1]]$y.range
pcoa2_limits <- ggplot_build(p_right)$layout$panel_params[[1]]$y.range
p_main <- p_main + coord_cartesian(xlim = pcoa1_limits, ylim = pcoa2_limits, expand = FALSE)
p_top <- p_top +
  scale_y_continuous(expand = expansion(mult = 0)) +
  coord_flip(ylim = pcoa1_limits)
p_right <- p_right +
  scale_y_continuous(expand = expansion(mult = 0)) +
  coord_cartesian(ylim = pcoa2_limits)

p_pcoa <- insert_xaxis_grob(p_main, p_top, grid::unit(0.2, "null"), position = "top")
p_pcoa <- insert_yaxis_grob(p_pcoa, p_right, grid::unit(0.2, "null"), position = "right")

# Panel B: alpha diversity.
alpha_df <- data.frame(
  sampleID = rownames(community_matrix),
  Shannon = diversity(community_matrix, index = "shannon"),
  Richness = rowSums(community_matrix > 0),
  Simpson = diversity(community_matrix, index = "simpson")
)

alpha_long <- metadata_clean %>%
  inner_join(alpha_df, by = "sampleID") %>%
  pivot_longer(cols = c(Shannon, Richness, Simpson), names_to = "Index", values_to = "Value") %>%
  mutate(Index = factor(Index, levels = c("Shannon", "Richness", "Simpson")))

p_alpha <- ggplot(alpha_long, aes(group, Value, fill = group)) +
  geom_jitter(width = 0.15, height = 0, size = 1, alpha = 0.4, color = "grey70") +
  geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA) +
  facet_wrap(~Index, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(
    method = "wilcox.test", label = "p.signif", comparisons = group_comparisons,
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
  ) +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Alpha Diversity", fill = "Group") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

# Save cohort alpha and beta diversity plots.
p_combined <- plot_grid(p_pcoa, p_alpha, labels = c("A", "B"), ncol = 2, rel_widths = c(1, 1))
ggsave(file.path(output_dir, "cohort_alpha_beta_diversity.pdf"), p_combined, width = 12, height = 5)
