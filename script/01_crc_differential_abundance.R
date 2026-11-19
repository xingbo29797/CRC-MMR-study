rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(Maaslin2)
library(MicrobiomeStat)
library(ANCOMBC)
library(ggplot2)
library(ggrepel)
library(eulerr)

# Run from the repository root.
input_dir <- file.path("data", "inputdata")
output_dir <- file.path("data", "output")
metadata_sheet <- "baseline"
profile_sheet <- "profile"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read inputs.
metadata <- read_excel(file.path(input_dir, "MMR_study_metadata.xlsx"), sheet = metadata_sheet)
profile <- read_excel(file.path(input_dir, "metaphlan4_profile.xlsx"), sheet = profile_sheet)

# Read non-human read counts.
file1_path <- file.path(input_dir, "filter_summary_n2578.txt")
file2_path <- file.path(input_dir, "filter_summary_n2432.txt")

summary_data1 <- read_tsv(
  file1_path, col_types = cols(.default = col_guess(), Library_ID = col_character(),
                              Nohomo_reads_count = col_double()),
  show_col_types = FALSE
)
summary_data2 <- read_tsv(
  file2_path, col_types = cols(.default = col_guess(), Library_ID = col_character(),
                              Nohomo_reads_count = col_double()),
  show_col_types = FALSE
)

summary_all <- bind_rows(summary_data1, summary_data2) %>%
  dplyr::select(Library_ID, Nohomo_reads_count) %>%
  distinct(Library_ID, .keep_all = TRUE)

# Standardize feature names.
standardize_feature <- function(x, dataset_name) {
  original <- as.character(x)
  standardized <- make.names(original)
  duplicated_names <- unique(standardized[
    duplicated(standardized) | duplicated(standardized, fromLast = TRUE)
  ])
  if (length(duplicated_names) > 0) {
    stop(dataset_name, ": duplicate names after make.names(): ", paste(duplicated_names, collapse = ", "))
  }
  standardized
}

# Compare CRC cases with Harbin controls.
run_crc_comparison <- function(case_cohort, output_prefix, metadata, profile, summary_all,
                               output_dir, min_prevalence = 0.10) {
  message("\nStarting: ", output_prefix)

  # Select cases; restrict SYSUCC to Yes_primary.
  case_metadata <- metadata %>%
    filter(cohort == case_cohort, disease == "CRC", type == "stool")

  if (case_cohort == "SYSUCC") {
    case_metadata <- case_metadata %>%
      filter(MMR_study == "Yes_primary")
  }

  # Use Harbin Normal as the control group.
  control_metadata <- metadata %>%
    filter(cohort == "Harbin", disease == "Normal", type == "stool")

  metadata_clean <- bind_rows(case_metadata, control_metadata) %>%
    distinct(sampleID, .keep_all = TRUE) %>%
    mutate(
      disease = factor(disease, levels = c("Normal", "CRC")),
      mmr_status = as.factor(mmr_status), Sex = as.factor(Sex), Age = as.numeric(Age),
      BMI = as.numeric(BMI), Tumor_Location = as.factor(Tumor_Location)
    )

  if (nrow(case_metadata) == 0 || nrow(control_metadata) == 0) {
    stop(output_prefix, ": no cases or Harbin Normal controls.")
  }
  if (anyDuplicated(metadata_clean$sampleID)) {
    stop(output_prefix, ": duplicate sample IDs in metadata.")
  }

  # Build the abundance matrix.
  profile_long <- profile %>%
    pivot_longer(cols = -Species, names_to = "sampleID", values_to = "Abundance")

  community_matrix_all <- profile_long %>%
    filter(sampleID %in% metadata_clean$sampleID) %>%
    pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
    column_to_rownames("sampleID") %>%
    as.matrix() / 100

  # Filter species by prevalence in CRC cases.
  case_sample_ids <- case_metadata$sampleID

  case_matrix <- profile_long %>%
    filter(sampleID %in% case_sample_ids) %>%
    pivot_wider(names_from = Species, values_from = Abundance, values_fill = 0) %>%
    column_to_rownames("sampleID") %>%
    as.matrix() / 100

  min_presence <- ceiling(nrow(case_matrix) * min_prevalence)
  core_species <- colnames(case_matrix)[colSums(case_matrix > 0) >= min_presence]

  if (length(core_species) == 0) {
    stop(output_prefix, ": no species meet the CRC prevalence threshold.")
  }

  community_matrix_filtered <- community_matrix_all[, core_species, drop = FALSE]
  metadata_clean <- metadata_clean %>%
    filter(sampleID %in% rownames(community_matrix_filtered))
  community_matrix_filtered <- community_matrix_filtered[metadata_clean$sampleID, , drop = FALSE]

  df_features <- as.matrix(community_matrix_filtered)
  df_metadata <- metadata_clean %>%
    column_to_rownames("sampleID")

  stopifnot(identical(rownames(df_features), rownames(df_metadata)))

  message(output_prefix, ": cases = ", nrow(case_metadata),
          "; controls = ", nrow(control_metadata),
          "; core species = ", length(core_species))

  # MaAsLin2.
  maaslin_output <- file.path(output_dir, paste0("maaslin2_", output_prefix))
  set.seed(123)

  logging::logReset()
  logging::basicConfig(level = "FINEST")

  # Validate model results if later MaAsLin2 steps fail.
  maaslin_result_file <- file.path(maaslin_output, "all_results.tsv")

  # Remove stale results before running MaAsLin2.
  if (file.exists(maaslin_result_file)) {
    file.remove(maaslin_result_file)
  }

  maaslin_error <- NULL
  fit_output <- tryCatch(
    Maaslin2(
      input_data = df_features,
      input_metadata = df_metadata,
      output = maaslin_output,
      fixed_effects = c("Age", "Sex", "disease"),
      reference = c("disease,Normal"),
      normalization = "TSS",
      transform = "LOG",
      correction = "BH",
      max_pngs = 10,
      min_abundance = 0,
      min_prevalence = 0
    ),
    error = function(e) {
      maaslin_error <<- conditionMessage(e)
      NULL
    }
  )

  # Check results before continuing.
  if (!file.exists(maaslin_result_file)) {
    if (!is.null(maaslin_error)) {
      stop(output_prefix, ": MaAsLin2 failed without all_results.tsv: ", maaslin_error)
    } else {
      stop(output_prefix, ": missing MaAsLin2 results: ", maaslin_result_file)
    }
  }

  maaslin_result_check <- tryCatch(
    read_tsv(maaslin_result_file, show_col_types = FALSE),
    error = function(e) {
      stop(output_prefix, ": unreadable MaAsLin2 results: ", conditionMessage(e))
    }
  )

  required_maaslin_cols <- c("feature", "metadata", "value", "coef", "stderr", "pval", "qval")
  if (nrow(maaslin_result_check) == 0 ||
      !all(required_maaslin_cols %in% names(maaslin_result_check))) {
    stop(output_prefix, ": empty or incomplete MaAsLin2 results.")
  }

  if (!is.null(maaslin_error)) {
    warning(
      output_prefix,
      ": validated MaAsLin2 results; continuing after a later error: ",
      maaslin_error
    )
  }

  # LinDA.
  df_metadata_linda <- df_metadata[rownames(df_features), , drop = FALSE]
  stopifnot(identical(rownames(df_features), rownames(df_metadata_linda)))
  df_metadata_linda$disease <- relevel(factor(df_metadata_linda$disease), ref = "Normal")

  set.seed(123)
  linda_result <- linda(
    feature.dat = t(df_features),
    meta.dat = df_metadata_linda,
    formula = "~ disease + Age + Sex",
    feature.dat.type = "proportion",
    p.adj.method = "BH",
    alpha = 0.05
  )

  linda_result_df <- linda_result$output$diseaseCRC %>%
    rownames_to_column("Species")

  linda_file <- file.path(output_dir, paste0("linda_", output_prefix, ".csv"))
  write.csv(linda_result_df, linda_file, row.names = FALSE)

  # ANCOM-BC2.
  current_samples <- rownames(df_metadata)

  # Match non-human read counts by sample ID.
  valid_samples <- intersect(current_samples, summary_all$Library_ID)

  reads_lookup <- summary_all %>%
    filter(Library_ID %in% valid_samples) %>%
    slice(match(valid_samples, Library_ID))

  reads_vector <- reads_lookup$Nohomo_reads_count

  rel_matrix_valid <- community_matrix_filtered[valid_samples, , drop = FALSE]
  stopifnot(identical(rownames(rel_matrix_valid), reads_lookup$Library_ID))
  stopifnot(all(is.finite(reads_vector)), all(reads_vector > 0))
  # Scale each sample by its matched non-human read count.
  abs_matrix_valid <- sweep(rel_matrix_valid, MARGIN = 1, STATS = reads_vector, FUN = "*")
  abs_matrix_integer <- round(abs_matrix_valid)

  df_features_abs <- abs_matrix_integer
  df_metadata_abs <- df_metadata[valid_samples, , drop = FALSE]
  df_metadata_abs$disease <- relevel(factor(df_metadata_abs$disease), ref = "Normal")
  stopifnot(identical(rownames(df_features_abs), rownames(df_metadata_abs)))

  res_ancombc2 <- ancombc2(
    data = t(df_features_abs),
    meta_data = df_metadata_abs,
    assay_name = "counts",
    fix_formula = "disease + Age + Sex",
    p_adj_method = "holm",
    pseudo_sens = TRUE,
    prv_cut = 0,
    lib_cut = 0,
    struc_zero = FALSE,
    neg_lb = FALSE,
    alpha = 0.05,
    verbose = TRUE
  )

  ancom_result_df <- res_ancombc2$res
  ancom_file <- file.path(output_dir, paste0("ancombc2_", output_prefix, ".csv"))
  write.csv(ancom_result_df, ancom_file, row.names = TRUE)

  # Combine differential abundance results.
  maaslin_result <- maaslin_result_check
  linda_result_read <- read_csv(linda_file, show_col_types = FALSE)
  ancom_result_read <- read_csv(ancom_file, show_col_types = FALSE)

  df_maas_filtered <- maaslin_result %>%
    filter(metadata == "disease", value == "CRC")

  df_maas <- df_maas_filtered %>%
    transmute(
      Feature = standardize_feature(feature, "MaAsLin2"),
      Maaslin_coef = coef,
      Maaslin_stderr = stderr,
      Maaslin_pval = pval,
      Maaslin_qval = qval
    )

  linda_feature_col <- names(linda_result_read)[1]
  df_linda <- linda_result_read %>%
    mutate(Feature = standardize_feature(.data[[linda_feature_col]], "LinDA")) %>%
    dplyr::select(Feature, log2FoldChange, pvalue, padj) %>%
    rename(
      Linda_log2FC = log2FoldChange, Linda_pval = pvalue, Linda_qval = padj
    )

  if (!"taxon" %in% names(ancom_result_read)) {
    stop(output_prefix, ": missing taxon column in ANCOM-BC2 results.")
  }

  df_ancom <- ancom_result_read %>%
    mutate(Feature = standardize_feature(taxon, "ANCOM-BC2")) %>%
    dplyr::select(Feature,
                  any_of(c("lfc_diseaseCRC", "p_diseaseCRC", "q_diseaseCRC", "passed_ss_diseaseCRC"))) %>%
    rename(
      Ancom_lfc = any_of("lfc_diseaseCRC"),
      Ancom_pval = any_of("p_diseaseCRC"),
      Ancom_qval = any_of("q_diseaseCRC"),
      Ancom_passed = any_of("passed_ss_diseaseCRC")
    )

  # Reject duplicate feature names before joining.
  if (anyDuplicated(df_maas$Feature) ||
      anyDuplicated(df_linda$Feature) ||
      anyDuplicated(df_ancom$Feature)) {
    stop(output_prefix, ": duplicate standardized features; cannot join results.")
  }

  merged_result <- df_maas %>%
    full_join(df_linda, by = "Feature") %>%
    full_join(df_ancom, by = "Feature") %>%
    mutate(
      Maaslin_CI_lower = Maaslin_coef - 1.96 * Maaslin_stderr,
      Maaslin_CI_upper = Maaslin_coef + 1.96 * Maaslin_stderr,
      Linda_direction = (Maaslin_coef * Linda_log2FC) > 0,
      Ancom_direction = (Maaslin_coef * Ancom_lfc) > 0,
      final = if_else(
        Maaslin_qval < 0.05 & abs(Maaslin_coef) >= 1 &
          ((Linda_qval < 0.05 & Linda_direction == TRUE) |
             (Ancom_qval < 0.05 & Ancom_direction == TRUE & Ancom_passed == TRUE)),
        1,
        0
      )
    )

  result_file <- file.path(output_dir, paste0("differential_abundance_", output_prefix, ".csv"))
  write.csv(merged_result, result_file, row.names = FALSE)

  counts_summary <- merged_result %>%
    summarise(
      Maaslin_filtered_count = sum(abs(Maaslin_coef) >= 1 & Maaslin_qval < 0.05, na.rm = TRUE),
      Maaslin_up_count = sum(Maaslin_coef >= 1 & Maaslin_qval < 0.05, na.rm = TRUE),
      Maaslin_down_count = sum(Maaslin_coef <= -1 & Maaslin_qval < 0.05, na.rm = TRUE),
      Final_is_1_count = sum(final == 1, na.rm = TRUE),
      Final_up_count = sum(final == 1 & Maaslin_coef >= 1, na.rm = TRUE),
      Final_down_count = sum(final == 1 & Maaslin_coef <= -1, na.rm = TRUE)
    )
  print(counts_summary)

  message("Completed: ", output_prefix)

  invisible(list(
    metadata = df_metadata,
    features = df_features,
    core_species = core_species,
    result = merged_result,
    counts_summary = counts_summary
  ))
}

# Run both cohort comparisons.
# SYSUCC CRC vs Harbin Normal.
result_sysucc <- run_crc_comparison(
  case_cohort = "SYSUCC",
  output_prefix = "sysucc_crc_vs_harbin_normal",
  metadata = metadata,
  profile = profile,
  summary_all = summary_all, output_dir = output_dir
)

# Harbin CRC vs Harbin Normal.
result_harbin <- run_crc_comparison(
  case_cohort = "Harbin",
  output_prefix = "harbin_crc_vs_harbin_normal",
  metadata = metadata,
  profile = profile,
  summary_all = summary_all, output_dir = output_dir
)

# Compare differential species across cohorts.
df_sysucc <- read.csv(file.path(output_dir, "differential_abundance_sysucc_crc_vs_harbin_normal.csv"), stringsAsFactors = FALSE)
df_harbin <- read.csv(file.path(output_dir, "differential_abundance_harbin_crc_vs_harbin_normal.csv"), stringsAsFactors = FALSE)

sysucc_filtered <- df_sysucc %>%
  filter(final == 1) %>%
  dplyr::select(Feature,
                Maaslin_coef_sysucc = Maaslin_coef,
                Maaslin_qval_sysucc = Maaslin_qval)

harbin_filtered <- df_harbin %>%
  filter(final == 1) %>%
  dplyr::select(Feature,
                Maaslin_coef_harbin = Maaslin_coef,
                Maaslin_qval_harbin = Maaslin_qval)

merged_comparison <- full_join(sysucc_filtered, harbin_filtered, by = "Feature") %>%
  mutate(
    In_SYSUCC = !is.na(Maaslin_coef_sysucc),
    In_Harbin = !is.na(Maaslin_coef_harbin),
    Direction_Consistent = if_else(
      In_SYSUCC & In_Harbin,
      (Maaslin_coef_sysucc * Maaslin_coef_harbin) > 0,
      NA
    ),
    Category = case_when(
      In_SYSUCC & In_Harbin & Direction_Consistent == TRUE ~ "Both (Consistent)",
      In_SYSUCC & In_Harbin & Direction_Consistent == FALSE ~ "Both (Inconsistent Direction)",
      In_SYSUCC & !In_Harbin ~ "SYSUCC_only",
      !In_SYSUCC & In_Harbin ~ "Harbin_only"
    )
  )

write.csv(merged_comparison,
          file.path(output_dir, "crc_differential_abundance_comparison.csv"),
          row.names = FALSE)

both_shared <- sum(merged_comparison$Category == "Both (Consistent)", na.rm = TRUE)
only_sysucc <- sum(merged_comparison$Category == "SYSUCC_only", na.rm = TRUE) +
  sum(merged_comparison$Category == "Both (Inconsistent Direction)", na.rm = TRUE)
only_harbin <- sum(merged_comparison$Category == "Harbin_only", na.rm = TRUE) +
  sum(merged_comparison$Category == "Both (Inconsistent Direction)", na.rm = TRUE)

euler_data <- c(
  "SYSUCC_CRC" = only_sysucc,
  "Harbin_CRC" = only_harbin,
  "SYSUCC_CRC&Harbin_CRC" = both_shared
)

pdf(file.path(output_dir, "crc_differential_species_overlap.pdf"), width = 5, height = 5)
fit <- euler(euler_data)
p_veen <- plot(
  fit,
  fills = c("#E16A86", "#3A86C8"),
  alpha = 0.55,
  labels = c("SYSUCC_CRC", "Harbin_CRC"),
  quantities = list(fontsize = 15, fontface = "bold"),
  edges = list(col = "white", lwd = 1.5)
)
print(p_veen)
dev.off()

# Correlate effect estimates for concordant shared species.
maaslin2_coef <- read_csv(file.path(output_dir, "crc_differential_abundance_comparison.csv"), show_col_types = FALSE)

metadata_plot <- metadata %>%
  filter(type == "stool") %>%
  filter(
    (cohort == "Harbin" & disease == "CRC") |
      (cohort == "SYSUCC" & MMR_study == "Yes_primary")
  )
sample_ids <- metadata_plot$sampleID

plot_df <- maaslin2_coef %>%
  filter(Category == "Both (Consistent)")

target_species <- plot_df$Feature

prevalence_df <- profile %>%
  filter(Species %in% target_species) %>%
  pivot_longer(cols = -Species, names_to = "sampleID", values_to = "Abundance") %>%
  filter(sampleID %in% sample_ids) %>%
  group_by(Species) %>%
  summarise(prevalence = mean(Abundance > 0), .groups = "drop")

plot_df <- plot_df %>%
  left_join(prevalence_df, by = c("Feature" = "Species")) %>%
  mutate(prevalence = replace_na(prevalence, 0))

write.csv(plot_df, file.path(output_dir, "crc_shared_species_coefficients.csv"), row.names = FALSE)

cor_res <- cor.test(
  plot_df$Maaslin_coef_sysucc,
  plot_df$Maaslin_coef_harbin,
  method = "spearman",
  exact = FALSE
)

rho <- unname(cor_res$estimate)
pval <- cor_res$p.value
cor_label <- paste0(
  "Spearman rho = ", sprintf("%.3f", rho),
  ", P = ", format.pval(pval, digits = 2, eps = 2e-16)
)

plot_df <- plot_df %>%
  mutate(
    mean_coef_abs = (abs(Maaslin_coef_harbin) + abs(Maaslin_coef_sysucc)) / 2,
    mean_coef = (Maaslin_coef_harbin + Maaslin_coef_sysucc) / 2
  ) %>%
  group_by(direction_group = if_else(mean_coef > 0, "Positive", "Negative")) %>%
  mutate(rank_in_group = rank(-mean_coef_abs, ties.method = "first")) %>%
  ungroup() %>%
  mutate(
    label = case_when(
      mean_coef > 0 & rank_in_group <= 20 ~ Feature,
      mean_coef < 0 & rank_in_group <= 10 ~ Feature,
      TRUE ~ NA_character_
    ),
    coef_direction = case_when(
      Maaslin_coef_harbin > 0 & Maaslin_coef_sysucc > 0 ~ "Positive",
      Maaslin_coef_harbin < 0 & Maaslin_coef_sysucc < 0 ~ "Negative",
      TRUE ~ "Mixed"
    )
  )

pdf(file.path(output_dir, "crc_shared_species_coefficient_correlation.pdf"), width = 11, height = 10)
p <- ggplot(plot_df, aes(x = Maaslin_coef_harbin, y = Maaslin_coef_sysucc)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.8, color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.8, color = "grey40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8, color = "grey50") +
  geom_point(aes(size = prevalence, color = coef_direction), alpha = 0.85) +
  scale_color_manual(
    name = "Direction",
    values = c(Positive = "#D95F02", Negative = "#1B9E77", Mixed = "grey70")
  ) +
  scale_size_continuous(range = c(2, 9), name = "Prevalence") +
  geom_text_repel(
    aes(label = label),
    size = 4.5,
    fontface = "italic",
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    na.rm = TRUE
  ) +
  annotate(
    "text", x = -Inf, y = Inf, label = cor_label,
    fontface = "italic", hjust = -0.05, vjust = 1.5, size = 5.5
  ) +
  labs(
    x = "coef in Harbin_CRC vs Harbin_Normal",
    y = "coef in SYSUCC_CRC vs Harbin_Normal"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    legend.position = "right",
    legend.box = "vertical"
  )

print(p)
dev.off()
