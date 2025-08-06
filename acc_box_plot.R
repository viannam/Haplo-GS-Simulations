rm(list = ls())
library(ggplot2)
library(dplyr)
library(purrr)

# Define SNP densities and model types
snp_densities <- c("50", "100", "500", "1000", "5000", "10000", "15000")
models <- c("SNP", "Haplo", "FW", "SW", "LD")

# Pattern naming function
get_pattern <- function(model, density) {
  if (model == "SNP" & density == "50") {
    return("^Pred_SNP_50_SNPs_")
  }
  if (model == "SNP" & density == "100") {
    return("^Pred_SNP_100_SNPs_")
  }
  if (model == "SNP" & density == "500") {
    return("^Pred_SNP_500_SNPs_")
  }
  if (model == "Haplo" & density == "500") {
    return("^Pred_Haplo_500_SNPs_")
  }
  prefix <- switch(model,
                   SNP   = paste0("^Pred_SNP_", density, "_"),
                   Haplo = paste0("^Pred_Haplo_", density, "_"),
                   FW    = paste0("^Pred_FW_", density, "_"),
                   SW    = paste0("^Pred_SW_", density, "_"),
                   LD    = paste0("^Pred_LD_", density, "_")
  )
  return(prefix)
}

# Build file path lists dynamically
base_path <- "/blue/mresende/share/viannam/Simulations/output"
file_list <- list()

for (density in snp_densities) {
  dir_path <- file.path(base_path, paste0(density, "_SNPs"))
  for (model in models) {
    pattern <- get_pattern(model, density)
    files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
    file_list[[paste(model, density, sep = "_")]] <- list(files = files, model = model, density = density)
  }
}

# Function to read and process prediction files
read_pred_files <- function(files, model_name, snp_density) {
  all_data <- lapply(files, function(file) {
    df <- read.table(file, header = TRUE)
    mean_trait1 <- mean(df$Trait1)
    replicate <- as.numeric(gsub("[^0-9]", "", basename(file)))
    data.frame(Mean_Trait1 = mean_trait1, Replicate = replicate)
  })
  final_data <- do.call(rbind, all_data)
  final_data$Model <- model_name
  final_data$SNP_Density <- snp_density
  return(final_data)
}

# Apply read_pred_files to all groups in the list
all_combined <- map_dfr(file_list, ~ read_pred_files(.x$files, .x$model, .x$density))

# Compute overall mean accuracy
overall_means <- all_combined %>%
  group_by(Model, SNP_Density) %>%
  summarize(Overall_Mean = mean(Mean_Trait1), .groups = "drop")

# Convert to factor with desired order
all_combined$SNP_Density <- factor(all_combined$SNP_Density, levels = c("50", "100", "500", "1000", "5000", "10000", "15000"))
all_combined$Model <- factor(all_combined$Model, levels = c("SNP", "Haplo", "FW", "SW", "LD"))

# Final polished plot
accuracy_plot <- ggplot(all_combined, aes(x = SNP_Density, y = Mean_Trait1, fill = Model)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.8, width = 0.6, outlier.size = 0.8, outlier.shape = 21, outlier.color = "black") +
  geom_jitter(aes(color = Model), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), size = 0.7, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = overall_means, aes(x = SNP_Density, y = Overall_Mean, label = sprintf("%.3f", Overall_Mean), group = Model), position = position_dodge2(width = 0.75, preserve = "single"), vjust = -0.8, size = 3, color = "black") +
  scale_fill_manual(values = c("SNP" = "#0072B2", "Haplo" = "#E69F00", "FW" = "#56B4E9", "SW" = "#009E73", "LD" = "#D55E00")) +
  scale_color_manual(values = c("SNP" = "#0072B2", "Haplo" = "#E69F00", "FW" = "#56B4E9", "SW" = "#009E73", "LD" = "#D55E00")) +
  theme_bw(base_size = 18) +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12),
    axis.text = element_text(color = "black", size = 15),
    axis.title = element_text(face = "bold", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Prediction Accuracy by SNP Density and Model",
    x = "SNP Density",
    y = "Prediction Accuracy (Trait1)",
    fill = "GS Model"
  ) +
  ylim(0, 1)

# Save as high-quality TIFF
ggsave(
  filename = "GS_accuracy_all_densities.tiff",
  plot = accuracy_plot,
  width = 18,
  height = 10,
  dpi = 600,
  units = "in",
  compression = "lzw"
)

print(accuracy_plot)
