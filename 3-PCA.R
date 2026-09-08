# ==============================================================
# PCA for optimal-scale (EOS)
# Input: All_species_EOS_for_PCA.xlsx
# Output: PCA scatter plot, loadings heatmap, species-score heatmap, csv result tables
# ==============================================================
library(readxl)
library(ggplot2)
library(pheatmap)

# --------------------------
# 1. Set working directory and import data
setwd(".")
df_scale_raw <- read_excel("All_species_EOS_for_PCA.xlsx")

# --------------------------
# 2. Data preprocessing
spe_labels <- df_scale_raw[["spe"]]
spe_pca_input <- df_scale_raw[, !colnames(df_scale_raw) %in% c("spe")]
rownames(spe_pca_input) <- spe_labels

# --------------------------
# 3. Perform PCA
# center + scale handled by prcomp
pca_result <- prcomp(spe_pca_input, center = TRUE, scale. = TRUE)
pca_summary <- summary(pca_result)

pca_scores      <- as.data.frame(pca_result$x)
pca_loadings    <- pca_result$rotation
pca_scores_top4 <- pca_scores[, 1:4]

PC1_lab <- paste0("PC1 (", round(pca_summary$importance[2,1]*100, 2), "%)")
PC2_lab <- paste0("PC2 (", round(pca_summary$importance[2,2]*100, 2), "%)")

# --------------------------
# 4. PCA Score Plot
png("PCA_Scale_Plot.png", width = 800, height = 600, bg = "white")
# Fix: assign ggplot output to an object p_pca
p_pca <- ggplot(data = pca_scores,
                aes(x = PC1, y = PC2, color = rownames(spe_pca_input))) +
  coord_cartesian(xlim = c(-7, 7), ylim = c(-7, 7)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_point(size = 3.5) +
  labs(x = PC1_lab, y = PC2_lab, color = "Species", title = "PCA Score Plot") +
  guides(fill = "none") +
  theme_bw() +
  scale_fill_manual(values = c("purple","orange","pink","blue","lightgreen","red","lightblue","darkgreen")) +
  scale_colour_manual(values = c("purple","orange","pink","blue","lightgreen","red","lightblue","darkgreen")) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),'cm'))

print(p_pca)  # explicit print to the png device
dev.off()

# Console output
print(p_pca)
print("=== PCA Summary ==="); print(pca_summary)
print("=== PCA Loadings (Env Variables) ==="); print(pca_loadings)

# Export tables for downstream analysis
write.csv(pca_scores,   "PCA_species_scores.csv", row.names = TRUE)
write.csv(pca_loadings, "PCA_variable_loadings.csv", row.names = TRUE)

# --------------------------
# 5. Loadings heatmap (PC1‑PC4)
loadings_4pc <- pca_result$rotation[, 1:4]
colnames(loadings_4pc) <- paste0("PC", 1:4, " (", round(pca_summary$importance[2, 1:4]*100, 1), "%)")

pheatmap(loadings_4pc,
         filename = "PCA_loadings_heatmap.png",
         width = 6, 
         height = 8,
         cluster_cols = FALSE,
         color = colorRampPalette(c("#71A2C1","white","#E63946"))(100),
         breaks = seq(-1, 1, length = 101),
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_col = 10,
         fontsize_row = 10)

# --------------------------
# 6. Species‑score heatmap (PC1‑PC4)
colnames(pca_scores_top4) <- paste0("PC", 1:4, " (", round(pca_summary$importance[2, 1:4]*100, 1), "%)")

pheatmap(pca_scores_top4,
         filename = "PCA_species_score_heatmap.png", 
         width = 6, 
         height = 8,
         cluster_cols = FALSE,
         color = colorRampPalette(c("#71A2C1","white","#E63946"))(100),
         breaks = seq(-5, 5, length = 101),
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_col = 10,
         fontsize_row = 10)
