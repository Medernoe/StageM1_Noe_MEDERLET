# Charger les bibliothèques nécessaires
library(limma)
library(ggplot2)
install.packages(pheatmap)
library(pheatmap)

# Importer les données de contrôle et de PCC
control_data <- read.csv("DATA19422_CONTROL.csv", row.names=1)
pcc_data <- read.csv("DATA19422_PCC.csv", row.names=1)

# Combiner les données de contrôle et de PCC
combined_data <- cbind(control_data, pcc_data)

# Créer un vecteur de groupes (contrôle vs PCC)
group_labels <- factor(c(rep("Control", ncol(control_data)), rep("PCC", ncol(pcc_data))))

# Créer une matrice de design
design <- model.matrix(~group_labels)

# Appliquer limma pour l'analyse de DEG
fit <- lmFit(combined_data, design)
fit <- eBayes(fit)

# Obtenir les résultats avec les p-values ajustées
results <- topTable(fit, coef=2, adjust="BH", number=Inf)

# Filtrer les gènes avec une p-value ajustée < 0.05
deg_results <- results[results$adj.P.Val < 0.05,]




#figure --------------------------------------------------------------------------------------------
# Visualisation : volcano plot pour DEG
volcano_data <- data.frame(logFC = deg_results$logFC, PValue = deg_results$P.Value)
ggplot(volcano_data, aes(x = logFC, y = -log10(PValue))) + 
  geom_point(aes(color = PValue < 0.05)) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log Fold Change", y = "-Log10 P-value") +
  scale_color_manual(values = c("gray", "red"))

# Visualisation : Heatmap des gènes différentiellement exprimés
# Filtrer les miARN avec un ajustement p-value FDR < 0.05
top_deg <- deg_results[deg_results$adj.P.Val < 0.05, ]
top_deg_genes <- rownames(top_deg)
heatmap_data <- data[top_deg_genes, ]

# Créer la heatmap
pheatmap(heatmap_data, scale = "row", clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")



#export ------------------------------------------------------------------------------------------
# Trier les résultats par logFC
deg_results_sorted <- deg_results[order(deg_results$logFC, decreasing=TRUE),]

# Afficher les résultats
write.csv(deg_results_sorted, "DEG_results.csv")
