# Charger les bibliothèques nécessaires
library(limma)
library(ggplot2)
library(pheatmap)
library(ggplot2)

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

#volcano plot 


# Ajouter une colonne 'Significance' 
results$Significance <- "Not significant"

# Mise à jour des valeurs selon les conditions
results$Significance[results$P.Value < 0.05 & results$logFC > 1] <- "Upregulated"
results$Significance[results$P.Value < 0.05 & results$logFC < -1] <- "Downregulated"

# Convertir en facteur 
results$Significance <- factor(results$Significance, levels = c("Not significant", "Upregulated", "Downregulated"))

# plot
ggplot(results, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change (logFC)",
       y = "-log10(P-value)",
       color = "Expression") +
  theme_minimal() +
  theme(legend.position = "right")




#heatmap - brouilon 

# Supposons que "results" est le tableau des résultats DEG (avec logFC et P.Value)
# On garde les gènes avec |logFC| > 1 et p-value < 0.05
deg_genes <- rownames(results[abs(results$logFC) > 1 & results$P.Value < 0.05, ])

# Filtrer la matrice DATA pour ne garder que ces gènes
mat_deg <- combined_data[deg_genes, ]

# Générer la heatmap
pheatmap(as.matrix(mat_deg), 
         scale = "row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap des gènes différentiellement exprimés")



#export ------------------------------------------------------------------------------------------
# Trier les résultats par logFC
deg_results_sorted <- deg_results[order(deg_results$logFC, decreasing=TRUE),]

# Afficher les résultats
write.csv(deg_results_sorted, "DEG_results_19422.csv")
