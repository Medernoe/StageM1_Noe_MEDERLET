#!/usr/bin/env Rscript

# Définir le répertoire d'input
setwd("/home/medernoe/Ubuntu_user/Stage_M1_ubuntu/StageM1_git/Data")


# Charger les bibliothèques nécessaires
library(limma)
library(ggplot2)
library(pheatmap)

# Importer les données de contrôle et de PCC
control_data <- read.csv("Pre-traitement/Output/miRNA_CONTROL.csv", row.names=1)
pcc_data <- read.csv("Pre-traitement/Output/miRNA_PCC.csv", row.names=1)

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
X11()
#volcano plot 


# Ajouter une colonne 'Significance' 
results$Significance <- "Not significant"

# Mise à jour des valeurs selon les conditions
results$Significance[results$P.Value < 0.05 & results$logFC > 1] <- "Upregulated"
results$Significance[results$P.Value < 0.05 & results$logFC < -1] <- "Downregulated"

# Convertir en facteur 
results$Significance <- factor(results$Significance, levels = c("Not significant", "Upregulated", "Downregulated"))

# plot
plot <- ggplot(results, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change (logFC)",
       y = "-log10(P-value)",
       color = "Expression") +
  theme_minimal() +
  theme(legend.position = "right")

print(plot)



#heatmap - brouillon 
dev.new()

# Assurer que combined_data est bien une matrice numérique
data_matrix <- as.matrix(combined_data)
mode(data_matrix) <- "numeric"

# Vérifier que les gènes DEG existent dans la matrice
common_genes <- intersect(rownames(deg_results), rownames(data_matrix))

# Extraire uniquement les gènes DEG pour la heatmap
heatmap_data <- data_matrix[common_genes, ]

# Vérifier que les données ne contiennent pas de valeurs NA
heatmap_data <- na.omit(heatmap_data)

# plot
plot2 <- pheatmap(heatmap_data, scale="row", main="Heatmap des DEG", color=colorRampPalette(c("blue", "white", "red"))(50))

print(plot2)
#export ------------------------------------------------------------------------------------------
# Trier les résultats par logFC
deg_results_sorted <- deg_results[order(deg_results$logFC, decreasing=TRUE),]

# Vérifier si le répertoire Output existe
if (!dir.exists("DEG/Output")) {
  # Créer le répertoire s'il n'existe pas
  dir.create("DEG/Output", recursive = TRUE)
}



# Enregistrer les résultats
write.csv(deg_results_sorted, "DEG/Output/DEG_results_miRNA.csv")

# Enregistrer le vplot 
ggsave("DEG/Output/volcano_plot_miRNA.png")




