#!/usr/bin/env Rscript

#Ce programme réalise la DEG des genes (GSE19422), et exporte uniquement les données significatif (pvalue.adj < 0.05)

#Input : Data/Pre-traitement/Output/    gene_PCC.csv gene_CONTROL.csv
#Output : Data/DEG/Output/     DEG_results_gene.csv volcano_plot_19422.png

#informations : une version possédant le code pour une heatmap est trouvable dans : Commits on Apr 15, 2025


# Charger les bibliothèques nécessaires
library(limma)
library(ggplot2)

# Importer les données de contrôle et de PCC
control_data <- read.csv("Data/Pre-traitement/Output/gene_CONTROL.csv", row.names=1)
pcc_data <- read.csv("Data/Pre-traitement/Output/gene_PCC.csv", row.names=1)

# Combiner les données de contrôle et de PCC
combined_data <- cbind(control_data, pcc_data)

# Créer un vecteur de groupes (contrôle vs PCC) avec control comme reference 
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

#export ------------------------------------------------------------------------------------------
# Trier les résultats par logFC
deg_results_sorted <- deg_results[order(deg_results$logFC, decreasing=TRUE),]


# Vérifier si le répertoire Output existe
if (!dir.exists("Data/DEG/Output")) {
  # Créer le répertoire s'il n'existe pas
  dir.create("Data/DEG/Output", recursive = TRUE)
}

# Enregistrer les résultats
write.csv(deg_results_sorted, "Data/DEG/Output/DEG_results_gene.csv")


# Enregistrer le vplot 
ggsave("Data/DEG/Output/volcano_plot_gene.png")


print("Fin de l'étape DEG pour les genes.")


