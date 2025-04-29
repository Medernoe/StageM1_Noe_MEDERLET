#!/usr/bin/env Rscript


#Ce programme réalise la DEG des miRNA (GSE29742), et exporte uniquement les données significatif (pvalue.adj < 0.05)

#Input : Data/Pre-traitement/Output/    miRNA_CONTROL.csv miRNA_PCC.csv
#Output : Data/DEG/Output/     DEG_results_miRNA.csv volcano_plot_miRNA.png
#         Data/Interaction/Input/     DEG_miRNA.csv

#informations : une version possédant le code pour une heatmap est trouvable dans : Commits on Apr 15, 2025


# Charger les bibliothèques nécessaires
library(limma)
library(ggplot2)

# Importer les données de contrôle et de PCC
control_data <- read.csv("Data/Pre-traitement/Output/miRNA_CONTROL.csv", row.names=1)
pcc_data <- read.csv("Data/Pre-traitement/Output/miRNA_PCC.csv", row.names=1)

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
write.csv(deg_results_sorted, "Data/DEG/Output/DEG_results_miRNA.csv")

# Enregistrer le vplot 
ggsave("Data/DEG/Output/volcano_plot_miRNA.png")


#export dans l'input d'interaction 
# Vérifier si le répertoire Output existe
if (!dir.exists("Data/Interaction/Input")) {
  # Créer le répertoire s'il n'existe pas
  dir.create("Data/Interaction/Input", recursive = TRUE)
}

deg_results_sorted$Symbol = rownames(deg_results_sorted)

# Enregistrer les résultats
write.csv(deg_results_sorted, "Data/Interaction/Input/DEG_miRNA.csv")

print("Fin de l'étape DEG pour les miRNA.")



