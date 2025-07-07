#!/usr/bin/env Rscript

#Ce programme réalise la DEG des genes (GSE19422), et exporte uniquement les données significatif (pvalue.adj < 0.05)

#Input : Data/Pre-traitement/Output/    gene_PCC.csv gene_CONTROL.csv
#Output : Data/DEG/Output/     DEG_results_gene.csv volcano_plot_19422.png

#informations : une version possédant le code pour une heatmap est trouvable dans : Commits on Apr 15, 2025


# Charger les packages
library(limma)
library(ggplot2)

#--------------------------------------------------------------------------
# 1) Calcul des résultats limma
deg_limma <- function(A, B, A_label, B_label, logFC_cutoff = 1, adjP_cutoff = 0.05) {
  combined_data <- cbind(A, B)
  group_labels  <- factor(
    c(rep(A_label, ncol(A)), rep(B_label, ncol(B))),
    levels = c(A_label, B_label)
  )
  design <- model.matrix(~ group_labels)
  fit    <- lmFit(combined_data, design)
  fit    <- eBayes(fit)
  topTable(fit, coef = 2, adjust.method = "BH", number = Inf)
}

#--------------------------------------------------------------------------
# 2) Volcano plot
plot_volcano <- function(results, logFC_cutoff = 1, adjP_cutoff = 0.05, title = "Volcano Plot") {
  # Préparer les labels
  results$Significance <- "Not significant"
  results$Significance[results$adj.P.Val < adjP_cutoff & results$logFC >  logFC_cutoff]  <- "Upregulated"
  results$Significance[results$adj.P.Val < adjP_cutoff & results$logFC < -logFC_cutoff]  <- "Downregulated"
  results$Significance <- factor(
    results$Significance,
    levels = c("Not significant", "Upregulated", "Downregulated")
  )
  
  # Seuils
  x_thresh <- c(-logFC_cutoff, logFC_cutoff)
  y_thresh <- -log10(adjP_cutoff)
  
  # Construction du graphique
  ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_vline(xintercept = x_thresh, linetype = "dashed") +
    geom_hline(yintercept = y_thresh, linetype = "dashed") +
    scale_color_manual(values = c(
      "Not significant" = "gray",
      "Upregulated"     = "red",
      "Downregulated"   = "blue"
    )) +
    labs(
      title = title,
      x     = "Log2 Fold Change",
      y     = "-log10(adj P-value)",
      color = "Significance"
    ) +
    theme_minimal()
}

#--------------------------------------------------------------------------
# 3) Filtre des DEG significatifs triés
significativity <- function(results, adjP_cutoff = 0.05, logFC_cutoff = 1) {
  deg_results <- subset(
    results,
    adj.P.Val < adjP_cutoff & abs(logFC) > logFC_cutoff
  )
  deg_results[order(deg_results$logFC, decreasing = TRUE), ]
}

#--------------------------------------------------------------------------
# 4) Pipeline 
run_analysis <- function(A,
                         B, 
                         A_label, 
                         B_label, 
                         logFC_cutoff = 1, 
                         adjP_cutoff = 0.05) {
  
  res    <- deg_limma(A, B, A_label, B_label, logFC_cutoff, adjP_cutoff)
  
  volcano <- plot_volcano(
    res,
    logFC_cutoff,
    adjP_cutoff,
    title = paste0("Volcano: ", B_label, " vs ", A_label)
  )
  
  deg_results <- significativity(res, adjP_cutoff, logFC_cutoff)
  
  list(deg  = deg_results, plot = volcano)
}



#--------------------------------------------------------------------------
# Import et fusion
Control <- read.csv("Data/Pre-traitement/Output/miRNA_CONTROL.csv", row.names = 1)
PCC     <- read.csv("Data/Pre-traitement/Output/miRNA_PCC.csv",     row.names = 1)
PGL     <- read.csv("Data/Pre-traitement/Output/miRNA_PGL.csv",     row.names = 1)
PPGL    <- cbind(PCC, PGL)


# Calculs
results_list <- list(
  PCC_vs_Control  = run_analysis(Control, PCC,  "Control", "PCC"),
  PGL_vs_Control  = run_analysis(Control, PGL,  "Control", "PGL"),
  PPGL_vs_Control = run_analysis(Control, PPGL, "Control", "PPGL"),
  PCC_vs_PGL      = run_analysis(PGL,    PCC,  "PGL",     "PCC")
)

# Automatisation de l’export
out_dir <- "Data/DEG/Output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (name in names(results_list)) {
  # Extraction
  deg_df  <- results_list[[name]]$deg
  volcano <- results_list[[name]]$plot
  
  # Fichiers de sortie
  csv_file  <- file.path(out_dir, paste0("DEM_results_", name, ".csv"))
  plot_file <- file.path(out_dir, paste0("miRNA_volcano_plot_",  name, ".png"))
  
  # Sauvegarde
  write.csv(deg_df, csv_file, row.names = TRUE)
  ggsave(plot_file, volcano)
}

message("Fin de l'étape DEG pour les genes. Les fichiers sont dans ", out_dir, ".")

