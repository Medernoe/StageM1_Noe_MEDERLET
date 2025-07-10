#!/usr/bin/env Rscript

#===============================================================================
# pre_traitement_miRNA.R — qualité GEO + extraction matrices miRNA pour DEG
#
# Usage:
#   Rscript pre_traitement_miRNA.R <pre_out_dir>
#
#   - pre_out_dir : répertoire où écrire les miRNA_*.csv
#===============================================================================

#!/usr/bin/env Rscript

#===============================================================================
# DEG_gene.R — exécute la DEG pour une paire Contrôle (A) / Tumor (B)
#
# Usage:
#   Rscript DEG_gene.R 
#     <control.csv> <tumor.csv> <out_results.csv> <out_plot.png>
#     [logFC_cutoff] [adjP_cutoff]
#
#   - control.csv     : Data/Pre-traitement/Output/miRNA_CONTROL.csv
#   - tumor.csv   : Data/Pre-traitement/Output/miRNA_PCC.csv
#   - out_results.csv : Data/DEG/Output/PCC_Sains/DEG_results_PCC_Sains.csv
#   - out_plot.png    : Data/DEG/Output/PCC_Sains/volcano_plot_PCC_Sains.png
#   - logFC_cutoff    : (optionnel, défaut = 1)
#   - adjP_cutoff     : (optionnel, défaut = 0.05)
#===============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
})

# 1) Lecture des arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript DEG_miRNA.R <control.csv> <tumor.csv> <out_results.csv> <out_plot.png> [logFC_cutoff] [adjP_cutoff]")
}
control_path   <- args[1]
tumor_path <- args[2]
out_results    <- args[3]
out_plot       <- args[4]

# seuils optionnels
logFC_cutoff <- if (length(args) >= 5) as.numeric(args[5]) else 1
adjP_cutoff  <- if (length(args) >= 6) as.numeric(args[6]) else 0.05



# 2) Fonctions LIMMA & plotting
deg_limma <- function(A, B, A_label, B_label, logFC_cutoff, adjP_cutoff) {
  comb <- cbind(A, B)
  grp  <- factor(c(rep(A_label, ncol(A)), rep(B_label, ncol(B))), levels=c(A_label, B_label))
  design <- model.matrix(~ grp)
  fit    <- lmFit(comb, design) |> eBayes()
  topTable(fit, coef=2, adjust.method="BH", number=Inf)
}


plot_volcano <- function(res, logFC_cutoff, adjP_cutoff, title="Volcano Plot") {
  res$Significance <- "Not significant"
  res$Significance[res$adj.P.Val < adjP_cutoff & res$logFC >  logFC_cutoff] <- "Upregulated"
  res$Significance[res$adj.P.Val < adjP_cutoff & res$logFC < -logFC_cutoff] <- "Downregulated"
  res$Significance <- factor(res$Significance, levels=c("Not significant","Upregulated","Downregulated"))
  
  xthr <- c(-logFC_cutoff, logFC_cutoff)
  ythr <- -log10(adjP_cutoff)
  
  ggplot(res, aes(x=logFC, y=-log10(adj.P.Val), color=Significance)) +
    geom_point(alpha=0.5, size=1) +
    geom_vline(xintercept=xthr, linetype="dashed") +
    geom_hline(yintercept=ythr, linetype="dashed") +
    scale_color_manual(values=c("gray","red","blue")) +
    labs(title=title, x="Log2 Fold Change", y="-log10(adj P-value)", color="Significance") +
    theme_minimal()
}


significativity <- function(res, adjP_cutoff, logFC_cutoff) {
  sel <- subset(res, adj.P.Val < adjP_cutoff & abs(logFC) > logFC_cutoff)
  sel[order(sel$logFC, decreasing=TRUE), ]
}


run_analysis <- function(A, B, lblA, lblB, logFC_cutoff, adjP_cutoff) {
  tbl <- deg_limma(A, B, lblA, lblB, logFC_cutoff, adjP_cutoff)
  
  plt <- plot_volcano(tbl, logFC_cutoff, adjP_cutoff,
                      title=paste0("Volcano: ", lblB, " vs ", lblA))
  
  deg <- significativity(tbl, adjP_cutoff, logFC_cutoff)
  list(deg=deg, plot=plt)
}



# 3) Lecture des fichiers et étiquettes
get_label <- function(path) {
  tools::file_path_sans_ext(basename(path))
}

ctrl_df <- read.csv(control_path,   row.names=1)
trt_df  <- read.csv(tumor_path, row.names=1)
lbl_ctrl <- get_label(control_path)
lbl_trt  <- get_label(tumor_path)



# 4) Analyse
res <- run_analysis(ctrl_df, trt_df, lbl_ctrl, lbl_trt, logFC_cutoff, adjP_cutoff)



# 5) Enregistrement

#Ajout de la colonne Symbol
deg_out <- res$deg
deg_out$Symbol <- rownames(deg_out)
write.csv(deg_out, out_results, row.names=FALSE)

# Enregistrement du graphique
ggsave(filename=out_plot, plot=res$plot)


message("DEG miRNA terminé pour : ", lbl_trt, " vs ", lbl_ctrl)

