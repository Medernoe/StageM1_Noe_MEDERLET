#!/usr/bin/env Rscript

#===============================================================================
# pre_traitement_miRNA.R — qualité GEO + extraction matrices miRNA pour DEG
#
# Usage:
#   Rscript pre_traitement_miRNA.R <pre_out_dir>
#
# Arguments :
#   pre_out_dir    : répertoire où écrire miRNA_*.csv
#===============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

# 1) Lecture des arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage:\n  Rscript pre_traitement_miRNA.R 
  <miRNA_PCC> 
  <miRNA_PGL>
  <miRNA_CONTROL> 
  <miRNA_PPGL>
  <pre_out_dir>")
}

#===============================================================================
# 2) Chargement des données GEO et plateforme
#===============================================================================
gset <- getGEO("GSE29742", GSEMatrix = TRUE, AnnotGPL = FALSE)
if (length(gset) > 1) {
  idx <- grep("GPL8227", attr(gset, "names"))
} else idx <- 1
gset <- gset[[idx]]
fvarLabels(gset) <- make.names(fvarLabels(gset))
miRNA <- exprs(gset)

#===============================================================================
# 3) Séparation des groupes
#===============================================================================
pgl_idx     <- c(2:12)
control_idx <- c(49:54)

miRNA_PCC     <- miRNA[,-c(pgl_idx, control_idx)]
miRNA_PGL     <- miRNA[, pgl_idx]
miRNA_CONTROL <- miRNA[, control_idx]
miRNA_PPGL    <- miRNA[,-control_idx]

#===============================================================================
# 4) Écriture des CSV dans pre_out_dir
#===============================================================================




# Répertoire de sortie pour les CSV
output_pre_out     <- args[5]  # ex. "Data/Pre-traitement/Output"
# Création des dossiers s’ils n’existent pas
if (!dir.exists(output_pre_out))   dir.create(output_pre_out,   recursive = TRUE)



output_miRNA_PCC   <- args[1] 
output_miRNA_PGL   <- args[2]  
output_miRNA_CONTROL <- args[3]  
output_miRNA_PPGL  <- args[4]  

write.csv(miRNA_PCC,     file = output_miRNA_PCC,     row.names = TRUE)
write.csv(miRNA_PGL,     file = output_miRNA_PGL,     row.names = TRUE)
write.csv(miRNA_CONTROL, file = output_miRNA_CONTROL, row.names = TRUE)
write.csv(miRNA_PPGL,    file = output_miRNA_PPGL,    row.names = TRUE)

message("Fin du pré-traitement miRNA")
