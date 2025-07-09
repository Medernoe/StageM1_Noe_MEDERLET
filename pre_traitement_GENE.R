#!/usr/bin/env Rscript

#===============================================================================
# pre_traitement_gene.R — qualité GEO + extraction matrices pour DEG et mapping
#
# Usage:
#   Rscript pre_traitement_gene.R <pre_out_dir> <conv_in_dir>
#
#   - pre_out_dir   : répertoire où écrire les CSV de pré-traitement (gene_*.csv)
#   - conv_in_dir   : répertoire où écrire merged_mapping_export.csv
#===============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

# 1) Lecture des arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage:\n  Rscript pre_traitement_gene.R 
  <gene_PCC> 
  <gene_PGL>
  <gene_CONTROL> 
  <gene_PPGL> 
  <merged_map> 
  <pre_out_dir> 
  <conv_in_dir>\n\n"
  )
}

#===============================================================================
# Chargement et annotation des GSE
#===============================================================================
# GSE19422
gset19422 <- getGEO("GSE19422", GSEMatrix = TRUE, AnnotGPL = FALSE)
if (length(gset19422) > 1) {
  idx <- grep("GPL6480", attr(gset19422, "names"))
} else idx <- 1
gset19422 <- gset19422[[idx]]
fvarLabels(gset19422) <- make.names(fvarLabels(gset19422))

# GSE51081
gset51081 <- getGEO("GSE51081", GSEMatrix = TRUE, AnnotGPL = FALSE)
if (length(gset51081) > 1) {
  idx2 <- grep("GPL4133", attr(gset51081, "names"))
} else idx2 <- 1
gset51081 <- gset51081[[idx2]]
fvarLabels(gset51081) <- make.names(fvarLabels(gset51081))

# Plateformes
GPL6480 <- getGEO("GPL6480", destdir = ".")
GPL4133 <- getGEO("GPL4133", destdir = ".")

# Mapping spots → gènes
mapping_6480 <- Table(GPL6480)[, c("SPOT_ID","REFSEQ","GENE","GENE_SYMBOL")]
mapping_4133 <- Table(GPL4133)[, c("SPOT_ID","REFSEQ","GENE","GENE_SYMBOL")]
merged_map   <- merge(mapping_4133, mapping_6480, by = "SPOT_ID",
                      suffixes = c("_4133","_6480"))
merged_map[merged_map==""] <- NA

#===============================================================================
# Extraction des matrices d’intensités
#===============================================================================
DATA1 <- exprs(gset19422)
DATA2 <- exprs(gset51081)

# Remplacement des rownames de DATA2 (ID → SPOT_ID)
new_rows <- mapping_4133$SPOT_ID[match(rownames(DATA2), mapping_4133$ID)]
if (any(is.na(new_rows))) {
  warning("Certains ID de GSE51081 n’ont pas de SPOT_ID correspondant")
} else {
  rownames(DATA2) <- new_rows
}

DATA2 <- unique(DATA2)
DATA1 <- DATA1[order(rownames(DATA1)), ]
DATA2 <- DATA2[order(rownames(DATA2)), ]

#===============================================================================
# Séparation des groupes dans DATA1
#===============================================================================
pgl_idx     <- c(6,7,14,16,18,19,21,25,26,31,34,36,37,41,42,50,59,62,64,75,80,82,84)
control_idx <- c(85,86,87,88,89,90)

DATA1_PCC     <- DATA1[,-c(pgl_idx, control_idx)]
DATA1_PGL     <- DATA1[, pgl_idx]
DATA1_CONTROL <- DATA1[, control_idx]
DATA1_PPGL    <- DATA1[,-control_idx]

#===============================================================================
# Écriture des CSV
#===============================================================================

output_pre_out   <- args[6]  # répertoire de sortie pour les CSV
output_conv_in   <- args[7]  # répertoire de sortie pour merged_mapping_export.csv
# Création des dossiers s’ils n’existent pas
if (!dir.exists(output_pre_out)) dir.create(output_pre_out, recursive = TRUE)
if (!dir.exists(output_conv_in)) dir.create(output_conv_in, recursive = TRUE)


output_gene_PCC   <- args[1]  
output_gene_PGL   <- args[2] 
output_gene_CONTROL <- args[3]
output_gene_PPGL  <- args[4] 

write.csv(DATA1_PCC,     file = output_gene_PCC ,     row.names = TRUE)
write.csv(DATA1_PGL,     file = output_gene_PGL,     row.names = TRUE)
write.csv(DATA1_CONTROL, file = output_gene_CONTROL, row.names = TRUE)
write.csv(DATA1_PPGL,    file = output_gene_PPGL,    row.names = TRUE)

# Mapping export

output_merged_map <- args[5]

merged_map <- merged_map[!(is.na(merged_map$GENE_4133) & is.na(merged_map$GENE_6480)), ]

# Écriture du fichier de mapping
write.csv(merged_map, file = output_merged_map,
          row.names = TRUE)

message("Fin du pré-traitement genes.")
