#!/usr/bin/env Rscript

#Ce programme vérifie la qualité dU GSE utilisé (Répetition, NA, ...) exporte les données pour la DEG 

#Input : GEOquery GSE29742 et GPL8227
#Output : Data/Pre-traitement/Output/    miRNA_PCC.csv miRNA_PGL.csv miRNA_CONTROL.csv #miRNA.csv

#informations : une version possédant l'étude des jeux de données et le formatage est disponible dans : Commits on Apr 15, 2025

# Charger les packages
library(GEOquery)
library(Biobase)

#=========== load series and platform data from GEO ===============================

#----- GSE29742 -------------------------------------------------------------------
gset <- getGEO("GSE29742", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL8227", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

miRNA <- exprs(gset)

#----- GPL -------------------------------------------------------------------

GPL8227 <- getGEO("GPL8227", destdir = ".")


#========= export ======================================================
#fusion des matrices 
#==========================================================================================================================
# Indices ajustés pour supprimer les PGL
pgl_indices <- c(2:12)
control_indices <- c(49:54)

# Supprimer les lignes correspondantes au indices 
miRNA_PCC <- miRNA[,-c(pgl_indices,control_indices)]

miRNA_PGL <- miRNA[,pgl_indices]

miRNA_CONTROL <- miRNA[,control_indices]


#-------------------------------------------------------------------------------------------------------------
#export des matrices 

# Vérifier si le répertoire Output existe
if (!dir.exists("Data/Pre-traitement/Output")) {
  # Créer le répertoire s'il n'existe pas
  dir.create("Data/Pre-traitement/Output", recursive = TRUE)
}


write.csv(miRNA_PCC, file = "Data/Pre-traitement/Output/miRNA_PCC.csv", row.names = T)
write.csv(miRNA_PGL, file = "Data/Pre-traitement/Output/miRNA_PGL.csv", row.names = T)
write.csv(miRNA_CONTROL, file = "Data/Pre-traitement/Output/miRNA_CONTROL.csv", row.names = T)
#write.csv(miRNA, file = "Data/Pre-traitement/Output/miRNA.csv", row.names = T)

print("fin du pre-traitement pour miRNA")