#!/usr/bin/env Rscript

#Ce programme vérifie la qualité dU GSE utilisé (Répetition, NA, ...) exporte les données pour la DEG 

#Input : GEOquery GSE29742 et GPL8227
#Output : Pre-traitement/Output/    miRNA_PCC.csv miRNA_PGL.csv miRNA_CONTROL.csv miRNA.csv


# Définir le répertoire de travail
setwd("/home/medernoe/Ubuntu_user/Stage_M1_ubuntu/StageM1_git/Data")


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

# Vérifier la structure de la table des formats de nom des lignes 
colnames(Table(GPL8227))
Table(GPL8227)[1:5, ]


#========= étude des données ======================================================

sum(is.na(miRNA)) # 0 Na 
sum(is.na(rownames(miRNA))) # 0 Na 
sum(rownames(miRNA) == "") # 0 lignes vides  

table(table(rownames(miRNA))) #aucune ligne avec le meme nom 


all(rownames(miRNA) %in% Table(GPL8227)[,1])  # TRUE si tous sont présents, FALSE sinon
#tout les noms de miRNA sont dans la table 


#aucune données ne doit etre filter, la matrice semble parfaite comme input 

#========= export ======================================================
#fusion des matrices 
#==========================================================================================================================
#suppresion des PGL du GSE19422 

# Indices ajustés pour supprimer les PGL
pgl_indices <- c(2:12)
control_indices <- c(49:54)



# Supprimer les lignes correspondantes au indices 
miRNA_PCC <- miRNA[,-c(pgl_indices,control_indices)]

miRNA_PGL <- miRNA[,pgl_indices]

miRNA_CONTROL <- miRNA[,control_indices]


# Vérifier le résultat
dim(miRNA_PCC) #37
dim(miRNA_PGL) #11
dim(miRNA_CONTROL) #6



#-------------------------------------------------------------------------------------------------------------
#export des matrices 

# Vérifier si le répertoire Output existe
if (!dir.exists("Pre-traitement/Output")) {
  # Créer le répertoire s'il n'existe pas
  dir.create("Pre-traitement/Output", recursive = TRUE)
}


write.csv(miRNA_PCC, file = "Pre-traitement/Output/miRNA_PCC.csv", row.names = T)
write.csv(miRNA_PGL, file = "Pre-traitement/Output/miRNA_PGL.csv", row.names = T)
write.csv(miRNA_CONTROL, file = "Pre-traitement/Output/miRNA_CONTROL.csv", row.names = T)
write.csv(miRNA, file = "Pre-traitement/Output/miRNA.csv", row.names = T)

print("fin du pre-traitement")