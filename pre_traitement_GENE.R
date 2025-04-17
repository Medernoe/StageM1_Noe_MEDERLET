#!/usr/bin/env Rscript

#Ce programme vérifie la qualité des GSE utilisés (Répetition, NA, ...) exporte les données pour la DEG et la table de mapping pour la conversion_spotID_geneID

#Input : GEOquery GSE19422/51081 et GPL6480
#Output : Data/Pre-traitement/Output/    gene_PCC.csv gene_PGL.csv gene_CONTROL.csv      #DATA51081.csv
#         Data/Conversion/Input/         merged_mapping_export.csv

#informations : une version possédant l'étude des jeux de données et le formatage est disponible dans : Commits on Apr 15, 2025

# Charger les packages
library(GEOquery)
library(Biobase)


# load series and platform data from GEO
#==========================================================================================================================
# GSE19422 ----------------------------------------------------------------------------------------------------------------------------------------------------------
gset19422 <- getGEO("GSE19422", GSEMatrix =TRUE, AnnotGPL=FALSE)
#ajout de la plateform GPL 
if (length(gset19422) > 1) idx <- grep("GPL6480", attr(gset19422, "names")) else idx <- 1
gset19422 <- gset19422[[idx]]

# make proper column names to match toptable 
fvarLabels(gset19422) <- make.names(fvarLabels(gset19422))


# GSE51081 -----------------------------------------------------------------------------------------------------------------------------------------------------------
 
gset51081 <- getGEO("GSE51081", GSEMatrix =TRUE, AnnotGPL=FALSE)
#ajout de la plateform GPL 
if (length(gset51081) > 1) idx <- grep("GPL4133", attr(gset51081, "names")) else idx <- 1
gset51081 <- gset51081[[idx]]

# make proper column names to match toptable 
fvarLabels(gset51081) <- make.names(fvarLabels(gset51081))

#==========================================================================================================================
# Télécharger l'annotation de la plateforme GPL6480-gset19422---------------------------------------------------------------------------------------------
GPL6480 <- getGEO("GPL6480", destdir = ".")

#créations de la table de mapping 
mapping_GPL6480 <- Table(GPL6480)[, c("ID", "SPOT_ID","REFSEQ","GENE","GENE_SYMBOL")]

# Télécharger l'annotation de la plateforme GPL4133-gset51081---------------------------------------------------------------------------------------------
GPL4133 <- getGEO("GPL4133", destdir = ".")

#créations de la table de mapping 
mapping_GPL4133 <- Table(GPL4133)[, c("ID", "SPOT_ID" , "REFSEQ","GENE","GENE_SYMBOL")]


# Fusionner les deux table de mapping sur la colonne SPOT_ID
merged_mapping <- merge(mapping_GPL4133, mapping_GPL6480, by = "SPOT_ID", suffixes = c("_4133", "_6480"))
merged_mapping[merged_mapping == ""] <- NA


#==========================================================================================================================
#matrice des microarray

#matrice des résultats, ligne Gènes, colonne GSE, cellule intensité lumineuse 
DATA1 = exprs(gset19422)
DATA2 = exprs(gset51081)


#-------------------------------------------------------------------------------------------------------------
#transformation des noms des lignes de "ID" vers "SPOT_ID" pour gse 51081 

# Créer un vecteur contenant les nouveaux noms de ligne (SPOT_ID)
new_rownames <- mapping_GPL4133$SPOT_ID[match(rownames(exprs(gset51081)), mapping_GPL4133$ID)]

# Appliquer le remplacement 
if (sum(is.na(new_rownames)) == 0) {
  rownames(DATA2) <- new_rownames
} else {
  warning("dimension incompatible.")
}


#-------------------------------------------------------------------------------------------------------------
#suppresion des lignes identiques dans datas 2 
DATA2 <- unique(DATA2)

#trie dans le meme ordre 
DATA2 <- DATA2[order(rownames(DATA2)), ]
DATA1 <- DATA1[order(rownames(DATA1)), ]



#fusion des matrices 
#==========================================================================================================================
#suppresion des PGL du GSE19422 

# Indices ajustés pour supprimer les PGL
pgl_indices <- c(6, 7, 14, 16, 18, 19, 21, 25, 26, 31, 34, 36, 37, 41, 42, 50, 59, 62, 64, 75, 80, 82, 84)
control_indices <- c(85, 86, 87, 88, 89, 90)


# Supprimer les lignes correspondantes
DATA1_PCC <- DATA1[,-c(pgl_indices,control_indices)]

DATA1_PGL <- DATA1[,pgl_indices]

DATA1_CONTROL <- DATA1[,control_indices]

#-------------------------------------------------------------------------------------------------------------
#fusion des 2 matrices 


data_combined <- cbind(DATA1_PCC, DATA2)

print(dim(data_combined))


#-------------------------------------------------------------------------------------------------------------
#export des matrices 

# Vérifier si le répertoire Output existe
if (!dir.exists("Data/Pre-traitement/Output")) {
  # Créer le répertoire s'il n'existe pas
  dir.create("Data/Pre-traitement/Output", recursive = TRUE)
}


write.csv(DATA1_PCC, file = "Data/Pre-traitement/Output/gene_PCC.csv", row.names = T)
write.csv(DATA1_PGL, file = "Data/Pre-traitement/Output/gene_PGL.csv", row.names = T)
write.csv(DATA1_CONTROL, file = "Data/Pre-traitement/Output/gene_CONTROL.csv", row.names = T)
#write.csv(DATA2, file = "Data/Pre-traitement/Output/DATA51081.csv", row.names = T)

#-------------------------------------------------------------------------------------------------------------
#export de merge_mapping 


if (!dir.exists("Data/Conversion/Input")) {
  # Créer le répertoire s'il n'existe pas
  dir.create("Data/Conversion/Input", recursive = TRUE)
}


merged_mapping_export <- merged_mapping[!(is.na(merged_mapping$GENE_4133) & is.na(merged_mapping$GENE_6480)), ]

# Écrire le fichier CSV
write.csv(merged_mapping_export, file = "Data/Conversion/Input/merged_mapping_export.csv", row.names = TRUE)

print("fin du pre-traitement pour les genes")

