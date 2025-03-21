# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

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
Table(GPL6480)[1:5, ]


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



# Supprimer les lignes correspondantes
miRNA_PCC <- miRNA[,-c(pgl_indices,control_indices)]

miRNA_PGL <- miRNA[,pgl_indices]

miRNA_CONTROL <- miRNA[,control_indices]


# Vérifier le résultat
dim(miRNA_PCC) #37
dim(miRNA_PGL) #11
dim(miRNA_CONTROL) #6



#-------------------------------------------------------------------------------------------------------------
#export des matrices 

write.csv(miRNA_PCC, file = "miRNA_PCC.csv", row.names = T)
write.csv(miRNA_PGL, file = "miRNA_PGL.csv", row.names = T)
write.csv(miRNA_CONTROL, file = "miRNA_CONTROL.csv", row.names = T)
write.csv(miRNA, file = "miRNA.csv", row.names = T)

