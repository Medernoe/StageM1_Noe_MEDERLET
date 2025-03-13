# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
library(GEOquery)
library(limma)
library(umap)
library(Biobase)
library(dplyr)


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

# Vérifier la structure de la table des formats de nom des lignes 
colnames(Table(GPL6480))
Table(GPL6480)[1:5, ]

#créations de la table de mapping 
mapping_GPL6480 <- Table(GPL6480)[, c("ID", "SPOT_ID","REFSEQ","GENE","GENE_SYMBOL")]
length(mapping_GPL6480$ID) 


# Télécharger l'annotation de la plateforme GPL4133-gset51081---------------------------------------------------------------------------------------------
GPL4133 <- getGEO("GPL4133", destdir = ".")

# Vérifier la structure de la table d'annotation
colnames(Table(GPL4133))
Table(GPL4133)[1:5, ]

#créations de la table de mapping 
mapping_GPL4133 <- Table(GPL4133)[, c("ID", "SPOT_ID" , "REFSEQ","GENE","GENE_SYMBOL")]
length(mapping_GPL4133$ID) 

#les tables d'annotation sont bien différentes, mais cela n'importe peut car elles permettent juste de connvertir l'annotation des GSE en mettant une annotation commune 
#l'annotation commune sélectionner est "SPOT_ID" 



#lecture des jeux de données 
#==========================================================================================================================
#-------------------------------------------------------------------------------------------------------------
# Affichage des metadata d'un GSM
pData(gset19422)

#matrice des résultats, ligne Gènes, colonne GSE, cellule intensité lumineuse 
DATA1 = exprs(gset19422)
DATA2 = exprs(gset51081)


# Vérification des dimensions
dim(exprs(gset19422))  # 41000 (Nombre de gène) x 90 (nombre de GSE)
dim(exprs(gset51081))  # 43376 x 9 (nombre de GSE)
#il y a 2376 lignes en plus dans le GSM51081 (possible lignes dupliqué, lignes vide, ... )

# Vérifier les noms des échantillons (GSE)
colnames(exprs(gset19422))
colnames(exprs(gset51081))

# Vérifier les noms des gènes
rownames(featureData(gset19422)) #Probe name version (= "Spot_ID" )
rownames(featureData(gset51081)) #Feature name version (="ID")
#attention "Id" n'est pas la meme chose dans les 2 PGL CF: structuration des données sur le liens du GSE 

#vérification de NA dans les data 
sum(is.na(DATA1))
sum(is.na(DATA2))
#aucun NA

#vérification que toute les lignes aie un nom 
sum(is.na(rownames(DATA1)))
sum(is.na(rownames(DATA2)))

sum(rownames(DATA1) == "")
sum(rownames(DATA2) == "")

all(!is.na(rownames(data)) & rownames(data) != "")
#aucun NA


#vérification des information  
#==========================================================================================================================
#vérifie le nombre de sport commun entre les mapping 
common_spots <- intersect(mapping_GPL4133$SPOT_ID, mapping_GPL6480$SPOT_ID)
length(common_spots)  # Vérifie combien de valeurs sont communes

#-------------------------------------------------------------------------------------------------------------
#vérification que les Spot_id des GSE sont lié au meme refseq / gene / symbole

# Fusionner les deux matrices sur la colonne SPOT_ID
merged_mapping <- merge(mapping_GPL4133, mapping_GPL6480, by = "SPOT_ID", suffixes = c("_4133", "_6480"))
merged_mapping[merged_mapping == ""] <- NA



#-------------------------------------------------------------------------------------------------------------
#supprimer les lignes contenant des NA 


#supprimer les lignes contenant des NA dans un des partie gene 
merged_mapping_gene <- merged_mapping[!is.na(merged_mapping$GENE_4133) & !is.na(merged_mapping$GENE_6480), ]
length(merged_mapping_gene)

#supprimer les lignes contenant des NA dans un des partie refseq 
merged_mapping_refseq <- merged_mapping[!is.na(merged_mapping$REFSEQ_4133) & !is.na(merged_mapping$REFSEQ_6480), ]
length(merged_mapping_refseq)

#supprimer les lignes contenant des NA dans un des partie refseq 
merged_mapping_symbol <- merged_mapping[!is.na(merged_mapping$GENE_SYMBOL_4133) & !is.na(merged_mapping$GENE_SYMBOL_6480), ]
length(merged_mapping_symbol)


#beaucoup plus de lignes avec des NA pour refseq que pour gene 
#avant retrait 45220, refseq 28940 , gene 31757 , symbole 31757
#peut etre sélectionner pour GSE19422 les id de GSE51081 quand ils sont vident 


#-------------------------------------------------------------------------------------------------------------
# Vérifier les différences entre les refseq/gene des deux matrices de mapping 
diff_rows_refseq <- na.omit( merged_mapping[merged_mapping$REFSEQ_4133 != merged_mapping$REFSEQ_6480, ] ) 
diff_rows_gene <- na.omit(merged_mapping[merged_mapping$GENE_4133 != merged_mapping$GENE_6480, ])
diff_rows_symbol <- na.omit(merged_mapping[merged_mapping$GENE_SYMBOL_4133 != merged_mapping$GENE_SYMBOL_6480, ])


# Afficher les lignes qui ne sont pas similaires
print(length(rownames(diff_rows_refseq))) #901
print(length(rownames(diff_rows_gene))) #219
print(length(rownames(diff_rows_symbol))) #1749

#possibilité à voir plus tard de prendre les information du GPL6480 

#-------------------------------------------------------------------------------------------------------------
#vérification de la présence de gene ou de refseq dulique 
sum(duplicated(merged_mapping[,1]))# Vérifier s'il y a des SPOT_ID dupliqués
sum(duplicated(merged_mapping[,2]))
sum(duplicated(merged_mapping[,3]))
sum(duplicated(merged_mapping[,4]))
sum(duplicated(merged_mapping[,5]))
sum(duplicated(merged_mapping[,6]))
sum(duplicated(merged_mapping[,7]))
table(table(merged_mapping[,1]))


#gros problème, beaucoup d'identifiant (tout type) sont dupliqué, que faire ceux si ? 
#si on utilise Spot_id dirige vers plusieurs gène 

#il est possible qu'ils s'agissent de spot_id non liés aux gènes (ex: (-)3xSLv1)


#vérification des duplication directement dans les tables gene/merged_mapping_gene 
sum(duplicated(merged_mapping_gene[,4]))
table(table(merged_mapping_gene[,4]))


sum(duplicated(merged_mapping_refseq[,3]))
table(table(merged_mapping_refseq[,3]))

merged_mapping_na <- na.omit(merged_mapping)
table(table(merged_mapping_na[,4]))

#==========================================================================================================================
#-------------------------------------------------------------------------------------------------------------
#vérification qu'il n'y ait pas de différences entre les noms de spot de GSE_19422 et la matrice 

# Extraire les noms de la colonne de référence
ref_names <- mapping_GPL6480[, "SPOT_ID"] 

# Extraire les noms des lignes de la matrice d'intensité
intensity_names <- rownames(DATA1)

# Vérifier si tous les noms de la matrice d'intensité sont dans la référence
all(intensity_names %in% ref_names)  # TRUE si tous sont présents, FALSE sinon


# Obtenir les noms absents dans la référence
noms_absents <- setdiff(intensity_names, ref_names)
print(noms_absents)  

#tout les noms du gse_19422 sont également présent dans la référence 


#-------------------------------------------------------------------------------------------------------------
#transformation des noms des lignes de "ID" vers "SPOT_ID" pour gse 51081 

# Créer un vecteur contenant les nouveaux noms de ligne (SPOT_ID)
new_rownames <- mapping_GPL4133$SPOT_ID[match(rownames(exprs(gset51081)), mapping_GPL4133$ID)]

# Vérifier s'il y a des valeurs NA (signifiant des ID sans correspondance)
sum(is.na(new_rownames))

#dimension avant échange 
length(new_rownames)
length(rownames(exprs(gset51081)))
#meme longueur 

# Appliquer le remplacement 
if (sum(is.na(new_rownames)) == 0) {
  rownames(DATA2) <- new_rownames
} else {
  warning("dimension incompatible.")
}

#-------------------------------------------------------------------------------------------------------------
#vérification qu'il n'y ait pas de différences entre les noms de spot de DATA2 (GSE_51081) et la matrice mapping_GPL6480

# Extraire les noms de la colonne de référence
ref_names <- mapping_GPL6480[, "SPOT_ID"] 

# Extraire les noms des lignes de la matrice d'intensité
intensity_names <- rownames(DATA2)

# Vérifier si tous les noms de la matrice d'intensité sont dans la référence
all(intensity_names %in% ref_names)  # TRUE si tous sont présents, FALSE sinon

# Obtenir les noms absents dans la référence
noms_absents <- setdiff(intensity_names, ref_names)
print(noms_absents)  


#les spot_id nouvellement implémenté à partir de la matrice mapping_GPL4133 sont identique au spot_id présent dans la matrice mapping_GPL6480 


#-------------------------------------------------------------------------------------------------------------
#vérification des duplications des noms des SPOt_ID dans DATA2/DATA1 

#vérification des duplications et de leur forme 

#data2 
dup_spot_ids <- rownames(DATA2)[duplicated(rownames(DATA2))]

sum(duplicated(rownames(DATA2))) # Vérifier s'il y a des SPOT_ID dupliqués
table(table(rownames(DATA2)))
table(rownames(DATA2))[table(rownames(DATA2)) > 1]
#présence uniquement de ligne répétée 10 fois 


# Extraire les lignes dupliquées pour faire une matrice
dup_spot_ids <- names(which(table(rownames(DATA2)) > 1))
duplicated_rows <- DATA2[(rownames(DATA2)) %in% dup_spot_ids, ]
duplicated_rows <- duplicated_rows[order(rownames(duplicated_rows)), ]

#dans la matrice on voit que les lignes sont répétée 10 fois, les lignes avec le meme spot_id sont identique 


#data1 

dup_spot_ids <- rownames(DATA1)[duplicated(rownames(DATA1))]

sum(duplicated(rownames(DATA1))) # Vérifier s'il y a des SPOT_ID dupliqués
table(table(rownames(DATA1)))
table(rownames(DATA1))[table(rownames(DATA1)) > 1]



#-------------------------------------------------------------------------------------------------------------
#suppresion des lignes identiques dans datas 2 
DATA2 <- unique(DATA2)

length(rownames(DATA1))
length(rownames(DATA2))

#il y a le meme nombre de ligne dans DATA1 et DATA 2 
#les noms de lignes de DATA1 et DATA2 sont-ils identique? 

DATA2 <- DATA2[order(rownames(DATA2)), ]
DATA1 <- DATA1[order(rownames(DATA1)), ]

#les noms de lignes de DATA1 et DATA2 sont-ils identique?  

#-------------------------------------------------------------------------------------------------------------
#vérification que les lignes soient bien identique 

data1_names = rownames(DATA1)
data2_names = rownames(DATA2)

all(data1_names %in% data2_names)

#les lignes sont identiques 


DATA1[DATA1 == ""]
DATA2[DATA2 == ""]
any(is.na(DATA1))
any(is.na(DATA2))

#aucune valeur manquante dans nos DATA1 et DATA2 



#fusion des matrices 
#==========================================================================================================================
#suppresion des PGL du GSE19422 

# Indices ajustés pour supprimer les PGL
pgl_indices <- c(6, 7, 14, 16, 18, 19, 21, 25, 26, 31, 34, 36, 37, 41, 42, 50, 59, 62, 64, 75, 80, 82, 84)
control_indices <- c(85, 86, 87, 88, 89, 90)

dim(DATA1) #90

# Supprimer les lignes correspondantes
DATA1_PCC <- DATA1[,-c(pgl_indices,control_indices)]

DATA1_PGL <- DATA1[,pgl_indices]

DATA1_CONTROL <- DATA1[,control_indices]


# Vérifier le résultat
dim(DATA1_PCC) #61
dim(DATA1_PGL) #23
dim(DATA1_CONTROL) #6


#-------------------------------------------------------------------------------------------------------------
#fusion des 2 matrices 


data_combined <- cbind(DATA1_PCC, DATA2)

print(dim(data_combined))


#-------------------------------------------------------------------------------------------------------------
#export des matrices 

write.csv(DATA1_PCC, file = "DATA19422_PCC.csv", row.names = T)
write.csv(DATA1_PGL, file = "DATA19422_PGL.csv", row.names = T)
write.csv(DATA1_CONTROL, file = "DATA19422_CONTROL.csv", row.names = T)
write.csv(DATA2, file = "DATA51081.csv", row.names = T)









