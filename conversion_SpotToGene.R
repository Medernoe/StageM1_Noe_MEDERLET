#!/usr/bin/env Rscript

# Définir le répertoire de travail
setwd("/home/medernoe/Ubuntu_user/Stage_M1_ubuntu/StageM1_git/Data")

updated_mapping <- read.csv("Conversion/Output/update_gene_mapping.csv")
DEG_19422 <- read.csv("DEG/Output/DEG_results_19422.csv")

# Créer un vecteur contenant les nouveaux noms de ligne (SPOT_ID)
ref_names <- updated_mapping[,c(2,11)]

# Vérifier s'il y a des valeurs NA (signifiant des ID sans correspondance)
sum(is.na(ref_names)) #beaucoup de NA 257
table(table(ref_names$UPDATED_GENE_ID)) #de nombreux gènes sont répétées plusieur fois 

#suppresion des lignes contenant des Na 
ref_names <- ref_names[!is.na(ref_names$UPDATED_GENE_ID), ]

df_final <- merge(DEG_19422, ref_names, by.x = "X", by.y = "SPOT_ID", all.x = TRUE)

sum(is.na(df_final$UPDATED_GENE_ID))#1211  

#beaucoup trops de NA 1/10, tentative de fusions avant DEG 
#En comparant avec les table des GPL, à la main il ne s'agit presque exclusivement de spot lié à des Na dans gene_id 

df_final_gene = df_final[!is.na(df_final$UPDATED_GENE_ID),]
sum(is.na(df_final_gene$UPDATED_GENE_ID))#0 la suppresion des na c'est bien passé 

# Écrire le fichier CSV
write.csv(df_final_gene, file = "Interaction/Input/DEG_updated_gene.csv", row.names = TRUE)


#=============================================================================================================
#Conversion des fichiers avant DEG - inutilisé pour le moment 

data_CONTROL_19422 <- read.csv("Pre-traitement/Output/DATA19422_CONTROL.csv")
data_PCC_19422 <- read.csv("Pre-traitement/Output/DATA19422_PCC.csv")

df_PCC_final <- merge(data_PCC_19422, ref_names, by.x = "X", by.y = "SPOT_ID", all.x = TRUE)
df_CONTROL_final <- merge(data_CONTROL_19422, ref_names, by.x = "X", by.y = "SPOT_ID", all.x = TRUE)


