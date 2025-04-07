#!/usr/bin/env Rscript

# Définir le répertoire de travail
setwd("/home/medernoe/Ubuntu_user/Stage_M1_ubuntu/StageM1_git/Data")

data = read.csv("Transcript/Output/DEG_complete_collecTRI.csv")

table(data$Type)

sum(data$Symbol == "")
#les 146 unknown ne correspondent pas à des gènes 


#====================================================================================================
#sépare les données pour leur export 
data = data[data$Symbol != "", ]

DEG_TF = data[data$Type == "TF", ]

DEG_GENE = data[data$Type == "protein_coding", ]

DEG_lncRNA = data[data$Type == 'lncRNA', ]

#===================================================================================================
#reécrit les fichiers 
write.csv(DEG_TF, file = "Interaction/Input/DEG_TF.csv", row.names = TRUE)

write.csv(DEG_GENE, file = "Interaction/Input/DEG_GENE.csv", row.names = TRUE)

write.csv(DEG_lncRNA, file = "Interaction/Input/DEG_lncRNA.csv", row.names = TRUE)


