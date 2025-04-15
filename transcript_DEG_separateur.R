#!/usr/bin/env Rscript


#Description :
#Ce script lit un fichier contenant des gènes différentiellement exprimés (DEG),
#classés par type (TF, protein_coding, lncRNA, etc.), puis sépare ces gènes
#selon leur biotype pour créer des fichiers dédiés à chaque catégorie.

#Input :
#- Data/Transcript/Output/DEG_complete_collecTRI.csv : contient les colonnes 'Symbol' et 'Type'

#Output :
#- Data/Interaction/Input/DEG_TF.csv        : facteurs de transcription (Type == 'TF')
#- Data/Interaction/Input/DEG_GENE.csv      : gènes codants classiques (Type == 'protein_coding')
#- Data/Interaction/Input/DEG_lncRNA.csv    : ARN longs non codants (Type == 'lncRNA')

#Remarque :
#Les entrées où le champ 'Symbol' est vide sont retirées (n = 146 ici).

#Utilité :
#Facilite l’analyse ciblée selon le type de transcrit (ex : interactions TF-lncRNA).


# =============================================================================


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

# Vérifier si le répertoire Output existe
if (!dir.exists("Interaction/Input/")) {
  # Créer le répertoire s'il n'existe pas
  dir.create("Interaction/Input/", recursive = TRUE)
}


write.csv(DEG_TF, file = "Interaction/Input/DEG_TF.csv", row.names = TRUE)

write.csv(DEG_GENE, file = "Interaction/Input/DEG_GENE.csv", row.names = TRUE)

write.csv(DEG_lncRNA, file = "Interaction/Input/DEG_lncRNA.csv", row.names = TRUE)


