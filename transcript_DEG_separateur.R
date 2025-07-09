#!/usr/bin/env Rscript


#Description :
#Ce script lit un fichier contenant des gènes différentiellement exprimés (DEG),
#classés par type (TF, protein_coding, lncRNA, etc.), puis sépare ces gènes
#selon leur biotype pour créer des fichiers dédiés à chaque catégorie.

#Input :
#- Data/Transcript/Output/DEG_final.csv : contient les colonnes 'Symbol' et 'Type'

#Output :
#- Data/Interaction/Input/DEG_TF.csv        : facteurs de transcription (Type == 'TF')
#- Data/Interaction/Input/DEG_GENE.csv      : gènes codants classiques (Type == 'protein_coding')
#- Data/Interaction/Input/DEG_lncRNA.csv    : ARN longs non codants (Type == 'lncRNA')


#Remarque :
#Les entrées où le champ 'Symbol' est vide sont retirées (n = 146 ici).
#possibilité de rajouter les autres types de transcripts trouvé 



# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
})
# 1) Lecture des arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript transcript_DEG_separateur.R 
        <DEG_final_gene.csv> 
        <DEG_GENE.csv>
        <DEG_TF.csv>   
        <DEG_lncRNA.csv>
        <output_directory>"
       )
}
#====================================================================================================
# Chemin d'entrée du fichier contenant les gènes différentiellement exprimés
input_DEG_final   <- args[1]
# Chemin de sortie pour les fichiers séparés
output_DEG_GENE   <- args[2] 
output_DEG_TF     <- args[3]
output_DEG_lncRNA <- args[4]
 
# Répertoire de sortie
output_directory <- args[5]

#====================================================================================================
# Ouvrir le fichier contenant les gènes différentiellement exprimés
data = read.csv(input_DEG_final)

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
if (!dir.exists(output_directory)) {
  # Créer le répertoire s'il n'existe pas
  dir.create(output_directory, recursive = TRUE)
}


write.csv(DEG_TF, file = output_DEG_TF, row.names = TRUE)

write.csv(DEG_GENE, file = output_DEG_GENE, row.names = TRUE)

write.csv(DEG_lncRNA, file = output_DEG_lncRNA, row.names = TRUE)

print("fin de la séparation des fichiers de DEG")

