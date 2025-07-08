#!/usr/bin/env Rscript

#Ce programme remplace les spot_id par les gene_id en usant des gene_id mis à jour à partir de conversion_oldID_to_newID

#Input : Data/Conversion/Output/   update_gene_mapping.csv
#        Data/DEG/Output/          DEG_results_gene.csv
#        Data/Conversion/Output/   dico_id.json
#Output : Data/Transcript/Input/   DEG_updated_gene.csv

#informations : une version possédant l'étude des jeux de données et le formatage est disponible dans : Commits on Apr 15, 2025

suppressPackageStartupMessages({
  library(jsonlite)
})
# 1) Lecture des arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript conversion_SpotToGene.R 
        <update_gene_mapping.csv> 
        <dico_id.json> 
        <DEG_results_gene.csv> 
        <DEG_updated_gene.csv>
        <output_directory>"
       )
}

#table de conversion contenant les gene id/spot id
in_updated_mapping   <- args[1]
#fichier JSON
in_data_symbole <- args[2]
#DEG
in_DEG_gene    <- args[3]
#chemin de sortie du fichier DEG mis à jour
out_DEG_updated   <- args[4]
#repertoire de sortie
output_directory <- args[5]


#ouvrir les fichier-------------------------------------------------------------------------------------------
# Lire le fichier JSON
data_symbole <- fromJSON(in_data_symbole)

#ouvrir la table contenant les gene id/spot id 
updated_mapping <- read.csv(in_updated_mapping)

#ouvrir la DEG
DEG_gene <- read.csv(in_DEG_gene)



#ajouter la colonne des gene_id----------------------------------------------------------------------------
# Créer un vecteur contenant les nouveaux noms de ligne (SPOT_ID)
ref_names <- updated_mapping[,c(2,11)]

#suppresion des lignes contenant des Na 
ref_names <- ref_names[!is.na(ref_names$UPDATED_GENE_ID), ]


#fusion de la DEG et de la table de mapping 
df_1 <- merge(DEG_gene, ref_names, by.x = "X", by.y = "SPOT_ID", all.x = TRUE)

df_2 <- df_1[!duplicated(df_1), ]

df_final_gene = df_2[!is.na(df_2$UPDATED_GENE_ID),]


#modification du nom de la colonne 
colnames(df_final_gene)[colnames(df_final_gene) == "UPDATED_GENE_ID"] <- "Gene_id"
colnames(df_final_gene)


#ajout des gene_symbol-----------------------------------------------------------------------------------------------------
# Ajouter une colonne "Symbol" remplie de NA
df_final_gene$Symbol <- NA

# Boucle sur chaque ligne du dataframe
for (i in 1:nrow(df_final_gene)) {
  gene_id <- as.character(df_final_gene$Gene_id[i])  # Récupérer le gene_id
  
  # Vérifier si gene_id existe dans le JSON
  if (gene_id %in% names(data_symbole)) {
    symbol_value <- data_symbole[[gene_id]]$symbole  # Extraire la valeur associée
    
    # Vérifier si la valeur est non vide
    if (!is.null(symbol_value) && symbol_value != "") {
      df_final_gene$Symbol[i] <- symbol_value  # Remplacer la valeur
    }
  } else {
    print(paste("Erreur : aucune clé associée au gene_id", gene_id)) #Aucun gene symbole n'existe pour ce gene_id 
  }
}

# Écrire le fichier CSV------------------------------------------------------------------------------------------
# Vérifier si le répertoire Output existe
if (!dir.exists(output_directory)) {
  # Créer le répertoire s'il n'existe pas
  dir.create(output_directory, recursive = TRUE)
}

write.csv(df_final_gene, file = out_DEG_updated, row.names = TRUE)

print("fin de la conversion des anciens id avec les nouveaux.")



