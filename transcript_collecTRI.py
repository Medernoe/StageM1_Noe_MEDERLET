"""
Ce script met à jour la classification des gènes dans un fichier de résultats (TFlink),
en utilisant les données de collecTRI pour distinguer les vrais facteurs de transcription (TF)
des autres gènes codants (protein_coding).

Input : Data/Transcript/Output/   DEG_complete_TFlink.csv
Output : Data/Transcript/Output/   DEG_complete_collecTRI.csv
"""

import decoupler as dc
import pandas as pd 

#%% Téléchargement des données collecTRI depuis decoupler (organisme humain)
df = dc.get_collectri(organism='human', split_complexes=False)

# Affiche les 5 premières lignes pour visualiser le contenu
print(df.head())

#%% Extraction des noms de gènes considérés comme TF dans collecTRI
TF = list(df["source"])  # Contient les symboles de gènes identifiés comme TF

#%% Chargement du fichier produit par TFlink
folder = "Data/Transcript/Output/"
DEG = pd.read_csv(folder + "DEG_complete_TFlink.csv", index_col=1)

# Suppression de la première colonne inutile (souvent un index doublon)
DEG.drop(DEG.columns[0], axis=1, inplace=True)

#%% Mise à jour du type de chaque gène selon collecTRI
i = 0
for Type in DEG["Type"]:  # Parcourt chaque ligne de la colonne "Type"
    
    # Ne traite que les types TF ou protein_coding
    if Type == "TF" or Type == "protein_coding": 
        
        symbol = DEG.iloc[i, 7]  # Colonne contenant le symbole du gène

        if symbol in TF:
            DEG.iloc[i, 8] = "TF"  # Si le symbole est dans collecTRI, c'est un vrai TF
        else:
            DEG.iloc[i, 8] = "protein_coding"  # Sinon, on le considère comme un gène codant

    i += 1

#%% Sauvegarde du nouveau fichier mis à jour
DEG.to_csv(folder + "DEG_complete_collecTRI.csv", index=1)


