#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description :
Ce script prédit les interactions TF → gène cible à partir des données CollecTRI (via decoupler),
en croisant :
- Une liste de TF différentiellement exprimés (DEG_TF.csv)
- Un fichier de gènes cibles (par ex. protein_coding, lncRNA…)

Les interactions sont conservées uniquement si la cible est présente dans le fichier d'entrée.


Input : Data/Interaction/Input/
- DEG_TF.csv : fichier CSV contenant les TF (colonne "Symbol")
- DEG_lncRNA.csv : fichier contenant les gene id ("Gene_id")
- DEG_miRNA.csv : fichier contenant les noms des miRNAs (colonne 1)
- DEG_GENE.csv : fichier contenant les gene id ("Gene_id") 


Output : Data/Interaction/Output/
- interaction_TF_TF.json 
- interaction_TF_lncRNA.json 
- interaction_TF_GENE.json 
- interaction_TF_miRNA.json   #fonction non implémenter pour le moment 


Remarques :
- Une version utlétrieur fonctionne avec sys.argv : Commits on Apr 15, 2025
=> Exemple :
    $ python predict_TF_gene_interactions.py Data/Interaction/Input/DEG_lncRNA.csv lncRNA
=> Résultat : 
    Data/Interaction/Output/interaction_TF_lncRNA.json
    
"""

import decoupler as dc
import pandas as pd
import sys
import json
import os 

#%%=========================== 1. Chargement du dataset CollecTRI ===========================
df = dc.get_collectri(organism='human', split_complexes=False)

#%%=========================== 2. Ouverture des fichiers ====================================
# ouverture du fichier des TF
input_DEG_TF = sys.argv[1]  # Chemin vers le fichier d'entrée des TF
TF = pd.read_csv(input_DEG_TF, dtype=str, na_filter=False, index_col=1)
#TF.drop(TF.columns[0], axis=1, inplace=True)  # Nettoyage de la première colonne (un index)


#ouverture du csv des miRNA, lncRNA, proteins_coding, etc  
input_DEG_GENE = sys.argv[2]  # Chemin vers le fichier d'entrée des gènes
input_DEG_lncRNA = sys.argv[3]  # Chemin vers le fichier d'entrée des lncRNAs
input_DEG_miRNA = sys.argv[4]  # Chemin vers le fichier d'entrée des TF

# Liste des fichiers d'entrée pour les gènes
input_files = [input_DEG_TF, input_DEG_GENE, input_DEG_lncRNA, input_DEG_miRNA] 


#ouverture des fichiers de sortie  # Chemin vers le fichier de sortie des miRNAs mis à jour
output_interaction_TF_TF = sys.argv[5]  # Chemin vers le fichier de sortie des interactions miRNA-TF
output_interaction_TF_GENE = sys.argv[6]  # Chemin vers le fichier de sortie des interactions miRNA-GENE
output_interaction_TF_lncRNA = sys.argv[7]  # Chemin vers le fichier de sortie des interactions miRNA-lncRNA
output_interaction_TF_miRNA = sys.argv[8]  # Chemin vers le fichier de sortie des interactions miRNA-lncRNA


# Liste des fichiers de sortie correspondants
output_files = [output_interaction_TF_TF, 
                 output_interaction_TF_GENE, 
                 output_interaction_TF_lncRNA, 
                 output_interaction_TF_miRNA]


# Création de la liste de tuples (entrée, sortie)
input_output_pairs = list(zip(input_files, output_files))
# input_output_pairs = [(input_DEG_TF, output_interaction_miRNA_TF), (input_DEG_GENE, output_interaction_miRNA_GENE), (input_DEG_lncRNA, output_interaction_miRNA_lncRNA)]


# Répertoire de sortie
folder_output = sys.argv[9] 

# Vérifier si le répertoire output existe
if not os.path.exists(folder_output):
    os.makedirs(folder_output)

#%%=========================== 3. Prédiction des interactions ==============================
for pairs in input_output_pairs : 
    
    input_file, output_file = pairs

    # ouverture du fichier CSV de DEG
    data = pd.read_csv(input_file, dtype=str, na_filter=False)
    
    list_TF = list(set(TF.iloc[:, -2]))  # Extraction unique des TF
    dico_interaction = {}
    
    # Création d’un ensemble contenant tous les symboles de gènes cibles
    symbol_set = set(data["Symbol"])
    
    for tf in list_TF:
        dico_interaction[tf] = {}
        
        # Sous-ensemble des interactions connues pour ce TF
        subset = df[df["source"] == tf][["target", "weight"]]
        
    
        for (index, row) in subset.iterrows():
            gene = row["target"]
            sens = row["weight"]
            
            if gene in symbol_set:
                dico_interaction[tf][gene] = sens
    
        # Retire les TF sans aucune interaction conservée
        if not dico_interaction[tf]:
            del dico_interaction[tf]
    
    #=========================== 4. Sauvegarde en JSON ========================================
    # Sauvegarde du dictionnaire d'interactions dans un fichier JSON
    with open(output_file, "w") as f:
        json.dump(dico_interaction, f)


print("Prédiction terminée pour TF! ")
