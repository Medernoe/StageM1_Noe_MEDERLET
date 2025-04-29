#!/usr/bin/env python3
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
- DEG_TF.csv : fichier contenant les gene id ("Gene_id")- DEG_miRNA.csv : fichier contenant les noms des miRNAs (colonne 1)
- DEG_GENE.csv : fichier contenant les gene id ("Gene_id") 
#- DEG_GENE.csv : fichier non implémenter pour le moment 


Output : Data/Interaction/Output/
- interaction_TF_TF.json 
- interaction_TF_lncRNA.json 
- interaction_TF_GENE.json 
- interaction_TF_miRNA.json   #fonction non implémenter pour le moment 


Output :
- Un fichier JSON nommé : interaction_TF_<nom>.json (où <nom> = sys.argv[2])
  enregistré dans le dossier : Data/Interaction/Output/


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
# Fichier contenant les TF
folder_input = "Data/Interaction/Input/"
file_data_TF = folder_input + "DEG_TF.csv"

TF = pd.read_csv(file_data_TF, dtype=str, na_filter=False, index_col=1)
#TF.drop(TF.columns[0], axis=1, inplace=True)  # Nettoyage de la première colonne (souvent un index)


#ouverture du csv des tf, lncRNA, proteins_coding, etc  
transcript_types = ["TF","GENE","lncRNA","miRNA"] #à changer ici pour rajouter les types, ecrit au format similaire que dans les fichier DEG, sauf miRNA ou interaction miR-miR impossible
#l'ouverture se fait dans la boucle plus bas 

#%%=========================== 3. Prédiction des interactions ==============================
for biotype in transcript_types : 
    
    file_data = f"DEG_{biotype}.csv"
    print(file_data)
    full_path = folder_input + file_data
    data = pd.read_csv(full_path, dtype=str, na_filter=False)
    
    list_TF = list(set(TF.iloc[:, -2]))  # Extraction unique des TF
    total_TF = len(list_TF)
    dico_interaction = {}
    
    # Création d’un ensemble contenant tous les symboles de gènes cibles
    symbol_set = set(data["Symbol"])
    
    for i, tf in enumerate(list_TF, start=1):
        dico_interaction[tf] = {}
        
        # Sous-ensemble des interactions connues pour ce TF
        subset = df[df["source"] == tf][["target", "weight"]]
        total_cible = len(subset)
    
        for ii, (index, row) in enumerate(subset.iterrows(), start=1):
            gene = row["target"]
            sens = row["weight"]
            print(f"TF : {tf}  - {i}/{total_TF} => cible potentielle parmis {biotype} {ii}/{total_cible}")
    
            if gene in symbol_set:
                dico_interaction[tf][gene] = sens
    
        # Retire les TF sans aucune interaction conservée
        if not dico_interaction[tf]:
            del dico_interaction[tf]
    
    #=========================== 4. Sauvegarde en JSON ========================================
    # Sauvegarde du dictionnaire d'interactions dans un fichier JSON
    full_name = f"interaction_TF_{biotype}.json"
    folder_output = "Data/Interaction/Output/"
      
    # Vérifier si le répertoire output existe
    if not os.path.exists(folder_output):
        os.makedirs(folder_output)
      
      
    with open(folder_output + full_name, "w") as f:
        json.dump(dico_interaction, f)


print("Prédiction terminée pour TF! ")
