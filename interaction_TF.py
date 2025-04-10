#!/usr/bin/env python3
"""
Script : predict_TF_gene_interactions.py

Description :
Ce script prédit les interactions **TF → gène cible** à partir des données CollecTRI (via decoupler),
en croisant :
- Une liste de TF différentiellement exprimés (DEG_TF.csv)
- Un fichier de gènes cibles (par ex. protein_coding, lncRNA…)

Les interactions sont conservées uniquement si la cible est présente dans le fichier d'entrée.

Entrées :
- "Data/Interaction/Input/DEG_TF.csv" : fichier CSV contenant les TF (colonne "Symbol")
- sys.argv[1] : fichier CSV contenant les cibles potentielles (colonne "Symbol")

Sortie :
- Un fichier JSON nommé : interaction_TF_<nom>.json (où <nom> = sys.argv[2])
  enregistré dans le dossier : Data/Interaction/Output/

Exemple :
$ python predict_TF_gene_interactions.py Data/Interaction/Input/DEG_lncRNA.csv lncRNA

Résultat : Data/Interaction/Output/interaction_TF_lncRNA.json
"""

import decoupler as dc
import pandas as pd
import sys
import json

#=========================== 1. Chargement du dataset CollecTRI ===========================
df = dc.get_collectri(organism='human', split_complexes=False)

#=========================== 2. Ouverture des fichiers ====================================
# Fichier contenant les TF
file_data_TF = "Data/Interaction/Input/DEG_TF.csv"
TF = pd.read_csv(file_data_TF, dtype=str, na_filter=False, index_col=1)
TF.drop(TF.columns[0], axis=1, inplace=True)  # Nettoyage de la première colonne (souvent un index)

# Fichier des gènes cibles passé en argument
file_data = sys.argv[1]
data = pd.read_csv(file_data, dtype=str, na_filter=False)
data.drop(data.columns[0], axis=1, inplace=True)

#=========================== 3. Prédiction des interactions ==============================
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
        print(f"TF : {tf}  - {i}/{total_TF} => cible potentielle {ii}/{total_cible}")

        if gene in symbol_set:
            dico_interaction[tf][gene] = sens

    # Retire les TF sans aucune interaction conservée
    if not dico_interaction[tf]:
        del dico_interaction[tf]

#=========================== 4. Sauvegarde en JSON ========================================
name = sys.argv[2]
full_name = f"interaction_TF_{name}.json"
folder = "Data/Interaction/Output/"

with open(folder + full_name, "w") as f:
    json.dump(dico_interaction, f)

print(f"✅ Prédiction terminée ! Fichier enregistré ici : {folder + full_name}")
