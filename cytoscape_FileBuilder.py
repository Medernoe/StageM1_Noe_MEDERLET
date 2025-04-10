#!/usr/bin/env python3
"""
Script : generate_cytoscape_nodes_edges.py

Description :
Ce script génère deux fichiers (`nodes.csv` et `edge.csv`) destinés à l'import dans Cytoscape
pour visualiser un réseau d'interactions régulatrices entre :
- miRNA → gènes, TF, lncRNA
- TF → gènes, TF, lncRNA

Les informations sont issues de fichiers CSV (DEG classés par type) et de fichiers JSON
contenant les interactions miRNA/TF générées en amont du pipeline.

Entrées :
- Data/Interaction/Input/DEG_TF.csv       ← TF exprimés
- Data/Interaction/Input/DEG_GENE.csv     ← gènes codants exprimés
- Data/Interaction/Input/DEG_lncRNA.csv   ← lncRNA exprimés
- Data/Interaction/Input/DEG_miRNA.csv    ← miRNA exprimés
- Data/Conversion/Output/dico_id.json     ← dictionnaire GeneID → Symbol
- Tous les fichiers JSON présents dans : Data/Interaction/Output/

Sorties :
- Data/Cytoscape/Input/nodes.csv  : liste des nœuds avec leur type (TF, miRNA, etc.)
- Data/Cytoscape/Input/edge.csv   : liste des arêtes (interactions), direction, sens, score

Utilisation :
Ce script s'intègre dans un pipeline de visualisation de réseaux régulateurs sous Cytoscape.
"""

import pandas as pd
import os
import json

#%% =========================== NODES ===========================

# Lecture des fichiers d'entrée (DEG par type)
folder_in = "Data/Interaction/Input/"
df_tf     = pd.read_csv(folder_in + "DEG_TF.csv",     usecols=["Symbol", "Type"])
df_gene   = pd.read_csv(folder_in + "DEG_GENE.csv",   usecols=["Symbol", "Type"])
df_lncRNA = pd.read_csv(folder_in + "DEG_lncRNA.csv", usecols=["Symbol", "Type"])

# Cas particulier : le fichier miRNA a une structure différente
df_miRNA = pd.read_csv(folder_in + "DEG_miRNA.csv")

# On renomme la colonne de position 1 en "Symbol"
df_miRNA.rename(columns={df_miRNA.columns[1]: "Symbol"}, inplace=True)

# Suppression du préfixe "hsa-miR-"
df_miRNA['Symbol'] = df_miRNA['Symbol'].str.replace("hsa-miR-", "", regex=False)

# Ajoute une colonne Type = "miRNA"
df_miRNA['Type'] = 'miRNA'

# Ne garder que les colonnes utiles
df_miRNA = df_miRNA[['Symbol', 'Type']]

# Concatène les 4 types de nœuds dans un seul DataFrame
df_concat = pd.concat([df_tf, df_gene, df_miRNA, df_lncRNA], ignore_index=True)

# Supprime les doublons par symbole
df_unique = df_concat.drop_duplicates(subset="Symbol")

# Export du fichier de nœuds pour Cytoscape
folder_out = "Data/Cytoscape/Input/"
df_unique.to_csv(folder_out + "nodes.csv", index=False)

#%% =========================== EDGE (1) : chargement des fichiers JSON ===========================

folder_json = "Data/Interaction/Output/"
json_contents = {}

# Chargement de tous les fichiers .json du dossier
for file_name in os.listdir(folder_json):
    if file_name.endswith(".json"):
        with open(os.path.join(folder_json, file_name), "r", encoding="utf-8") as f:
            json_contents[file_name] = json.load(f)

#%% =========================== EDGE (2) : interactions miRNA → autres ===========================

miRNA_json_files = [
    "interaction_miRNA_gene.json",
    "interaction_miRNA_lncRNA.json",
    "interaction_miRNA_TF.json"
]

df_list = []
for json_file in miRNA_json_files:
    dico = json_contents.get(json_file, {})
    data = []
    for mir_name, targets in dico.items():
        for gene_id, infos in targets.items():
            data.append({
                "Source": mir_name,
                "Target": gene_id,
                "score": infos["score"],
                "sens" : -1,
                "validation": infos["val"]
            })
    df_list.append(pd.DataFrame(data))

# Concatène toutes les interactions miRNA
df_interaction_miR = pd.concat(df_list, ignore_index=True)

# Nettoyage du nom des miRNA
df_interaction_miR['Source'] = df_interaction_miR['Source'].str.replace("hsa-miR-", "", regex=False)

# Conversion des gene_id en symboles via dico_id
with open("Data/Conversion/Output/dico_id.json", "r") as f:
    dico_id = json.load(f)

def convert_id_to_symbol(gene_id):
    return dico_id.get(str(gene_id), {}).get("symbole", "NA")

df_interaction_miR["Target"] = df_interaction_miR["Target"].apply(convert_id_to_symbol)

#%% =========================== EDGE (3) : interactions TF → autres ===========================

TF_json_files = [
    "interaction_TF_gene.json",
    "interaction_TF_lncRNA.json",
    "interaction_TF_TF.json"
]

df_list = []
for json_file in TF_json_files:
    dico = json_contents.get(json_file, {})
    data = []
    for TF_name, targets in dico.items():
        for gene_id, sens in targets.items():
            data.append({
                "Source": TF_name,
                "Target": gene_id,
                "score": None,
                "sens" : sens,
                "validation": None
            })
    df_list.append(pd.DataFrame(data))

df_interaction_TF = pd.concat(df_list, ignore_index=True)

#%% =========================== EDGE (4) : fusion finale et export ===========================

# Concatène interactions TF et miRNA
df_concat = pd.concat([df_interaction_miR, df_interaction_TF], ignore_index=True)

# Export du fichier final d’arêtes
df_concat.to_csv(folder_out + "edge.csv", index=False)

print("✅ Fichiers nodes.csv et edge.csv générés avec succès dans :", folder_out)
