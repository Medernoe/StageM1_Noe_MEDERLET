#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script : generate_cytoscape_nodes_edges.py

Description :
Ce script génère deux fichiers (`nodes.csv` et `edge.csv`) destinés à l"import dans Cytoscape
pour visualiser un réseau d"interactions régulatrices entre :
- miRNA → gènes, TF, lncRNA
- TF → gènes, TF, lncRNA

Les informations sont issues de fichiers CSV (DEG classés par type) et de fichiers JSON
contenant les interactions miRNA/TF générées en amont du pipeline.

Entrées :
- Data/Interaction/Input/DEG_TF.csv        TF exprimés
- Data/Interaction/Input/DEG_GENE.csv      gènes codants exprimés
- Data/Interaction/Input/DEG_lncRNA.csv    lncRNA exprimés
- Data/Interaction/Input/DEG_miRNA.csv     miRNA exprimés
- Data/Conversion/Output/dico_id.json      dictionnaire GeneID => Symbol
- Tous les fichiers JSON présents dans : Data/Interaction/Output/

Sorties :
- Data/Cytoscape/Input/nodes.csv  : liste des nœuds avec leur type (TF, miRNA, etc.)
- Data/Cytoscape/Input/edge.csv   : liste des arêtes (interactions), direction, sens, score


"""

import pandas as pd
import os
import json

#%% =========================== NODES ===========================

# Lecture des fichiers d"entrée (DEG par type)
folder_in = "Data/Interaction/Input/"
df_tf     = pd.read_csv(folder_in + "DEG_TF.csv",     usecols=["Symbol", "Type", "logFC"])
df_gene   = pd.read_csv(folder_in + "DEG_GENE.csv",   usecols=["Symbol", "Type", "logFC"])
#df_lncRNA = pd.read_csv(folder_in + "DEG_lncRNA.csv", usecols=["Symbol", "Type", "logFC"])

# Cas particulier : le fichier miRNA a une structure différente
df_miRNA = pd.read_csv(folder_in + "DEG_miRNA.csv")
# Suppression du préfixe "hsa-miR-"
df_miRNA["Symbol"] = df_miRNA["Symbol"].str.replace("hsa-miR-", "", regex=False)
# Ajoute une colonne Type = "miRNA"
df_miRNA["Type"] = "miRNA"
# Ne garder que les colonnes utiles
df_miRNA = df_miRNA[["Symbol", "Type", "logFC"]]



# Concatène les 3/4 types de nœuds dans un seul DataFrame
#df_concat = pd.concat([df_tf, df_gene, df_miRNA, df_lncRNA], ignore_index=True)
df_nodes_concat = pd.concat([df_tf, df_gene, df_miRNA], ignore_index=True)


#Selectionner les lignes duplliqué en fonction de la plus grosse dérégulation 
# ajouter une colonne temporaire avec la valeur absolue
df_nodes_concat['abs_logFC'] = df_nodes_concat['logFC'].abs()
# Pour chaque symbole, récupérer l'index de la plus petite valeur absolue, idxmax ou idxmin
idx = df_nodes_concat.groupby('Symbol')['abs_logFC'].idxmax()
# ne conserver que ces lignes et supprimer la colonne temporaire
df_nodes_unique = df_nodes_concat.loc[idx]
df_nodes_unique = df_nodes_unique.drop(columns='abs_logFC')
df_nodes_final = pd.DataFrame.reindex(df_nodes_unique, columns=["Symbol", "Type", "logFC"])


# Export du fichier de nœuds pour Cytoscape
folder_out = "Data/Cytoscape/Input/"
if not os.path.exists(folder_out):
    os.makedirs(folder_out)
     
     
df_nodes_final.to_csv(folder_out + "nodes.csv", index=False)


#%%% Creation des différent tableau pour l'enrichissement 


# Répertoire de l'enrichissement fonctionnel 
folder_GO = "Data/Enrichissement/Input/"
if not os.path.exists(folder_GO):
    os.makedirs(folder_GO)
    
# Sélectionner les lignes où la colonne Symbol vaut 'TF', 'miRNA, ou/et 'protein_coding'
df_GO = df_nodes_final[(df_nodes_final['Type'] == 'TF') 
                             | 
                             (df_nodes_final['Type'] == 'protein_coding') 
                             ]

df_GO_TF = df_nodes_final[df_nodes_final['Type'] == 'TF'] 
df_GO_GENE = df_nodes_final[df_nodes_final['Type'] == 'protein_coding'] 
df_GO_miRNA = df_nodes_final[df_nodes_final['Type'] == 'miRNA'] 


#export
df_GO.to_csv(folder_GO + "GO.csv", index=False)
df_GO_TF.to_csv(folder_GO + "GO_TF.csv", index=False)
df_GO_GENE.to_csv(folder_GO + "GO_GENE.csv", index=False)
df_GO_miRNA.to_csv(folder_GO + "GO_miRNA.csv", index=False)


#%% =========================== EDGE (1) : chargement des fichiers JSON ===========================

folder_json = "Data/Interaction/Output/"
json_contents = {}

# Chargement de tous les fichiers .json du dossier
for file_name in os.listdir(folder_json):
    if file_name.endswith(".json"):
        with open(os.path.join(folder_json, file_name), "r", encoding="utf-8") as f:
            json_contents[file_name] = json.load(f)

#%% =========================== EDGE (2) : interactions miRNA → autres ===========================

#list des fichiers d'interactions (json) 
miRNA_json_files = [
    #"interaction_miRNA_lncRNA.json",
    "interaction_miRNA_GENE.json",
    "interaction_miRNA_TF.json"
]

df_list = []
#parcour les fichiers d'interactions 
for json_file in miRNA_json_files:
    dico = json_contents.get(json_file, {})
    data = []
    #un mir et sa cible 
    for mir_name, targets in dico.items():
        
        for gene_id, infos in targets.items():
            data.append({
                "Target": gene_id,
                "Source": mir_name,
                #"score": infos["score"],
                "sens" : -1,
                "validation": infos["val"]
            })
    df_list.append(pd.DataFrame(data))

# Concatène toutes les interactions miRNA
df_interaction_miR = pd.concat(df_list, ignore_index=True)

# Nettoyage du nom des miRNA
df_interaction_miR["Source"] = df_interaction_miR["Source"].str.replace("hsa-miR-", "", regex=False)



"""
Peut etre modifier ici le processus de conversion ou plus tot et ne travailler que sur les gene symbol 

"""
# Conversion des gene_id en symboles via dico_id
with open("Data/Conversion/Output/dico_id.json", "r") as f:
    dico_id = json.load(f)

def convert_id_to_symbol(gene_id):
    return dico_id.get(str(gene_id), {}).get("symbole", "NA")

df_interaction_miR["Target"] = df_interaction_miR["Target"].apply(convert_id_to_symbol)

#%% =========================== EDGE (3) : interactions TF → autres ===========================

TF_json_files = [
    #"interaction_TF_lncRNA.json",
    "interaction_TF_GENE.json",
    "interaction_TF_TF.json", 
    "interaction_TF_miRNA.json"
]

df_list = []
for json_file in TF_json_files:
    dico = json_contents.get(json_file, {})
    data = []
    for TF_name, targets in dico.items():
        for gene_id, sens in targets.items():
            data.append({
                "Target": gene_id,
                "Source": TF_name,
                #"score": None,
                "sens" : sens,
                "validation": 1
            })
    df_list.append(pd.DataFrame(data))

df_interaction_TF = pd.concat(df_list, ignore_index=True)

#%% =========================== EDGE (4) : fusion finale et export ===========================

# Concatène interactions TF et miRNA
df_edges_concat = pd.concat([df_interaction_miR, df_interaction_TF], ignore_index=True)

# Export du fichier final d’arêtes
df_edges_concat.to_csv(folder_out + "edges.csv", index=False)

print("Fichiers nodes.csv et edges.csv générés avec succès dans :", folder_out)
