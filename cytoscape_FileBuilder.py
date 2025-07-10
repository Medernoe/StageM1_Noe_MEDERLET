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
#- Data/Interaction/Input/DEG_lncRNA.csv    lncRNA exprimés
- Data/Interaction/Input/DEG_miRNA_updated.csv     miRNA exprimés
- Data/Conversion/Output/dico_id.json      dictionnaire GeneID => Symbol
- Tous les fichiers JSON présents dans : Data/Interaction/Output/

Sorties :
- Data/Cytoscape/Input/nodes.csv  : liste des nœuds avec leur type (TF, miRNA, etc.)
- Data/Cytoscape/Input/edge.csv   : liste des arêtes (interactions), direction, sens, score

-Data/Enrichissement/Input
    GO.csv
    GO_TF.csv
    GO_GENE.csv
    GO_miRNA.csv
"""

import pandas as pd
import os
import json
import sys

#%% =========================== NODES ===========================

input_DEG_miRNA_updated = sys.argv[1]
input_DEG_TF = sys.argv[2]  
input_DEG_GENE = sys.argv[3]
input_dico_id = sys.argv[4]



df_tf     = pd.read_csv(input_DEG_TF,     usecols=["Symbol", "Type", "logFC","Gene_id"])
df_gene   = pd.read_csv(input_DEG_GENE ,   usecols=["Symbol", "Type", "logFC","Gene_id"])
#df_lncRNA = pd.read_csv(folder_in + "DEG_lncRNA.csv", usecols=["Symbol", "Type", "logFC","Gene_id"])

# Cas particulier : le fichier miRNA a une structure différente
df_miRNA = pd.read_csv(input_DEG_miRNA_updated)
# Ajoute une colonne Type = "miRNA"
df_miRNA["Type"] = "miRNA"
# Ne garder que les colonnes utiles
df_miRNA = df_miRNA[["Symbol", "Type", "logFC"]]


# Concatène les 3/4 types de nœuds dans un seul DataFrame
#df_concat = pd.concat([df_tf, df_gene, df_miRNA, df_lncRNA], ignore_index=True)
df_nodes_concat = pd.concat([df_tf, df_gene, df_miRNA], ignore_index=True)


#Selectionner les lignes dupliqué en fonction de la plus grosse dérégulation
# Ce processus permet de ne garder qu'une seule ligne par symbole, en prenant la ligne avec la plus grande valeur absolue de logFC
# Peut etre remplacé par la moyenne, la médiane, ou le minimum 
# ajouter une colonne temporaire avec la valeur absolue
df_nodes_concat['abs_logFC'] = df_nodes_concat['logFC'].abs()
# Pour chaque symbole, récupérer l'index de la plus petite valeur absolue, idxmax ou idxmin
idx = df_nodes_concat.groupby('Symbol')['abs_logFC'].idxmax()
# ne conserver que ces lignes et supprimer la colonne temporaire
df_nodes_unique = df_nodes_concat.loc[idx]
df_nodes_unique = df_nodes_unique.drop(columns='abs_logFC')
df_nodes_final = pd.DataFrame.reindex(df_nodes_unique, columns=["Symbol", "Type", "logFC","Gene_id"])


#retrait de la colonne gene_id inutile pour l'export des nodes 
df_nodes_export = df_nodes_final.iloc[:,0:3]
# Suppression du préfixe "hsa-miR-"
df_nodes_export["Symbol"] = df_nodes_export["Symbol"].str.replace("hsa-miR-", "", regex=False)


#export
#verification de l'existence du dossier de sortie
output_cyto = sys.argv[18]
if not os.path.exists(output_cyto):
    os.makedirs(output_cyto)

# Export du fichier de nœuds pour Cytoscape    
output_nodes = sys.argv[12]    
df_nodes_export.to_csv(output_nodes, index=False)


#%%% ==================ENRICHISSEMENT==================================

#Creation des différent csv pour l'enrichissement     
# Sélectionner les lignes où la colonne Symbol vaut 'TF', 'miRNA, ou/et 'protein_coding'
df_GO = df_nodes_final[(df_nodes_final['Type'] == 'TF') 
                             | 
                             (df_nodes_final['Type'] == 'protein_coding') 
                             ]

df_GO_TF = df_nodes_final[df_nodes_final['Type'] == 'TF'] 
df_GO_GENE = df_nodes_final[df_nodes_final['Type'] == 'protein_coding'] 
df_GO_miRNA = df_nodes_final[df_nodes_final['Type'] == 'miRNA'] 


#export
output_GO = sys.argv[14]
output_GO_TF = sys.argv[15]
output_GO_GENE = sys.argv[16]
output_GO_miRNA = sys.argv[17]

# vérification existence du répertoire de l'enrichissement fonctionnel 
output_enri = sys.argv[19]
if not os.path.exists(output_enri):
    os.makedirs(output_enri)

#contient les gènes et TF pour l'enrichissement fonctionnel
df_GO.to_csv(output_GO, index=False)
#contient les TF pour l'enrichissement fonctionnel
df_GO_TF.to_csv(output_GO_TF, index=False)
#contient les gènes codants pour l'enrichissement fonctionnel
df_GO_GENE.to_csv(output_GO_GENE, index=False)
#contient les miRNA pour l'enrichissement fonctionnel basés sur leurs cibles
df_GO_miRNA.to_csv(output_GO_miRNA, index=False)


#%% =========================== EDGE (1) : chargement des fichiers JSON ===========================

input_interaction_miRNA_TF = sys.argv[5]
input_interaction_miRNA_GENE = sys.argv[6]      
input_interaction_miRNA_lncRNA = sys.argv[7]

input_interaction_TF_TF = sys.argv[8]       
input_interaction_TF_GENE = sys.argv[9]
input_interaction_TF_lncRNA = sys.argv[10]
input_interaction_TF_miRNA = sys.argv[11]

#%% =========================== EDGE (2) : interactions miRNA → autres ===========================

#list des fichiers d'interactions (json) 
miRNA_json_files = [
    input_interaction_miRNA_lncRNA,
    input_interaction_miRNA_GENE,
    input_interaction_miRNA_TF
]

df_list = []
#parcour les fichiers d'interactions 
for json_file in miRNA_json_files:
    with open(json_file, "r") as f:
        dico = json.load(f)

    data = []

    #un mir et sa cible 
    for mir_name, targets in dico.items():
        # Pour chaque cible du miRNA, on crée une ligne dans le DataFrame
        # Le score est inversé pour que 1 soit le meilleur score (plus on est proche de 1, mieux c'est)
        for gene_id, infos in targets.items():
            data.append({
                "Source": mir_name,
                "Target": gene_id,
                "score": 1-infos["score"], #inverse la tendance, plus on est proche de 1, plus le score est bon 
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
input_dico_id = sys.argv[4]
with open(input_dico_id, "r") as f:
    dico_id = json.load(f)

def convert_id_to_symbol(gene_id):
    return dico_id.get(str(gene_id), {}).get("symbole", "NA")

df_interaction_miR["Target"] = df_interaction_miR["Target"].apply(convert_id_to_symbol)

#%% =========================== EDGE (3) : interactions TF → autres ===========================

TF_json_files = [
    input_interaction_TF_TF,
    input_interaction_TF_GENE,
    input_interaction_TF_lncRNA,
    input_interaction_TF_miRNA
]

df_list = []
for json_file in TF_json_files:
    with open(json_file, "r") as f:
        dico = json.load(f)
    data = []
    for TF_name, targets in dico.items():
        for gene_id, sens in targets.items():
            data.append({
                "Source": TF_name,
                "Target": gene_id,
                "score": 1,
                "sens" : sens,
                "validation": 1
            })
    df_list.append(pd.DataFrame(data))

df_interaction_TF = pd.concat(df_list, ignore_index=True)

#%% =========================== EDGE (4) : fusion finale et export ===========================

# Concatène interactions TF et miRNA
df_edges_concat = pd.concat([df_interaction_miR, df_interaction_TF], ignore_index=True)

# Export du fichier final d’arêtes
output_edges = sys.argv[13] 

df_edges_concat.to_csv(output_edges, index=False)

print("Fichiers nodes.csv et edges.csv générés avec succès dans :", output_cyto)
