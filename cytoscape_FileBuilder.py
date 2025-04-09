import pandas as pd
import os
import json

#%%NODES
#%%%ouverture des fichier
folder_in = "Data/Interaction/Input/" 
# Charger les CSV en sélectionnant uniquement les colonnes 'Symbole' et 'type'
df_tf = pd.read_csv(folder_in + "DEG_TF.csv", usecols=["Symbol", "Type"])
df_gene = pd.read_csv(folder_in + "DEG_GENE.csv", usecols=["Symbol", "Type"])
df_lncRNA = pd.read_csv(folder_in + "DEG_lncRNA.csv", usecols=["Symbol", "Type"])

#cas spécial pour miRNA qui possèdent un DF strucutré différemment 
df_miRNA = pd.read_csv(folder_in + "DEG_miRNA.csv")

#Renomme la colonne 0 en "Symbol" 
position = 1
old_name = df_miRNA.columns[position]
df_miRNA.rename(columns={old_name: 'Symbol'}, inplace=True)

#Retirer le préfixe "hsa-miR-" 
df_miRNA['Symbol'] = df_miRNA['Symbol'].str.replace("hsa-miR-", "", regex=False)

#Créer la colonne "type" et lui attribuer la valeur "miRNA"
df_miRNA['Type'] = 'miRNA'

#garder que ces deux colonnes (Symbole et type)
df_miRNA = df_miRNA[['Symbol', 'Type']]

print(df_miRNA.head())
print(df_tf.head())
print(df_gene.head())
print(df_lncRNA.head())
#%%%création du dataframe de nodes 

# Concaténer tous les DataFrames dans un seul
df_concat = pd.concat([df_tf, df_gene, df_miRNA, df_lncRNA], ignore_index=True)

# Supprimer les doublons en se basant sur la colonne 'Symbol'
df_unique = df_concat.drop_duplicates(subset="Symbol")

folder_out = "Data/Cytoscape/Input/"
# Sauvegarder le DataFrame résultant dans un nouveau CSV
df_unique.to_csv(folder_out + "nodes.csv", index=False)



#%%EDGE 
#%%%ouverture des Json 

# Chemin vers le dossier où se trouvent tes JSON
folder_in = "Data/Interaction/Output/"

# Créer un dictionnaire pour stocker le contenu de chaque fichier
json_contents = {}

# Parcourir tous les fichiers du dossier
for file_name in os.listdir(folder_in):
    # Vérifier si le nom de fichier se termine par ".json"
    if file_name.endswith(".json"):
        file_path = os.path.join(folder_in, file_name)
        # Ouvrir et charger le fichier JSON
        with open(file_path, "r", encoding="utf-8") as f:
            data = json.load(f)
            # Placer le contenu dans le dictionnaire, indexé par le nom du fichier
            json_contents[file_name] = data

#%%%création d'un tableau pour miRNA 

#fichiers JSON
miRNA_gene = "interaction_miRNA_gene.json"
miRNA_lncRNA = "interaction_miRNA_lncRNA.json"
miRNA_TF     = "interaction_miRNA_TF.json" 

types = [miRNA_gene, miRNA_lncRNA, miRNA_TF]

# Liste pour stocker les DataFrame intermédiaires
df_list = []

for json_file in types:
    # Récupère le dictionnaire associé au fichier
    dico = json_contents.get(json_file, "")
    
    # On construit une liste de lignes (chaque ligne = un dict)
    data = []
    for mir_name, gene_dict in dico.items():
        for gene_id, gene_info in gene_dict.items():
            data.append({
                "Source": mir_name,
                "Target": gene_id,
                "score": gene_info["score"],
                "sens" : (-1), 
                "validation": gene_info["val"]
            })
    
    # Création du DataFrame temporaire
    df_temp = pd.DataFrame(data)
    # Ajout du DataFrame à la liste
    df_list.append(df_temp)

# Concaténation de tous les DataFrame en un seul
df_interaction_miR = pd.concat(df_list, ignore_index=True)

#Retirer le préfixe "hsa-miR-" 
df_interaction_miR['Source'] = df_interaction_miR['Source'].str.replace("hsa-miR-", "", regex=False)

# Affichage du DataFrame final
print(df_interaction_miR)

#%%%%conversion des gene_id en symbol 
with open("Data/Conversion/Output/dico_id.json", "r") as f:
    dico_id = json.load(f)

# Fonction de conversion avec le dictionnaire dico_id
def convert_id_to_symbol(gene_id):
    # On s'assure que gene_id est une chaîne
    gene_id_str = str(gene_id)
    return dico_id.get(gene_id_str, {}).get("symbole", "NA")

# Application de la conversion sur chaque ligne de la colonne 'Target'
df_interaction_miR["Target"] = df_interaction_miR["Target"].apply(convert_id_to_symbol)

    
#%%%création d'un tableau pour TF 
TF_gene   = "interaction_TF_gene.json"
TF_lncRNA = "interaction_TF_lncRNA.json"
TF_TF     = "interaction_TF_TF.json"

types = [TF_gene, TF_lncRNA, TF_TF]

# Liste pour stocker les DataFrame intermédiaires
df_list = []

for json_file in types:
    # Récupère le dictionnaire associé au fichier
    dico = json_contents.get(json_file, "")
    
    # On construit une liste de lignes (chaque ligne = un dict)
    data = []
    for TF_name, gene_dict in dico.items():
        # Ici, gene_dict est un dictionnaire du type {"BCL2": "-1", "CCND2": "1", ...}
        for gene_id, sens in gene_dict.items():
            data.append({
                "Source": TF_name,
                "Target": gene_id,
                "score": None, 
                "sens" : sens, 
                "validation": None
            })
    
    
    # Création du DataFrame temporaire
    df_temp = pd.DataFrame(data)
    # Ajout du DataFrame à la liste
    df_list.append(df_temp)

# Concaténation de tous les DataFrames en un seul
df_interaction_TF = pd.concat(df_list, ignore_index=True)

# Affichage du DataFrame final
print(df_interaction_TF)

duplicates = df_interaction_TF[df_interaction_TF.duplicated()]
print(duplicates)
#aucune ligne n'est dupliqué 

#%%%création et export du df final de EDGE 

# Concaténer tous les DataFrames dans un seul
df_concat = pd.concat([df_interaction_miR, df_interaction_TF], ignore_index=True)

folder_out = "Data/Cytoscape/Input/"
# Sauvegarder le DataFrame résultant dans un nouveau CSV
df_concat.to_csv(folder_out + "edge.csv", index=False)

