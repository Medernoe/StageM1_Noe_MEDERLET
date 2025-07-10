#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Description :
Ce script interroge une base de données MySQL locale (mirabel) pour récupérer les interactions connues
entre des miRNAs différentiellement exprimés (DEG) et une liste de gènes (DEG aussi). Il utilise les IDs
"mimat" associés aux miRNAs pour retrouver les interactions dans les tables SQL, puis trie les résultats 
pour ne garder que les meilleures interactions (top 5% par score).
Met à jour les noms des miR avec la version 21 de MirBase. 


Input : Data/Interaction/Input/
- DEG_miRNA.csv : fichier contenant les noms des miRNAs (colonne 1)
- DEG_lncRNA.csv : fichier contenant les gene id ("Gene_id")
- DEG_TF.csv : fichier contenant les gene id ("Gene_id")- DEG_miRNA.csv : fichier contenant les noms des miRNAs (colonne 1)
- DEG_GENE.csv : fichier contenant les gene id ("Gene_id")


Output : Data/Interaction/Output/
- interaction_miRNA_TF.json 
- interaction_miRNA_lncRNA.json 
- interaction_miRNA_GENE.json 
# Data/Interaction/Input/  DEG_miRNA.csv : version mise à jours des miR

Fonctions principales :
- `miR_name_mimat()` : récupère le mimat (identifiant Ensembl) à partir du nom du miRNA
- `gene_mirna_interaction()` : récupère les scores d'interaction entre un gène (gene_id) et un mimat
- `pourcent_trie()` : trie et sélectionne le top 5% des interactions les plus fortes


Remarques :
- Ce script suppose que la base de données MySQL "mirabel" est installée localement.
- Utilisation prévue dans un pipeline d'analyse d'interactions miRNA-gène.
- Une version utlétrieur fonctionne avec sys.argv : Commits on Apr 15, 2025
=> Exemple :
$ python predict_miRNA_gene_interactions.py Data/Interaction/Input/DEG_lncRNA.csv lncRNA
"""

import mysql.connector
import pandas as pd 
import json 
import os 
import sys 

USER = sys.argv[10]  # Identifiant de l'utilisateur MySQL
PASSWORD = sys.argv[11]  # Mot de passe de l'utilisateur MySQL

#%% Connexion à la base de données MySQL locale
mydb = mysql.connector.connect(
  host = "localhost",
  user = USER, #ces informations sont propres à la création d'un utilisateur mysql 
  password= PASSWORD, 
)


#%% Fonction : récupération du mimat à partir d'un nom de miRNA
def miR_name_mimat(miR_name):
    mycursor = mydb.cursor()
    mycursor.execute("USE mirabel") 
    
    mycursor.execute(f"""SELECT  miR_mimat.number
                     FROM miR_mimat
                     JOIN miR ON miR_mimat.mir_id = miR.id
                     WHERE miR.nameR = "{miR_name}"
                     LIMIT 1;""")
                    
    mimat = mycursor.fetchall()[0][0]
    
    mycursor.close()
    return mimat

#%% Fonction : de mije à jours du mimat avec la version 21 de MirBase
def mirname_updated_by_mimat(mimat) :
    mycursor = mydb.cursor()
    mycursor.execute("USE mirabel") 
    
    mycursor.execute(f"""SELECT miR.nameR
                     FROM miR
                     JOIN miR_mimat ON miR_mimat.mir_id = miR.id
                     WHERE miR_mimat.number = "{mimat}"
                     AND miR_mimat.MBversion = "21";""") #21 est la derniere version de MirBase utilisé dans mirabel. 
    
    myresult = mycursor.fetchall()[0][0]
    
    mycursor.close()
    
    return myresult

#%% Fonction : récupération d'une interaction gene-miRNA
def gene_mirna_interaction(gene_id, miR_mimat):
    mycursor = mydb.cursor()
    mycursor.execute("USE mirabel") 
    
    mycursor.execute(f"""SELECT predictions2.score, predictions2.VALID  
                    FROM predictions2 
                    JOIN genes3 ON genes3.ncbi_gene_id = predictions2.id_g 
                    WHERE genes3.ncbi_gene_id = "{gene_id}"
                    AND predictions2.mimat = "{miR_mimat}"
                    LIMIT 12;""")
    
    myresult = mycursor.fetchall()
    mycursor.close()
    return myresult

#%% Fonction : retourne les X% meilleurs scores d'un dictionnaire d'interactions
def pourcent_trie(dico, pourcentage = 0.05): 
    dico_trie = dict(sorted(dico.items(), key=lambda item: item[1]['score'], reverse=False)) 
    n = len(dico_trie)
    top = int(n * pourcentage) or 1
    return dict(list(dico_trie.items())[:top])  


#========================================
#%% ouverture des fichiers d'entrée et de sortie

# ouverture du fichier des miRNAs
input_DEG_miRNA = sys.argv[1]  # Chemin vers le fichier d'entrée des miRNAs
miRNA = pd.read_csv(input_DEG_miRNA, dtype=str, na_filter=False)


#ouverture du csv des tf, lncRNA, proteins_coding, etc  
input_DEG_TF = sys.argv[2]  # Chemin vers le fichier d'entrée des TF
input_DEG_GENE = sys.argv[3]  # Chemin vers le fichier d'entrée des gènes
input_DEG_lncRNA = sys.argv[4]  # Chemin vers le fichier d'entrée des lncRNAs
# Liste des fichiers d'entrée pour les gènes
input_files = [input_DEG_TF, input_DEG_GENE, input_DEG_lncRNA] 


#ouverture des fichiers de sortie  # Chemin vers le fichier de sortie des miRNAs mis à jour
output_interaction_miRNA_TF = sys.argv[6]  # Chemin vers le fichier de sortie des interactions miRNA-TF
output_interaction_miRNA_GENE = sys.argv[7]  # Chemin vers le fichier de sortie des interactions miRNA-GENE
output_interaction_miRNA_lncRNA = sys.argv[8]  # Chemin vers le fichier de sortie des interactions miRNA-lncRNA
# Liste des fichiers de sortie correspondants
output_files = [output_interaction_miRNA_TF, output_interaction_miRNA_GENE, output_interaction_miRNA_lncRNA]


# Création de la liste de tuples (entrée, sortie)
input_output_pairs = list(zip(input_files, output_files))
# input_output_pairs = [(input_DEG_TF, output_interaction_miRNA_TF), (input_DEG_GENE, output_interaction_miRNA_GENE), (input_DEG_lncRNA, output_interaction_miRNA_lncRNA)]


# Répertoire de sortie
folder_output = sys.argv[9] 

# Vérifier si le répertoire output existe
if not os.path.exists(folder_output):
    os.makedirs(folder_output)



#%%=========================================
# Boucle principale : pour chaque type de gène (TF, GENE, lncRNA)

for pairs in input_output_pairs : 
    # Décomposition des tuples en variables distinctes
    input_file, output_file = pairs

    # ouverture du fichier CSV de DEG
    data = pd.read_csv(input_file, dtype=str, na_filter=False)

    list_mirname = []
    list_mimat = []
    dico_interaction = {}
    
    for mirname in miRNA["Symbol"]: 
        # On récupère le mimat à partir du nom du miRNA
        mimat = miR_name_mimat(mirname)
        mirname = mirname_updated_by_mimat(mimat)
        
        if not mimat:
            print(f"Aucun mimat trouvé pour {mirname}")
            continue
        
        list_mirname.append(mirname)
        #list_mimat.append(mimat)
        dico_interaction[mirname] = {}

        for gene in data["Gene_id"]:
            
            if gene not in dico_interaction[mirname]:
                #prédiction est une liste de tuples (score, val)
                prediction = gene_mirna_interaction(gene, mimat)
                
                if prediction:
                    dico_interaction[mirname][gene] = {
                        "score": prediction[0][0],
                        "val": str(prediction[0][1])
                    }
    
        # On garde seulement le top X% des interactions
        dico_interaction[mirname] = pourcent_trie(dico_interaction[mirname])
    
    with open(output_file, "w") as f:
        json.dump(dico_interaction, f)

#%%==========================================
# enregistrement du fichier DEG avec les noms des miR à jours 

output_DEG_miRNA_updated = sys.argv[5]
miRNA["Symbol"] = list_mirname 
miRNA.to_csv(output_DEG_miRNA_updated, index=False)            


