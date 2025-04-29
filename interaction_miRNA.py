
"""
Description :
Ce script interroge une base de données MySQL locale (mirabel) pour récupérer les interactions connues
entre des miRNAs différentiellement exprimés (DEG) et une liste de gènes (DEG aussi). Il utilise les IDs
"mimat" associés aux miRNAs pour retrouver les interactions dans les tables SQL, puis trie les résultats 
pour ne garder que les meilleures interactions (top 5% par score).

Input : Data/Interaction/Input/
- DEG_miRNA.csv : fichier contenant les noms des miRNAs (colonne 1)
- DEG_lncRNA.csv : fichier contenant les gene id ("Gene_id")
- DEG_TF.csv : fichier contenant les gene id ("Gene_id")- DEG_miRNA.csv : fichier contenant les noms des miRNAs (colonne 1)
- DEG_GENE.csv : fichier contenant les gene id ("Gene_id")

Output : Data/Interaction/Output/
- interaction_miRNA_TF.json 
- interaction_miRNA_lncRNA.json 
- interaction_miRNA_GENE.json 

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

#%% Connexion à la base de données MySQL locale
mydb = mysql.connector.connect(
  host="localhost",
  user="medernoe", #ces informations sont propres à la création d'un utilisateur mysql 
  password="1234", 
)

#%% Fonction : récupération du mimat à partir d'un nom de miRNA
def miR_name_mimat(miR_name):
    mycursor = mydb.cursor()
    mycursor.execute("USE mirabel") 
    
    mycursor.execute(f"""SELECT miR_mimat.number, mimat.number, predictions2.mimat
                    FROM miR
                    JOIN miR_mimat ON miR.id = miR_mimat.mir_id
                    JOIN mimat ON miR_mimat.number = mimat.number
                    JOIN predictions2 ON mimat.number = predictions2.mimat
                    WHERE miR.nameR = "{miR_name}"
                    LIMIT 1;""")
    
    myresult = mycursor.fetchall()
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

#%% Chargement des fichiers CSV
folder_input = "Data/Interaction/Input/"

#ouverture du csv des miRNA 
file_data_miRNA = folder_input + "DEG_miRNA.csv"
miRNA = pd.read_csv(file_data_miRNA, dtype=str, na_filter=False)


#ouverture du csv des tf, lncRNA, proteins_coding, etc  
transcript_types = ["TF","GENE","lncRNA"] #à changer ici pour rajouter les types, ecrit au format similaire que dans les fichier DEG, sauf miRNA ou interaction miR-miR impossible

#%% Boucle principale : pour chaque miRNA, associer les interactions aux gènes

for biotype in transcript_types : 
    
    file_data = f"DEG_{biotype}.csv"
    full_path = folder_input + file_data
    data = pd.read_csv(full_path, dtype=str, na_filter=False)


    list_mimat = []
    total_miR = len(miRNA["Symbol"])
    total = len(data["Gene_id"])
    i = 0
    dico_interaction = {}
    
    for mirname in miRNA["Symbol"]: 
        i += 1
        mimat = miR_name_mimat(mirname)
        
        if not mimat:
            print(f"Aucun mimat trouvé pour {mirname}")
            continue
    
        mimat = str(mimat[0][0])  # on extrait juste l'identifiant
        list_mimat.append(mimat)
        dico_interaction[mirname] = {}
    
        ii = 0
        for gene in data["Gene_id"]:
            ii += 1
            print(f"miRNA : {mirname}  - {i}/{total_miR} => {biotype} {ii}/{total}")
    
            if gene not in dico_interaction[mirname]:
                prediction = gene_mirna_interaction(gene, mimat)
                
                if prediction:
                    dico_interaction[mirname][gene] = {
                        "score": prediction[0][0],
                        "val": str(prediction[0][1])
                    }
    
        # On garde seulement le top X% des interactions
        dico_interaction[mirname] = pourcent_trie(dico_interaction[mirname])
    
    #%% Sauvegarde du dictionnaire d'interactions dans un fichier JSON
    full_name = f"interaction_miRNA_{biotype}.json"
    folder_output = "Data/Interaction/Output/"
    
    # Vérifier si le répertoire output existe
    if not os.path.exists(folder_output):
        os.makedirs(folder_output)
    
    
    with open(folder_output + full_name, "w") as f:
        json.dump(dico_interaction, f)
  
                     
            