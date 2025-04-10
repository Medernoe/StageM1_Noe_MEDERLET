"""
Script : get_gene_biotype_from_ensembl.py

Ce script interroge l'API REST d'Ensembl pour déterminer le biotype du transcrit principal 
associé à chaque gène, à partir de son symbole (ex : TP53, NEAT1, MIR21, etc.).

Il utilise aussi une liste externe de facteurs de transcription humains pour reclassifier 
les gènes codants en "TF" si applicable.

Input :
- DEG_updated_gene.csv : fichier contenant une colonne "Symbol" avec les noms des gènes
- Homo_sapiens_TF.txt : fichier tabulé contenant une liste de symboles de facteurs de transcription

Output :
- DEG_complete.csv : même tableau que l'entrée, avec une colonne ajoutée "Type" indiquant :
    • "TF" pour les facteurs de transcription
    • "protein_coding" pour les gènes codants
    • "miRNA", "lncRNA", etc., selon la classification Ensembl
    • "Unknown" si le symbole n’est pas trouvé

Utilité :
Ce script permet de catégoriser automatiquement des gènes différentiellement exprimés selon 
leur fonction ou type de transcrit principal.

Attention :
Il est necessaire de modifier ce code en le fusionnant avec transcript_collecTRI.py CF: ligne 90 
"""


import requests
import pandas as pd 

#%% Fonction pour interroger l'API d'Ensembl à partir d'un symbole de gène
def requete_data(symbole = ""):
    """
    Fait une requête à l'API REST d'Ensembl avec un symbole de gène (ex: TP53),
    et retourne un dictionnaire contenant les informations associées.
    """
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/homo_sapiens/{symbole}?expand=1"
    url = server + ext

    attempts = 0
    max_attempts = 5
    dico = None 
    
    # Essaye jusqu'à 5 fois en cas d'échec
    while attempts < max_attempts:
        if dico is None:
            attempts += 1 
            try:
                r = requests.get(url, headers={"Content-Type": "application/json"})
                r.raise_for_status()  # Lève une erreur si code != 200
                dico = r.json()  # Réponse au format JSON (dictionnaire)
            except requests.exceptions.RequestException as e:
                print(f"Erreur lors de la récupération de données : {e}")
                dico = None
        else:
            attempts = max_attempts  # Sortie anticipée si la requête a réussi
            
    return dico 

#%% Fonction pour isoler le type de gène ("biotype") à partir de l'API
def isole_data(symbole):
    """
    Isole le biotype du gène depuis les données récupérées par l’API.
    Par ex : "protein_coding", "miRNA", "lncRNA", etc.
    """
    dico_API = requete_data(symbole)
    
    if dico_API:
        biotype = dico_API.get("biotype", None)
    else: 
        biotype = "Unknown"  # Si la requête a échoué
    
    return biotype


#%%
# Dossiers d'entrée et sortie
folder = "Data/Transcript/Input/"
folder_output = "Data/Transcript/Output/"

# Lecture du fichier des gènes différentiellement exprimés (DEG)
gene = pd.read_csv(folder + "DEG_updated_gene.csv")

#%% Lecture du fichier contenant les facteurs de transcription humains (TF)

"""
Il est necessaire de modifier ici la symbol_list par les TF porvenant de collecTRI

"""
with open(folder + "Homo_sapiens_TF.txt", "r", encoding="cp1252") as file:
    lines = file.readlines()

# Création d'une liste des symboles de TF à partir de la deuxième colonne (index 1)
symbol_list = [line.split("\t")[1].strip() for line in lines[1:] if line.split("\t")[1].strip()]

#%% Création d'une colonne "Type" en déterminant le biotype pour chaque symbole
biotypes = []
i = 0
total = len(gene["Symbol"])

for symbole in gene["Symbol"]:
    print(f"Symbole {i}/{total}")
    biotype = isole_data(symbole)

    # Vérifie si c'est un TF dans notre liste
    if biotype == "protein_coding":
        if symbole in symbol_list:
            biotype = "TF"

    biotypes.append(biotype)
    print("=>" + biotypes[i])
    i += 1

#%% Enregistre les résultats dans un nouveau fichier CSV
gene["Type"] = biotypes
gene.to_csv(folder_output + "DEG_complete.csv", index=False)
