#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Ce script interroge l'API REST d'Ensembl pour determiner le biotype du transcrit principal 
associe Ã  chaque gene, a  partir de son symbole (ex : TP53, NEAT1, MIR21, etc.).

Il utilise aussi une liste externe provenant de collecTRI de facteurs de transcription humains 
pour reclassifier les genes codants en "TF" si applicable.

Input :
- Data/Transcript/Input/  DEG_updated_gene.csv : fichier contenant une colonne "Symbol" avec les noms des genes
- requete decoupler collecTRI

Output :
- Data/Transcript/Output/ DEG_final.csv : meme tableau que l'entree, avec une colonne ajoutee "Type" indiquant :
     "TF" pour les facteurs de transcription
     "protein_coding" pour les genes codants
     "miRNA", "lncRNA", etc., selon la classification Ensembl
     "Unknown" si le symbole nâest pas trouve

"""


import requests #permet de faire les requetes API ensembl
import pandas as pd 
import decoupler as dc #permet de recuperer les datas de collecTRI
import os 
import sys


#%% Fonction pour interroger l'API d'Ensembl Ã  partir d'un symbole de gene
def requete_data(symbole = ""):
    """
    Fait une requete Ã  l'API REST d'Ensembl avec un symbole de gene (ex: TP53),
    et retourne un dictionnaire contenant les informations associees.
    """
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/homo_sapiens/{symbole}?expand=1"
    url = server + ext

    attempts = 0
    max_attempts = 5
    dico = None 
    
    # Essaye jusqu'Ã  5 fois en cas d'echec
    while attempts < max_attempts:
        if dico is None:
            attempts += 1 
            try:
                r = requests.get(url, headers={"Content-Type": "application/json"})
                r.raise_for_status()  # Leve une erreur si code != 200
                dico = r.json()  # Reponse au format JSON (dictionnaire)
            except requests.exceptions.RequestException as e:
                print(f"Erreur lors de la recuperation de donnees : {e}")
                dico = None
        else:
            attempts = max_attempts  # Sortie anticipee si la requete a reussi
            
    return dico 

#%% Fonction pour isoler le type de gene ("biotype") Ã  partir de l'API
def isole_data(symbole):
    """
    Isole le biotype du gene depuis les donnees recuperees par l'API.
    Par ex : "protein_coding", "miRNA", "lncRNA", etc.
    """
    dico_API = requete_data(symbole)
    
    if dico_API:
        biotype = dico_API.get("biotype", None)
    else: 
        biotype = "Unknown"  # Si la requete a echoue
    
    return biotype



# Lecture du fichier des genes differentiellement exprimes (DEG)
input_file = sys.argv[1]
gene = pd.read_csv(input_file)


#%% Utilisation de collecTRI afin de determiner les TF humain 
#Telechargement des donnees collecTRI depuis decoupler (organisme humain)
df = dc.get_collectri(organism='human', split_complexes=False)

# Extraction des noms de genes consideres comme TF dans collecTRI
TF_list = list(df["source"]) 

#%% Creation d'une colonne "Type" en determinant le biotype pour chaque symbole

print("Debut de la determination des types de transcript.")
#initialise la list qui servira Ã  creer la colonne type 
biotypes = []
i = 0
total = len(gene["Symbol"])

for symbole in gene["Symbol"]:
    print(f"Symbole {i}/{total}")
    biotype = isole_data(symbole)

    # Verifie si c'est un TF dans notre liste
    if biotype == "protein_coding":
        if symbole in TF_list:
            biotype = "TF"
            
    #remplie la list
    biotypes.append(biotype)
    i += 1

#%% Enregistre les resultats dans un nouveau fichier CSV
gene["Type"] = biotypes

#%%
folder_output = sys.argv[3] 
# Vérifier si le répertoire output existe
if not os.path.exists(folder_output):
    os.makedirs(folder_output)

#%% Export final
output_file = sys.argv[2]
gene.to_csv(output_file, index=False)
