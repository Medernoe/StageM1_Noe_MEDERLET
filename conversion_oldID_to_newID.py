#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Entrez
import json 
import pandas as pd
import os


"""
Ce programme met à jour une liste d'identifiants de gènes (GeneID) en interrogeant la base NCBI.
Il vérifie si les identifiants sont encore valides, les met à jour si nécessaire, et récupère les symboles génétiques associés.

Input : Data/Conversion/Input/ 
- merged_mapping_export.csv

Output : Data/Conversion/Output/ 
- update_gene_mapping.csv : contient les spot_id, gene_id à jour 
- dico_id.json : contient pour un gene_id, son symbole ainsi que ses transcripts 
- id_updated.json : list complète de tous les id à jours 
- dico_controle.json : metadatas du nombre de version obsolète, à jours, ... Difficile à lire, correspond aux différentes voies empruntées 

"""

# Fonction récursive pour extraire les accession RefSeq d'une structure JSON
def extract_accessions(data):
    accessions = []
    if isinstance(data, dict):
        for key, value in data.items():
            if key == 'Gene-commentary_accession':  
                accession = value 
                if isinstance(accession, str) and accession.startswith(('NM_', 'XM_', 'NR_', 'XR_')):
                    if accession not in accessions: 
                        accessions.append(accession)  
            elif isinstance(value, (dict, list)):
                accessions.extend(extract_accessions(value))  
    elif isinstance(data, list):
        for item in data:
            accessions.extend(extract_accessions(item))  
    return accessions

# Fonction principale de mise à jour d'un GeneID
def update(gene_id, dico_uptade):
    if gene_id in dico_uptade:  # Évite de répéter une requête déjà faite
        print(f"Le GeneID {gene_id} et son symbole ont déjà été mis à jour.")
        return gene_id, dico_uptade[gene_id]["symbole"], dico_uptade[gene_id]["refseq"]
    
    Entrez.email = "noe.mederlet@univ-rouen.fr"
    max_attempts = 5
    attempts = 0
    records = None
    
    # Tentatives de récupération des données depuis NCBI
    while attempts < max_attempts:
        try:
            handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            break  
        except Exception as e:
            attempts += 1
            print(f"Erreur lors de la récupération (tentative {attempts}/{max_attempts}) : {e}")
    
    if not records:
        print(f"Aucune correspondance pour l'ID {gene_id}.")
        return gene_id, None, None 

    gene_data = records[0]
    gene_track = gene_data.get("Entrezgene_track-info", {}).get("Gene-track", {})
    status = gene_track.get("Gene-track_status", None)
    status_value = getattr(status, "attributes", {}).get("value", None)  # Extrait si le gène est obsolète

    # Cas où l'ID est obsolète et remplacé
    if status_value == "secondary":
        current_ids = gene_track.get("Gene-track_current-id", [])
        for dbtag in current_ids:
            if dbtag.get("Dbtag_db") == "GeneID":
                new_gene_id = dbtag["Dbtag_tag"]["Object-id"]["Object-id_id"]
                break 

        if new_gene_id != gene_id:
            print(f"Le GeneID {gene_id} est obsolète et a été remplacé par {new_gene_id}.")
            gene_id = new_gene_id
            print("Récupération des informations mises à jour ...")
            
            attempts = 0
            while attempts < max_attempts:
                try:
                    handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
                    records = Entrez.read(handle)
                    handle.close()
                    break  
                except Exception as e:
                    attempts += 1
                    print(f"Erreur lors de la récupération (tentative {attempts}/{max_attempts}) : {e}")
            
            if not records:
                return gene_id, None, None 
            else: 
                gene_data = records[0]
                symbole = str(gene_data.get("Entrezgene_gene", {})
                              .get("Gene-ref", {})
                              .get("Gene-ref_formal-name", {})
                              .get("Gene-nomenclature", {})
                              .get("Gene-nomenclature_symbol", "N/A"))
                refseq = extract_accessions(gene_data)
        else: 
            print(f"Le GeneID {gene_id} est valide.")
    else:
        print(f"Le GeneID {gene_id} est probablement valide.")
        symbole = str(gene_data.get("Entrezgene_gene", {})
                      .get("Gene-ref", {})
                      .get("Gene-ref_formal-name", {})
                      .get("Gene-nomenclature", {})
                      .get("Gene-nomenclature_symbol", "N/A"))
        refseq = extract_accessions(gene_data)
       
    print(f"Symbole du gène : {symbole}")
    return gene_id, symbole, refseq

# Fonction d'évaluation de la "qualité" d'une entrée
def scoring_values(x, y, z):
    list_data = [x, y, z]
    return sum(1 for data in list_data if data is not None)

#%% Chargement des données initiales
print("Début du programme de conversion des anciens ID au nouveaux.")

#input
folder_intput = "Data/Conversion/Input/"
file_data = "merged_mapping_export.csv"
df = pd.read_csv(folder_intput + file_data, dtype=str, na_filter=False)
total = len(df["GENE_4133"])

# Dictionnaires de contrôle/statistiques
controle = {
    "arret": 0, "identique": 0, "different": 0,
    "old": 0, "new": 0, "old_maj": 0, "new_maj": 0,
    "none": 0, "none_old": 0, "none_new": 0, "ID/ID diff": 0
}

id_updated = []
dico_id = {}

#%% Boucle principale de traitement
for i in range(0, total): 
    print("-"*50)
    print(f"gene id numéro {i}/{total}")
    controle["arret"] = i
    
    if df["GENE_4133"].iloc[i] == df["GENE_6480"].iloc[i]: 
        controle["identique"] += 1 
        print("Identiques")
        gene_id, symbole, refseq = update(df["GENE_4133"].iloc[i], dico_id) 
        id_updated.append(gene_id if gene_id is not None else None)

    elif df["GENE_4133"].iloc[i] != df["GENE_6480"].iloc[i]:
        controle["different"] += 1  
        print("Différents")

        if df["GENE_4133"].iloc[i] == "NA" and df["GENE_6480"].iloc[i] != "NA": 
            print("NA/ID")
            gene_id, symbole, refseq = update(df["GENE_6480"].iloc[i], dico_id)
            if gene_id: 
                controle["old"] += 1
            else:
                controle["none_old"] += 1
            id_updated.append(gene_id)

        elif df["GENE_4133"].iloc[i] != "NA" and df["GENE_6480"].iloc[i] == "NA":
            print("ID/NA")
            gene_id, symbole, refseq = update(df["GENE_4133"].iloc[i], dico_id)
            if gene_id: 
                controle["new"] += 1
            else:
                controle["none_new"] += 1
            id_updated.append(gene_id)

        else:
            id_4133, symbole_4133, refseq_4133 = update(df["GENE_4133"].iloc[i], dico_id)
            id_6180, symbole_6180, refseq_6180 = update(df["GENE_6480"].iloc[i], dico_id)

            refseq_4133 = refseq_4133 if refseq_4133 else None
            refseq_6180 = refseq_6180 if refseq_6180 else None

            score_4133 = scoring_values(id_4133, symbole_4133, refseq_4133)
            score_6180 = scoring_values(id_6180, symbole_6180, refseq_6180)
            
            print("ID/ID")

            if id_4133 == id_6180:
                controle["ID/ID diff"] += 1
                id_updated.append(id_4133)
            elif score_4133 > score_6180:
                controle["new_maj"] += 1
                id_updated.append(id_4133)
            elif score_4133 < score_6180:
                controle["old_maj"] += 1
                id_updated.append(id_6180)
            else:
                controle["none"] += 1
                id_updated.append(None)

    # Mise à jour du dictionnaire si le GeneID n'était pas présent
    if gene_id not in dico_id:
        dico_id[gene_id] = {"symbole": symbole, "refseq": refseq}

    # Sauvegarde temporaire tous les 1000 gènes mis en commentaire pour le pipeline nextflow
    """
    if i % 1000 == 0:
        print(f"Sauvegarde de la liste id_updated après {i} gènes.")
        with open(folder_output + f"id_updated_{i}.json", "w") as f:
            json.dump(id_updated, f)
        with open(folder_output + f"dico_controle_{i}.json", "w") as f:
            json.dump(controle, f)
        with open(folder_output + f"dico_id_{i}.json", "w") as f:
            json.dump(dico_id, f)
    """
    
#%% Export final

#output
folder_output = "Data/Conversion/Output/"
# Vérifier si le répertoire output existe
if not os.path.exists(folder_output):
    os.makedirs(folder_output)
with open(folder_output + "dico_id.json", "w") as f:
    json.dump(dico_id, f)

with open(folder_output + "id_updated.json", "w") as f:
    json.dump(id_updated, f)

with open(folder_output + "dico_controle.json", "w") as f:
    json.dump(controle, f)

df["UPDATED_GENE_ID"] = id_updated
df.to_csv(folder_output + "update_gene_mapping.csv", index=False)

print("fin du programme de conversion !!")
