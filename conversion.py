from Bio import Entrez
import json 
import pandas as pd

"""
Ce programme est à usage unique, il sert à mettre à jour une list d'identifiant gene, et à récuperer les symboles associés 
"""

def extract_accessions(data):
    """
    Récupère les numéros d'accessions refseq (NM_, XM_, NR_, XR_).
    """
    accessions = []  # Liste pour stocker les accessions récupérées
    
    # Vérifier si l'argument 'data' est un dictionnaire
    if isinstance(data, dict):
        # Parcours de toutes les clés et valeurs dans le dictionnaire
        for key, value in data.items():
            # Si la clé est 'Gene-commentary_accession', cela indique un accession
            if key == 'Gene-commentary_accession':  
                
                # Récupérer l'accession de la valeur associée à la clé
                accession = value 
                
                # Vérifier si l'accession est une chaîne de caractères valide et commence par l'un des préfixes valides
                if isinstance(accession, str) and (accession.startswith('NM_') 
                                                   or accession.startswith('XM_') 
                                                   or accession.startswith('NR_') 
                                                   or accession.startswith('XR_')):
                    
                    # Vérifier si l'accession n'est pas déjà présente dans la liste
                    if accession not in accessions: 
                        print(accession)
                        accessions.append(accession)  # Ajouter l'accession à la liste
                
            # Si la valeur est un dictionnaire ou une liste, on appelle la fonction récursivement pour chercher dans les sous-données
            elif isinstance(value, (dict, list)):
                accessions.extend(extract_accessions(value))  # Ajouter les résultats obtenus récursivement
                
    # Si la donnée 'data' est une liste
    elif isinstance(data, list):
        # Parcourir chaque élément de la liste et appeler la fonction récursivement
        for item in data:
            accessions.extend(extract_accessions(item))  # Ajouter les résultats obtenus de chaque élément

    # Retourner la liste des accessions récupérées
    return accessions


def update(gene_id, dico_uptade):
    """
    Réalise une requête API au NCBI pour vérifier la version du gene_id.
    Si l'ID est obsolète, il le met à jour et récupère les informations du gène,
    y compris le symbole en limitant les requêtes.
    """
    
    #-------------------------------------------------------------------------------------------------- 
    #test qui évite de faire une requete inutile 
    if gene_id in dico_uptade: 
        print(f"Le GeneID {gene_id} et son symbole ont déjà été mis à jour.")
        symbole = dico_uptade[gene_id]
        return gene_id, symbole
    
    #-------------------------------------------------------------------------------------------------- 
    #variables
    Entrez.email = "noe.mederlet@univ-rouen.fr"
    max_attempts = 5 # fixe le seuil après lequel les tentatives de requetes s'arretent 
    attempts = 0
    records = None
      
    #-------------------------------------------------------------------------------------------------- 
    #Récupération des données 
    
    while attempts < max_attempts:
        try:
            handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            break  
        
        except Exception as e:
            attempts += 1
            print(f"Erreur lors de la récupération (tentative {attempts}/{max_attempts}) : {e}")
            if attempts < max_attempts:
                print("Nouvelle tentative...")
            else:
                print("Toutes les tentatives ont échoué.")
                records = None

    if not records:
        print(f"Aucune correspondance pour l'ID {gene_id}.")
        return None


    #-------------------------------------------------------------------------------------------------- 
    # Récupération des données du gène
    
    gene_data = records[0]
    track_info = gene_data.get("Entrezgene_track-info", {})
    gene_track = track_info.get("Gene-track", {})
    
    # Extraction du status du gène
    status = gene_track.get("Gene-track_status", None)
    status_value = getattr(status, "attributes", {}).get("value", None) #attributes == la ou sont est stockée le status du gène valide/obsolète 
    
    # Extraction du symbole 
    symbole = str(gene_data.get("Entrezgene_gene", {})
                  .get("Gene-ref", {})
                  .get("Gene-ref_formal-name", {})
                  .get("Gene-nomenclature", {})
                  .get("Gene-nomenclature_symbol", "N/A"))
    
    #-------------------------------------------------------------------------------------------------- 
    # Récupération des données du gène
    if status_value == "secondary":
        # L'ID a été remplacé par un nouvel ID
        current_ids = gene_track.get("Gene-track_current-id", [])
        for dbtag in current_ids:
            if dbtag.get("Dbtag_db") == "GeneID":
                new_gene_id = dbtag["Dbtag_tag"]["Object-id"]["Object-id_id"]
                break #arret après la récupération du nouveau gene id 

        if new_gene_id != gene_id:
            print(f"Le GeneID {gene_id} est obsolète et a été remplacé par {new_gene_id}.")
            # Mise à jour du GeneID et récupération des informations du gène (symbole, description, etc.)
            gene_id = new_gene_id
            
            print("Récupération des informations misent à jour ...")
            
            #--------------------------------------------------------------------------------------------------
            #refait une requete avec le gene_id à jour pour obtenir le symbole à jour 
            
            attempts = 0 #mise à jour de la variable attempts 
            
            while attempts < max_attempts:
                try:
                    handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
                    records = Entrez.read(handle)
                    handle.close()
                    break  
                
                except Exception as e:
                    attempts += 1
                    print(f"Erreur lors de la récupération (tentative {attempts}/{max_attempts}) : {e}")
                    if attempts < max_attempts:
                        print("Nouvelle tentative...")
                    else:
                        print("Toutes les tentatives ont échoué.")
                        records = None
            
            if records : 
                gene_data = records[0]

                # Extraction du symbole 
                symbole = str(gene_data.get("Entrezgene_gene", {})
                              .get("Gene-ref", {})
                              .get("Gene-ref_formal-name", {})
                              .get("Gene-nomenclature", {})
                              .get("Gene-nomenclature_symbol", "N/A"))
                

            else: 
                symbole = None 
          
                
        else: 
            print(f"Le GeneID {gene_id} est valide.")
            
            
    else:
        print(f"Le GeneID {gene_id} est valide.")
    
    print(f"Symbole du gène : {symbole}")
    
    return gene_id, symbole 

    

#%% Import de la table de jointure et création des variables 

folder = "C:/Users/noeme/OneDrive/Desktop/Stage_M1/Data/"
file_data = "merged_mapping_export.csv"

df = pd.read_csv(folder + file_data, dtype=str, na_filter=False)
 


# Conversion des colonnes en string
total = len(df["GENE_4133"])


#dictionnaire de controle 
controle = { "arret" : 0  , #peut etre utilisé pour reprendre si le programme plante 
             "identique" : 0,
             "different" : 0,
             "old" : 0,
             "new" : 0, 
             "old_maj" : 0, 
             "new_maj" : 0,
             "none" : 0
             }

#list utiliser pour stocker les nouelles version de gene_id 
id_updated = []  

#id et symbole sert à eviter de repeter des requetes API 
dico_id_symbole = {}


#%%boucle qui parcours lignes à lignes le dataframe 

for i in range(0,total): 

    print("-"*50)
    print(f"gene id numéro {i}/{total}")
    controle["arret"] = i
    
    #-----------------------------------------------------------------------------------
    #cas ou les 2 gènes sont identique 
    
    if df["GENE_4133"].iloc[i] == df["GENE_6480"].iloc[i]: 
        
        controle["identique"] +=1 
        print("Identiques")
        
        gene_id, symbole = update(df["GENE_4133"].iloc[i],dico_id_symbole)
        
        if gene_id != None : 
            id_updated.append(gene_id)
        
        else:
            id_updated.append(None)
        
    
    #-----------------------------------------------------------------------------------
    #cas ou les gènes sont différents 
    
    elif df["GENE_4133"].iloc[i] != df["GENE_6480"].iloc[i]:
        
        controle["different"] += 1  
        print("Différents")

        
        #-----------------------------------------------------------------------------------
        #GENE_4133 vaut NA et GENE_6480 à une valeur 
        
        if df["GENE_4133"].iloc[i] == "NA" and df["GENE_6480"].iloc[i] != "NA": 
            
            print("NA/ID")

            gene_id, symbole = update(df["GENE_6480"].iloc[i],dico_id_symbole)
          
            if gene_id != None : 
                controle["old"] += 1
                id_updated.append(gene_id)
            
            else:
                controle["none"] += 1
                id_updated.append(None)
                
        
        #-----------------------------------------------------------------------------------
        #GENE_4133 à une valeur vaut NA et GENE_6480 vaut NA 

        elif df["GENE_4133"].iloc[i] != "NA" and df["GENE_6480"].iloc[i] == "NA":
            
            gene_id, symbole = update(df["GENE_4133"].iloc[i],dico_id_symbole)
            
            print("ID/NA")

            if gene_id != None : 
                controle["new"] += 1
                id_updated.append(gene_id)
            
            else:
                controle["none"] += 1
                id_updated.append(None)
                
        
        #-----------------------------------------------------------------------------------
        #cas ou les 2 gene_id ont une valeur, mais différente
        
        else: 
            id_4133, symbole_4133 = update(df["GENE_4133"].iloc[i],dico_id_symbole)
            id_6180, symbole_6180 = update(df["GENE_6480"].iloc[i],dico_id_symbole) 
            
            print("ID/ID")

            if id_4133 == id_6180:
                id_updated.append(id_4133)
                
                gene_id = id_4133 
                symbole = symbole_4133
                
                
            elif id_4133 and id_6180 is None:
                controle["old_maj"] += 1
                id_updated.append(id_4133)
                
                gene_id = id_4133 
                symbole = symbole_4133
                
                
            elif id_4133 is None and id_6180:
                controle["new_maj"] += 1
                id_updated.append(id_6180)
                
                gene_id = id_6180 
                symbole = symbole_6180
                
                
            else:
                controle["none"] += 1
                id_updated.append(None)
                symbole = None
        
    #ajoute les valeurs au dico 
    dico_id_symbole[gene_id] = symbole
    
    
    #-----------------------------------------------------------------------------------
    #Sauvegarde tout les 1000 gene_id 
    
    if i % 1000 == 0:
        print(f"Sauvegarde de la liste id_updated après {i} gènes.")
        
        # Sauvegarder des variables temporaires dans des fichiers afin de garder une trace en cas d'erreur 
        with open(folder + f"id_updated_{i}.json", "w") as f:
            json.dump(id_updated, f)
            
        with open(folder + f"dico_controle_{i}.json", "w") as f:
            json.dump(controle, f)
    
        with open(folder + f"dico_id_symbole_{i}.json", "w") as f:
            json.dump(dico_id_symbole, f)
            
   
        
#%%export des datas 

#ajouter la list dans une nouvelle colonne du tableau  
df["UPDATED_GENE_ID"] = id_updated 

#Exporter le nouveau dataframe 
output_path = "update_gene_mapping.csv"
df.to_csv(folder + output_path, index=False)
    







