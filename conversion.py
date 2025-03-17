from Bio import Entrez

import pandas as pd



#%%

def update(gene_id, list_uptade):
    """
    Réalise une requete API au NCBI, afin de vérifier la version du gene_id 
    """
    
    if gene_id in list_uptade : 
        print(f"Le GeneID {gene_id} a déja été mis à jour.")
        return gene_id
    
    Entrez.email = "noe.mederlet@univ-rouen.fr"
    
    max_attempts = 3
    attempts = 0
    records = None

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

    gene_data = records[0]
    track_info = gene_data.get("Entrezgene_track-info", {})
    gene_track = track_info.get("Gene-track", {})

    status = gene_track.get("Gene-track_status", None)
    status_value = getattr(status, "attributes", {}).get("value", None)

    if status_value == "secondary":
        current_ids = gene_track.get("Gene-track_current-id", [])
        new_gene_id = None
        for dbtag in current_ids:
            if dbtag.get("Dbtag_db") == "GeneID":
                new_gene_id = dbtag["Dbtag_tag"]["Object-id"]["Object-id_id"]
                break

        if new_gene_id:
            print(f"Le GeneID {gene_id} est obsolète et a été remplacé par {new_gene_id}.")
            return new_gene_id
        else:
            print(f"Le GeneID {gene_id} est obsolète, et aucun nouvel ID n'a été trouvé.")
            return None
    else:
        print(f"Le GeneID {gene_id} est valide.")
        return gene_id

    

#%% Import de la table de jointure 


# Vous pouvez aussi préciser le type dès la lecture du CSV pour certaines colonnes via le paramètre dtype
df = pd.read_csv("C:/Users/noeme/OneDrive/Desktop/Stage_M1/Data/merged_mapping_export.csv", dtype=str, na_filter=False)


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


id_updated = []  
#-----------------------------------------------------------------------------------
#boucle qui parcours lignes à lignes le dataframe 
for i in range(0,total): 
    
    print(f"gene id numéro {i}/{total}")
    controle["arret"] = i
    
    #-----------------------------------------------------------------------------------
    #cas ou les 2 gènes sont identique 
    
    if df["GENE_4133"].iloc[i] == df["GENE_6480"].iloc[i]: 
        
        controle["identique"] +=1 
        
        id = update(df["GENE_4133"].iloc[i],id_updated)
        
        if id != None : 
            id_updated.append(df["GENE_4133"].iloc[i])
        
        else:
            id_updated.append(None)
        
    
    #-----------------------------------------------------------------------------------
    #cas ou les gènes sont différents 
    
    elif df["GENE_4133"].iloc[i] != df["GENE_6480"].iloc[i]:
        
        controle["different"] += 1  
        
        
        #-----------------------------------------------------------------------------------
        #GENE_4133 vaut NA et GENE_6480 à une valeur 
        
        if df["GENE_4133"].iloc[i] == "NA" and df["GENE_6480"].iloc[i] != "NA": 
            
            
            id = update(df["GENE_6480"].iloc[i],id_updated)
          
            if id != None : 
                controle["old"] += 1
                id_updated.append(df["GENE_6480"].iloc[i])
            
            else:
                controle["none"] += 1
                id_updated.append(None)
                
        
        #-----------------------------------------------------------------------------------
        #GENE_4133 à une valeur vaut NA et GENE_6480 vaut NA 

        elif df["GENE_4133"].iloc[i] != "NA" and df["GENE_6480"].iloc[i] == "NA":
            
            id = update(df["GENE_4133"].iloc[i],id_updated)
          
            if id != None : 
                controle["new"] += 1
                id_updated.append(df["GENE_4133"].iloc[i])
            
            else:
                controle["none"] += 1
                id_updated.append(None)
                
        
        #-----------------------------------------------------------------------------------
        #cas ou les 2 gene_id ont une valeur, mais différente
        
        else: 
            id_4133 = update(df["GENE_4133"].iloc[i],id_updated)
            id_6180 = update(df["GENE_6480"].iloc[i],id_updated) 
            
            if id_4133 == id_6180:
                id_updated.append(id_4133)
                
            elif id_4133 and id_6180 is None:
                controle["old_maj"] += 1
                id_updated.append(id_4133)
                
                
            elif id_4133 is None and id_6180:
                controle["new_maj"] += 1
                id_updated.append(id_6180)
                
                
            else:
                controle["none"] += 1
                id_updated.append(None)
                

#ajouter la list dans une nouvelle colonne du tableau et exporter le nouveau dataframe  
df["UPDATED_GENE_ID"] = id_updated 
       
output_path = "C:/Users/noeme/OneDrive/Desktop/Stage_M1/Data/updated_merged_mapping_export.csv"
df.to_csv(output_path, index=False)
    
