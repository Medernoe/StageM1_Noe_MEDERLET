from Bio import Entrez

import pandas as pd
#%% Import de la table de jointure 


# Vous pouvez aussi préciser le type dès la lecture du CSV pour certaines colonnes via le paramètre dtype
df = pd.read_csv("C:/Users/noeme/OneDrive/Desktop/Stage_M1/Data/merged_mapping_export.csv", dtype=str, na_filter=False)


# Conversion des colonnes en string
gene_4133 = df["GENE_4133"]
gene_6480 = df.iloc[:, 8]

nouveaux = []
coord = []
old = []
new = []
coord_double = []
for i in range(len(gene_4133)): 
    if gene_4133.iloc[i] == gene_6480.iloc[i]:
        nouveaux.append(gene_4133.iloc[i])
        
    elif gene_4133.iloc[i] != gene_6480.iloc[i]:
        
        coord.append(i)
        
        print(gene_4133.iloc[i], "//", gene_6480.iloc[i])
        
        if gene_4133.iloc[i] == "NA" and gene_6480.iloc[i] != "NA": 
            old.append(gene_6480.iloc[i])
          
        elif gene_4133.iloc[i] != "NA" and gene_6480.iloc[i] == "NA":
            new.append(gene_4133.iloc[i])
        
        else: 
            coord_double.append(i)
            












#%%

def update(gene_id):
    Entrez.email = "noe.mederlet@univ-rouen.fr"
    
    try:
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print("Erreur lors de la récupération :", e)
        records = []
    
    if not records:
        print(f"Aucune correspondance pour l'ID {gene_id}.")
    else:
        gene_data = records[0]
        
        # Récupération des informations de suivi du gène
        track_info = gene_data.get("Entrezgene_track-info", {})
        gene_track = track_info.get("Gene-track", {})
        
        # Statut du gene : 'current' ou 'secondary' (obsolète / remplacé)
        status = gene_track.get("Gene-track_status", None)
        # Dans le dictionnaire Biopython, status est souvent un StringElement 
        # qui peut avoir un attribut 'value'.
        status_value = getattr(status, "attributes", {}).get("value", None)
        
        if status_value == "secondary":
            # L'ID a été remplacé ou discontinué
            current_ids = gene_track.get("Gene-track_current-id", [])
            
            # On cherche un GeneID de remplacement
            new_gene_id = None
            for dbtag in current_ids:
                if dbtag.get("Dbtag_db") == "GeneID":
                    # Récupération du nouvel ID
                    new_gene_id = dbtag["Dbtag_tag"]["Object-id"]["Object-id_id"]
                    break
            
            if new_gene_id:
                print(f"Le GeneID {gene_id} est obsolète et a été remplacé par {new_gene_id}.")
                
            else:
                print(f"Le GeneID {gene_id} est obsolète, et aucun nouvel ID n'a été trouvé.")
        else:
            # Si le statut n’est pas “secondary”, on considère que l’ID est valide
            print(f"Le GeneID {gene_id} est valide.")
    
            # Vous pouvez extraire et afficher plus de détails, par exemple :
            gene_ref = gene_data.get("Entrezgene_gene", {}).get("Gene-ref", {})
            locus = gene_ref.get("Gene-ref_locus", "N/A")
            desc = gene_ref.get("Gene-ref_desc", "N/A")
            print(f"Locus       : {locus}")
            print(f"Description : {desc}")
