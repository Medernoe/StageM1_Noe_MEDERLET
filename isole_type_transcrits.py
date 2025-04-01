import requests
import pandas as pd 

def requete_data(symbole = ""):
    """
    Réalise une requete NCBI à partir d'un symbole et sort un dico 
    """
    
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/homo_sapiens/{symbole}?expand=1"
    url = server + ext



    attempts = 0 #mise à jour de la variable attempts 
    max_attempts = 5
    dico = None 
    
    #réalise max 5 fois la requete API 
    while attempts < max_attempts:
        #continue les requete tant que le dico est vide 
        if dico == None : 
            attempts += 1 
            
            #essaye la requete et si cela fontionne sort le dico 
            try:
                r = requests.get(url, headers={"Content-Type": "application/json"})
                r.raise_for_status() 
                dico = r.json()
            
            #sinon sort un dico vide 
            except requests.exceptions.RequestException as e:
                print(f"Erreur lors de la récupération de données : {e}")
                dico = None
        
        #stop les requete si une requete à fonctionnée 
        else : 
            attempts = max_attempts
            
    return dico 





def isole_data(symbole):
    """
    Isole le biotype, possibilité de récupérer d'autre informations 
    """

    # Récupérer les données du gène via l'API
    dico_API = requete_data(symbole)
    
    if dico_API : 
        #isole le type (mRNA, miRNA, lncRNA, ...)
        biotype = dico_API.get("biotype",None)
        #ensembl = dico_API.get("id",None)
        
    else : 
        #mets unknown si aucun dico n'est obtenue 
        biotype = "Unknown"
    
    
    return biotype 



#test à la main des fonctions 
"""
miRNA = "MIR21"
lncRNA = "NEAT1"

biotype = isole_data(lncRNA)
"""





#%%
#répertoire de travail 
folder = "Data/Transcript/Input/"
folder_output = "Data/Transcript/Output/"
#fichier de DEG 
gene  = pd.read_csv(folder + "DEG_updated_gene.csv")


#%% ouvrir et lire les facteur de transcrption 

# Lire le fichier 
with open(folder + "Homo_sapiens_TF.txt", "r", encoding="cp1252") as file:
    lines = file.readlines()

# Extraire la colonne "Symbol" 
symbol_list = [line.split("\t")[1].strip() for line in lines[1:] if line.split("\t")[1].strip()]

# Extraire la colonne "ensembl_id" 
#ensembl_id = [line.split("\t")[2].strip() for line in lines[1:] if line.split("\t")[1].strip()]

#%%creer la colonne de biotype 

#initialisation des variables 
biotype = []
i = 0
total = len(gene["Symbol"])

for symbole in gene["Symbol"]:  
    print(f"Symbole {i}/{total}")
    
    if symbole in symbol_list : 
        biotype.append("TF")
        
    else: 
        biotype.append(isole_data(symbole)) 
    
    print("=>" + biotype[i])
    i += 1
    
gene["Type"] = biotype


gene.to_csv(folder_output + "DEG_complete.csv", index=False)