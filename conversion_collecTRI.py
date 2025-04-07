import decoupler as dc
import pandas as pd 

#%%ouvre le dataset de collecTRI

df = dc.get_collectri(organism='human', split_complexes=False)
print(df.head())  # Affiche les premières lignes pour voir à quoi ça ressemble

#%% isole la liste des TF 

TF = list(df["source"])

#%%ouvre le fichier réaliser avec TFlink 

folder = "Data/Transcript/Output/"

DEG = pd.read_csv(folder + "DEG_complete_TFlink.csv", index_col=1)

DEG.drop(DEG.columns[0], axis=1, inplace=True)

#%%Parcour de nos valeurs 

i = 0
#parcour les lignes de la colonne type 
for Type in DEG["Type"]: 
    #regarde uniquement TF et protein coding 
    if Type == "TF" or Type == "protein_coding" : 
        
        #récupère le symbole 
        symbol = DEG.iloc[i,7]
        
        #si le symbol est dans la list de TF on change son type par TF
        if symbol in TF :
            DEG.iloc[i,8] = "TF"
         
        #sinon on change par protein coding 
        else : 
            DEG.iloc[i,8] = "protein_coding"
    i += 1
    
    
DEG.to_csv(folder + "DEG_complete_collecTRI.csv", index=1)


print("salut" if "ZBED1" in TF else "")
