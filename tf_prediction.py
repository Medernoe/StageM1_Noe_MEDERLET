import decoupler as dc
import pandas as pd
import sys  

#%% 1. Chargement du dataset CollecTRI

# On récupère le dataset via la fonction get_collectri du module decoupler
df = dc.get_collectri(organism='human', split_complexes=False)


#%% 2. Ouverture des fichiers CSV

# Fichier contenant les informations des TF
file_data_TF = "Data/Interaction/Input/DEG_TF.csv"

# Chargement du CSV 
TF = pd.read_csv(file_data_TF, dtype=str, na_filter=False, index_col=1)
# Suppression de la première colonne 
TF.drop(TF.columns[0], axis=1, inplace=True)


# Lecture temporaire à la main du fichier de données avec les gènes
#data = pd.read_csv("Data/Interaction/Input/DEG_GENE.csv", dtype=str, na_filter=False, index_col=1)
#data.drop(data.columns[0], axis=1, inplace=True)


file_data = sys.argv[1]
data = pd.read_csv(file_data, dtype=str, na_filter=False)
data.drop(data.columns[0], axis=1, inplace=True)
#%% 3. Prédiction des interactions

# Récupère la liste de nos TF 
# La conversion en set permet de retirer les doublons, puis on retransforme en liste
list_TF = list(set(TF.iloc[:, -2]))
total_TF = len(list_TF)  # Nombre total de TF uniques

# Initialisation des compteurs et du dictionnaire qui contiendra les interactions
i = 0  # Compteur pour suivre l'avancement sur les TF
dico_interaction = {}

# Prépare un ensemble (set) des symboles présents dans le fichier de data 
symbol_set = set(data["Symbol"])

# Parcourt chaque TF dans notre  de TF
for tf in list_TF:
    i += 1
    
    # Initialise un dico pour chaque TF 
    dico_interaction[tf] = {}
    
    # Filtre le df de colecTRI pour obtenir les lignes correspondante au tf de notre jeux de données 
    subset = df[df["source"] == tf][["target", "weight"]]
    
    
    total_cible = len(subset)  # Nombre total de cibles potentielles pour ce TF
    ii = 0  # Compteur pour suivre l'avancement des cibles pour le TF actuel
    
    # Parcourt chaque ligne du sous-ensemble filtré (cible potentiel et sens)
    for index, row in subset.iterrows():
        gene = row["target"]    # Le gène cible
        sens = row["weight"]    # Le sens de l'interaction
        ii += 1

        # Affiche l'avancement de l'analyse pour ce TF
        print(f"TF : {tf}  - {i}/{total_TF} => cible potentielle {ii}/{total_cible}")
        
        # Vérifie si le gène ciblé est présent dans notre jeu de données data (symbol_set)
        if gene in symbol_set:
            # Ajoute l'interaction dans le dictionnaire sous forme de dictionnaire (clé: target, valeur: sens)
            dico_interaction[tf][gene] = sens
    
    #supprime le dico de la list si il n'a pas d'interaction
    if not dico_interaction[tf]:
        del dico_interaction[tf]

#%% 4. Sauvegarde des interactions 

import json

# Utilise un argument de la ligne de commande pour définir le nom de fichier de sortie
name = sys.argv[2]
full_name = f"interaction_TF_{name}.json"
folder = "Data/Interaction/Output/"

with open(folder + full_name, "w") as f:
    json.dump(dico_interaction, f)

print(f"Prédiction terminée, fichier enregister ici : {folder + full_name}")