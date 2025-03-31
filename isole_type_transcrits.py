import requests

def requete_data(symbole = ""):
    """
    Réalise pour une espèce et un symbole une recherche dans l'API Ensembl
    afin de récupérer les informations liées (Gene, Transcripts, etc.).
    """
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/homo_sapiens/{symbole}?expand=1"
    url = server + ext

    try:
        r = requests.get(url, headers={"Content-Type": "application/json"})
        r.raise_for_status() 
        return r.json()
    except requests.exceptions.RequestException as e:
        print(f"Erreur lors de la récupération de données : {e}")
        return None



def isole_data(symbole):
    """
    Isole les données depuis la requête API et les stocke dans un dictionnaire.
    """


    # Récupérer les données du gène via l'API
    dico_API = requete_data(symbole)
    
    biotype = dico_API.get("biotype","NA")
    
    return biotype

miRNA = "MIR21"
lncRNA = "NEAT1"

biotype = isole_data(lncRNA)
print(biotype)