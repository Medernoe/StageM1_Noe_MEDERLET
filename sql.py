import mysql.connector
import pandas as pd 
import json 
import sys 


#%%log sql 

mydb = mysql.connector.connect(
  host="localhost",
  user="medernoe",
  password="1234", 
)

#%%fonction 

def miR_name_mimat(miR_name):
    """
    Fais la requete pour avoir un mimat à partir d'un mirNAME
    """
    mycursor = mydb.cursor()
    mycursor.execute("USE mirabel") 
    
    mycursor.execute(f"""SELECT miR_mimat.number, mimat.number, predictions2.mimat
                    FROM miR
                    JOIN miR_mimat ON miR.id = miR_mimat.mir_id
                    JOIN mimat ON miR_mimat.number = mimat.number
                    JOIN predictions2 ON mimat.number = predictions2.mimat
                    WHERE miR.nameR = "{miR_name}"
                    LIMIT 1;""")
                    
    myresult = mycursor.fetchall()
    
    mycursor.close()

    return myresult




def gene_mirna_interaction(gene_id, miR_mimat):
    """
    Fais la requete pour avoir les information d'interactoion entre un gene_id et un mimat 
    """
    mycursor = mydb.cursor()
    mycursor.execute("USE mirabel") 
    
    mycursor.execute(f"""SELECT predictions2.score, predictions2.VALID  
                    FROM predictions2 
                    JOIN genes3 ON genes3.ncbi_gene_id = predictions2.id_g 
                    WHERE genes3.ncbi_gene_id = "{gene_id}"
                    AND predictions2.mimat = "{miR_mimat}"
                    LIMIT 12;""")
                    
    myresult = mycursor.fetchall()
    
    mycursor.close()

    return myresult



def pourcent_trie(dico, pourcentage = 0.05): 
    """
    Retourne les 5% des meilleurs scrores du dico 
    """
    dico_trie = dict(sorted(dico.items(), key=lambda item: item[1]['score'], reverse=False)) 
    
    n = len(dico_trie)
    top = int(n * pourcentage) or 1
     
    top_dico = dict(list(dico_trie.items())[:top])  
    
    return top_dico

#%% ouverture des CSV 
file_data_miRNA = "Data/DEG/Output/DEG_results_miRNA.csv"
file_data = sys.argv[1]
 
    
miRNA = pd.read_csv(file_data_miRNA, dtype=str, na_filter=False)
miRNA.to_csv("Data/Interaction/Input/DEG_miRNA.csv", index=1)


data = pd.read_csv(file_data, dtype=str, na_filter=False)


#%% Convertir les mirnames en mimat 

list_mimat = []
total_miR = len(miRNA.iloc[:,0])
total = len(data["Gene_id"]) 
i = 0
ii = 0
dico_interaction = {}
for mirname in miRNA.iloc[:,0]: 
    i += 1 
    
    mimat = miR_name_mimat(mirname)
    
    #transforme le tuples mimat en str unique 
    mimat = str(mimat[0][0]) 
    
    list_mimat.append(mimat)
    
    dico_interaction[mirname] = {}
    
    ii = 0
    
    for gene in data["Gene_id"] : 
        
        ii += 1 
        print(f"miRNA : {mirname}  -{i}/{total_miR} => gene {ii}/{total}") 
        
        if gene not in dico_interaction[mirname] : 

            prediction = gene_mirna_interaction(gene, mimat)
            
            if prediction : 
                dico_interaction[mirname][gene] = {"score": prediction[0][0] , 
                                             "val" : str(prediction[0][1]) } 
                
    dico_interaction[mirname] = pourcent_trie(dico_interaction[mirname])

name = sys.argv[2]
full_name = f"interaction_miRNA_{name}.json"
folder = "Data/Interaction/Output/"        
with open(folder + full_name, "w") as f:
    json.dump(dico_interaction, f)            
                     
            