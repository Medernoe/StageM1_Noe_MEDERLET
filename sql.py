import mysql.connector
import pandas as pd 
import json 


#%%log sql 

mydb = mysql.connector.connect(
  host="localhost",
  user="medernoe",
  password="1234", 
)

#%%fonction 

def miR_name_mimat(miR_name):
    """
    Fais la requete pour avoir un mimat 
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
    Fais la requete pour avoir les information d'interactoion entre les un gene_id et un mimat 
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


#%% ouverture des CSV 
file_data_miRNA = "Data/DEG/Output/DEG_results_miRNA.csv"
file_data_gene = "Data/Transcript/Output/DEG_complete_old.csv"

miRNA = pd.read_csv(file_data_miRNA, dtype=str, na_filter=False)
Gene = pd.read_csv(file_data_gene, dtype=str, na_filter=False)


#%% Convertir les mirnames en mimat 

list_mimat = []
total_miR = len(miRNA.iloc[:,0])
total_gene = len(Gene["Gene_id"]) 
i = 0
ii = 0
dico_interaction = {}
for mirname in miRNA.iloc[:,0]: 
    i += 1 
    
    mimat = miR_name_mimat(mirname)
    
    #transforme le tuples mimat en str unique 
    mimat = str(mimat[0][0]) 
    
    list_mimat.append(mimat)
    
    dico_interaction[mimat] = {}
    
    ii = 0
    
    for gene in Gene["Gene_id"] : 
        
        ii += 1 
        print(f"miRNA : {mirname}  -{i}/{total_miR} => gene {ii}/{total_gene}") 
        
        if gene not in dico_interaction[mimat] : 

            prediction = gene_mirna_interaction(gene, mimat)
            
            if prediction : 
                dico_interaction[mimat][gene] = {"score ": prediction[0][0] , 
                                             "val" : str(prediction[0][1]) } 
    
  
folder = "Data/Interaction/Output/"        
with open(folder + "interaction_miRNA_gene.json", "w") as f:
    json.dump(dico_interaction, f)            
            

