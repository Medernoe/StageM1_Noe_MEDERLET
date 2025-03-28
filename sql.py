import mysql.connector
import pandas as pd 

#%%log sql 

mydb = mysql.connector.connect(
  host="localhost",
  user="medernoe",
  password="m23n25s22n24"
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
                    WHERE miR.nameR = "{miR_name}";""")
                    
    myresult = mycursor.fetchall()
    
    mycursor.close()

    return myresult




def gene_mirna_interaction(gene_id, miR_mimat):
    """
    Fais la requete pour avoir les information d'interactoion entre les un gene_id et un mimat 
    """
    mycursor = mydb.cursor()
    mycursor.execute("USE mirabel") 
    
    mycursor.execute(f"""SELECT predictions2.mimat, predictions2.id_g, predictions2.score, predictions2.VALID  
                    FROM predictions2 
                    JOIN genes3 ON genes3.ncbi_gene_id = predictions2.id_g 
                    WHERE genes3.ncbi_gene_id = "{gene_id}"
                    AND predictions2.mimat = "{miR_mimat}"
                    LIMIT 12;""")
                    
    myresult = mycursor.fetchall()
    
    mycursor.close()

    return myresult

inter = (gene_mirna_interaction("400451", "89") ) 
myresult = miR_name_mimat("hsa-let-7a")
#%% ouverture des CSV 
file_data_miRNA = "Data/DEG/Output/DEG_results_miRNA.csv"
file_data_gene = "Data/DEG/Output/DEG_results_19422.csv"

miRNA = pd.read_csv(file_data_miRNA, dtype=str, na_filter=False)
Gene = pd.read_csv(file_data_gene, dtype=str, na_filter=False)


#%% COnvertir les mirnames en mimat 

mimat = []
total = len(miRNA.iloc[:,0])
i = 0
for mirname in miRNA.iloc[:,0]: 
    i += 1 
    print(f"miRNA : {mirname}  -{i}/{total}")
    mimat = miR_name_mimat(mirname)
    if mimat : 
        mimat.append(mimat)
    else: 
        mimat.append("indisponible")
    
for mirname in miRNA.iloc[:,0]: 
    mimat = miR_name_mimat(mirname)
    mimat.append(mimat)


