# Define directory variables for cleaner paths

# Répertoires principaux
data = "Data/"
output = "Output/"
input = "Input/"

# Sous-dossiers de traitement
pre = "Pre-traitement/"
deg = "DEG/"
conv = "Conversion/"
trans = "Transcript/"
inter = "Interaction/"
cyto = "Cytoscape/"
enri = "Enrichissement/"

# Dossiers spécifiques Input/Output pour chaque étape
pre_out    = data + pre + output
deg_out    = data + deg + output

conv_in    =  data + conv + input
conv_out   =  data + conv + output

trans_in   = data + trans + input
trans_out  = data + trans + output

inter_in   = data + inter + input
inter_out  = data + inter + output

cyto_in = data + cyto + input 

enri_in = data + enri + input
enri_out = data + enri + output

#base du répertoire du projet 
script_dir = workflow.basedir


rule all:
    input:
        # Pre-traitement outputs
        expand(pre_out + "{file}", file=[
            "gene_PCC.csv",
            "gene_PGL.csv",
            "gene_CONTROL.csv",
            "miRNA_PCC.csv",
            "miRNA_PGL.csv",
            "miRNA_CONTROL.csv",
        ]),

        # DEG outputs
        expand(deg_out + "{file}", file=[
            "DEG_results_gene.csv",
            "volcano_plot_gene.png",
            "DEG_results_miRNA.csv",
            "volcano_plot_miRNA.png",
        ]),

        # Interaction input from DEG_miRNA
        inter_in + "DEG_miRNA.csv",

        # Conversion outputs
        expand(conv_out + "{file}", file=[
            "update_gene_mapping.csv",
            "dico_id.json",
            "id_updated.json",
            "dico_controle.json",
        ]),

        # Transcript conversion & separation
        trans_in + "DEG_updated_gene.csv",
        trans_out + "DEG_final.csv",

        # Transcript separation outputs
        expand(inter_in + "{file}", file=[
            "DEG_TF.csv",
            "DEG_GENE.csv",
            "DEG_lncRNA.csv",
        ]),
        
        # Prediction output
        expand(inter_out + "{file}", file =[
            "interaction_miRNA_TF.json",
            "interaction_miRNA_GENE.json",
            "interaction_miRNA_lncRNA.json", 
            "interaction_TF_TF.json",
            "interaction_TF_GENE.json", 
            "interaction_TF_lncRNA.json",
            "interaction_TF_miRNA.json",
        ]),

        # Prediction output from DEG_miRNA
        inter_in + "DEG_miRNA_updated.csv",

        # Fichiers cytoscapes 
        cyto_in + "edges.csv", 
        cyto_in + "nodes.csv",

        expand(enri_in + "{file}", file =[
            "GO.csv",
            "GO_TF.csv",
            "GO_GENE.csv",
            "GO_miRNA.csv",
        ]),

        #fichier enrichissement
        expand(enri_out + "{file}", file =[
            "GO_all.png",
            "GO_down_all.png",
            "GO_down_genes.png",
            "GO_down_tfs.png",
            "GO_genes.png",
            "GO_mir.png",
            "GO_tfs.png",
            "GO_up_all.png",
            "GO_up_genes.png",
            "GO_up_tfs.png",
        ])

# 1) Pré-traitements
rule pre_traitement_gene:
    output:
        gene_PCC = pre_out + "gene_PCC.csv",
        gene_PGL = pre_out + "gene_PGL.csv",
        gene_control = pre_out + "gene_CONTROL.csv",
        merged_mapping_export = conv_in + "merged_mapping_export.csv",
    shell:
        """
        mkdir -p {pre_out} {conv_in}
        Rscript {script_dir}/pre_traitement_GENE.R
        """

rule pre_traitement_miRNA:
    output:
        miRNA_PCC = pre_out + "miRNA_PCC.csv",
        miRNA_PGL = pre_out + "miRNA_PGL.csv",
        miRNA_CONTROL = pre_out + "miRNA_CONTROL.csv",
    shell:
        """
        mkdir -p {pre_out}
        Rscript {script_dir}/pre_traitement_miRNA.R
        """

# 2) DEG
rule DEG_gene:
    input:
        gene_PCC = pre_out + "gene_PCC.csv",
        gene_control = pre_out + "gene_CONTROL.csv",
    output:
        DEG_results_gene = deg_out + "DEG_results_gene.csv",
        volcano_plot_gene = deg_out + "volcano_plot_gene.png",
    shell:
        """
        mkdir -p {deg_out}
        Rscript {script_dir}/DEG_gene.R
        """

rule DEG_miRNA:
    input:
        miRNA_CONTROL = pre_out + "miRNA_CONTROL.csv",
        miRNA_PCC = pre_out + "miRNA_PCC.csv",
    output:
        DEG_results_miRNA = deg_out + "DEG_results_miRNA.csv",
        volcano_plot_miRNA = deg_out + "volcano_plot_miRNA.png",
        DEG_miRNA = inter_in + "DEG_miRNA.csv",
    shell:
        """
        mkdir -p {deg_out} {inter_in}
        Rscript {script_dir}/DEG_miRNA.R
        """

# 3) Conversion d'IDs
rule conversion_oldID_to_newID:
    input:
        merged_mapping_export = conv_in + "merged_mapping_export.csv",
    output:
        update_gene_mapping = conv_out + "update_gene_mapping.csv",
        dico_id = conv_out + "dico_id.json",
        id_updated = conv_out + "id_updated.json",
        dico_controle = conv_out + "dico_controle.json",
    shell:
        """
        mkdir -p {conv_out}
        python3 {script_dir}/conversion_oldID_to_newID.py 
        """

rule conversion_spot_to_gene:
    input:
        update_gene_mapping = conv_out + "update_gene_mapping.csv",
        dico_id = conv_out + "dico_id.json",
        DEG_results_gene = deg_out + "DEG_results_gene.csv",
    output:
        DEG_updated_gene = trans_in + "DEG_updated_gene.csv",
    shell:
        """
        mkdir -p {trans_in}
        Rscript {script_dir}/conversion_SpotToGene.R
        """

# 4) Types de transcript et séparation
rule transcript_types:
    input:
        DEG_updated_gene = trans_in + "DEG_updated_gene.csv",
    output:
        DEG_final = trans_out + "DEG_final.csv",
    shell:
        """
        mkdir -p {trans_out}
        python3 {script_dir}/transcript_types.py
        """

rule transcript_separate:
    input:
        DEG_final = trans_out + "DEG_final.csv",
    output:
        DEG_TF = inter_in + "DEG_TF.csv",
        DEG_GENE = inter_in + "DEG_GENE.csv",
        DEG_lncRNA = inter_in + "DEG_lncRNA.csv",
    shell:
        """
        mkdir -p {inter_in}
        Rscript {script_dir}/transcript_DEG_separateur.R
        """

# 5) Prediction des interactions 
rule Interaction_miRNA:
    input:
        DEG_TF = inter_in + "DEG_TF.csv",
        DEG_GENE = inter_in + "DEG_GENE.csv",
        DEG_lncRNA = inter_in + "DEG_lncRNA.csv",
        DEG_miRNA = inter_in + "DEG_miRNA.csv",
    output:
        DEG_miRNA_updated = inter_in + "DEG_miRNA_updated.csv",
        interaction_miRNA_TF = inter_out + "interaction_miRNA_TF.json",
        interaction_miRNA_GENE = inter_out + "interaction_miRNA_GENE.json",
        interaction_miRNA_lncRNA = inter_out + "interaction_miRNA_lncRNA.json",
    shell:
        """
        mkdir -p {inter_out}
        python3 {script_dir}/interaction_miRNA.py
        """

rule Interaction_TF:
    input:
        DEG_miRNA_updated = inter_in + "DEG_miRNA_updated.csv",
        DEG_TF = inter_in + "DEG_TF.csv",
        DEG_GENE = inter_in + "DEG_GENE.csv",
        DEG_lncRNA = inter_in + "DEG_lncRNA.csv",
    output:
        interaction_TF_TF = inter_out + "interaction_TF_TF.json",
        interaction_TF_GENE = inter_out + "interaction_TF_GENE.json",
        interaction_TF_lncRNA = inter_out + "interaction_TF_lncRNA.json",
        interaction_TF_miRNA = inter_out + "interaction_TF_miRNA.json",
    shell:
        """
        mkdir -p {inter_out}
        python3 {script_dir}/interaction_TF.py
        """

# 6) Constructions des fichiers cytoscape 
rule cytoscape: 
    input: 
        DEG_miRNA_updated = inter_in + "DEG_miRNA_updated.csv",
        DEG_TF = inter_in + "DEG_TF.csv",
        DEG_GENE = inter_in + "DEG_GENE.csv",
        dico_id = conv_out + "dico_id.json",

        interaction_miRNA_TF = inter_out + "interaction_miRNA_TF.json",
        interaction_miRNA_GENE = inter_out + "interaction_miRNA_GENE.json",
        interaction_miRNA_lncRNA = inter_out + "interaction_miRNA_lncRNA.json",

        interaction_TF_TF = inter_out + "interaction_TF_TF.json",
        interaction_TF_GENE = inter_out + "interaction_TF_GENE.json",
        interaction_TF_lncRNA = inter_out + "interaction_TF_lncRNA.json",
        interaction_TF_miRNA = inter_out + "interaction_TF_miRNA.json",

    output: 
        nodes = cyto_in + "nodes.csv",
        edges = cyto_in + "edges.csv",

        GO = enri_in + "GO.csv",
        GO_TF = enri_in + "GO_TF.csv",
        GO_GENE = enri_in + "GO_GENE.csv",
        GO_miRNA = enri_in + "GO_miRNA.csv",    
    shell: 
        """
        mkdir -p {cyto_in} {enri_in}
        python3 {script_dir}/cytoscape_FileBuilder.py

        """

# 7) Enrichissement fonctionnel 
rule enrichissement: 
    input: 
        GO = enri_in + "GO.csv",
        GO_TF = enri_in + "GO_TF.csv",
        GO_GENE = enri_in + "GO_GENE.csv",
        GO_miRNA = enri_in + "GO_miRNA.csv",   
    output: 
        GO_all = enri_out + "GO_all.png", 
        GO_down_all = enri_out + "GO_down_all.png", 
        GO_down_genes = enri_out + "GO_down_genes.png", 
        GO_down_tfs = enri_out + "GO_down_tfs.png", 
        GO_genes = enri_out + "GO_genes.png", 
        GO_mir = enri_out + "GO_mir.png", 
        GO_tfs = enri_out + "GO_tfs.png", 
        GO_up_all = enri_out + "GO_up_all.png", 
        GO_up_genes = enri_out + "GO_up_genes.png", 
        GO_up_tfs = enri_out + "GO_up_tfs.png", 
    shell: 
        """
        mkdir -p {enri_out}         
        Rscript {script_dir}/Enrichissement.R
        """
