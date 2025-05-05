# Define directory variables for cleaner paths
pre_out    = "Data/Pre-traitement/Output/"
deg_out    = "Data/DEG/Output/"
conv_in    = "Data/Conversion/Input/"
conv_out   = "Data/Conversion/Output/"
trans_in   = "Data/Transcript/Input/"
trans_out  = "Data/Transcript/Output/"
inter_in   = "Data/Interaction/Input/"
inter_out  = "Data/Interaction/Output/"
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

        

# 1) Pré-traitements
rule pre_traitement_gene:
    output:
        pre_out + "gene_PCC.csv",
        pre_out + "gene_PGL.csv",
        pre_out + "gene_CONTROL.csv",
        conv_in + "merged_mapping_export.csv",
    shell:
        """
        mkdir -p {pre_out} {conv_in}
        Rscript {script_dir}/pre_traitement_GENE.R
        """

rule pre_traitement_miRNA:
    output:
        pre_out + "miRNA_PCC.csv",
        pre_out + "miRNA_PGL.csv",
        pre_out + "miRNA_CONTROL.csv",
    shell:
        """
        mkdir -p {pre_out}
        Rscript {script_dir}/pre_traitement_miRNA.R
        """

# 2) DEG
rule DEG_gene:
    input:
        pre_out + "gene_PCC.csv",
        pre_out + "gene_CONTROL.csv",
    output:
        deg_out + "DEG_results_gene.csv",
        deg_out + "volcano_plot_gene.png",
    shell:
        """
        mkdir -p {deg_out}
        Rscript {script_dir}/DEG_gene.R
        """

rule DEG_miRNA:
    input:
        miRNA_CONTROL=pre_out + "miRNA_CONTROL.csv",
        miRNA_PCC=pre_out + "miRNA_PCC.csv",
    output:
        deg_out + "DEG_results_miRNA.csv",
        deg_out + "volcano_plot_miRNA.png",
        inter_in + "DEG_miRNA.csv",
    shell:
        """
        mkdir -p {deg_out} {inter_in}
        Rscript {script_dir}/DEG_miRNA.R
        """

# 3) Conversion d'IDs
rule conversion_oldID_to_newID:
    input:
        conv_in + "merged_mapping_export.csv",
    output:
        conv_out + "update_gene_mapping.csv",
        conv_out + "dico_id.json",
        conv_out + "id_updated.json",
        conv_out + "dico_controle.json",
    shell:
        """
        mkdir -p {conv_out}
        python3 {script_dir}/conversion_oldID_to_newID.py
        """

rule conversion_spot_to_gene:
    input:
        update_mapping=conv_out + "update_gene_mapping.csv",
        dico_id=conv_out + "dico_id.json",
        DEG_results_gene=deg_out + "DEG_results_gene.csv",
    output:
        trans_in + "DEG_updated_gene.csv",
    shell:
        """
        mkdir -p {trans_in}
        Rscript {script_dir}/conversion_SpotToGene.R
        """

# 4) Types de transcript et séparation
rule transcript_types:
    input:
        trans_in + "DEG_updated_gene.csv",
    output:
        trans_out + "DEG_final.csv",
    shell:
        """
        mkdir -p {trans_out}
        python3 {script_dir}/transcript_types.py
        """

rule transcript_separate:
    input:
        trans_out + "DEG_final.csv",
    output:
        inter_in + "DEG_TF.csv",
        inter_in + "DEG_GENE.csv",
        inter_in + "DEG_lncRNA.csv",
    shell:
        """
        mkdir -p {inter_in}
        Rscript {script_dir}/transcript_DEG_separateur.R
        """

# 5) Prediction des interactions 
rule Interaction_miRNA:
    input:
        inter_in + "DEG_miRNA.csv",
        inter_in + "DEG_TF.csv",
        inter_in + "DEG_GENE.csv",
        inter_in + "DEG_lncRNA.csv",
    output:
        inter_in + "DEG_miRNA_updated.csv",
        inter_out + "interaction_miRNA_TF.json",
        inter_out + "interaction_miRNA_GENE.json",
        inter_out + "interaction_miRNA_lncRNA.json",
    shell:
        """
        mkdir -p {inter_out}
        python3 {script_dir}/interaction_miRNA.py
        """

rule Interaction_TF:
    input:
        inter_in + "DEG_TF.csv",
        inter_in + "DEG_GENE.csv",
        inter_in + "DEG_lncRNA.csv",
    output:
        inter_out + "interaction_TF_TF.json",
        inter_out + "interaction_TF_GENE.json",
        inter_out + "interaction_TF_lncRNA.json",
        inter_out + "interaction_TF_miRNA.json",
    shell:
        """
        mkdir -p {inter_out}
        python3 {script_dir}/interaction_TF.py
        """

