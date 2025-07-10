# -*- Snakefile -*-

# Pour pouvoir passer data_dir / input_dir / output_dir via --config
DATA_DIR   = config.get("data_dir",   "Data3/")
INPUT_DIR  = config.get("input_dir",  "Input/")
OUTPUT_DIR = config.get("output_dir", "Output/")


# Identifiant connexion à mySQL
USER = "medernoe"
PASSWORD = "1234"


# Sous-dossiers de traitement
PRE  = "Pre-traitement/"
DEG  = "DEG/"
CONV = "Conversion/"
TRANS= "Transcript/"
INTER= "Interaction/"
CYTO = "Cytoscape/"
ENRI = "Enrichissement/"

# Chemins construits
pre_out   = DATA_DIR + PRE   + OUTPUT_DIR
deg_out   = DATA_DIR + DEG   + OUTPUT_DIR
conv_in   = DATA_DIR + CONV  + INPUT_DIR
conv_out  = DATA_DIR + CONV  + OUTPUT_DIR
trans_in  = DATA_DIR + TRANS + INPUT_DIR
trans_out = DATA_DIR + TRANS + OUTPUT_DIR
inter_in  = DATA_DIR + INTER + INPUT_DIR
inter_out = DATA_DIR + INTER + OUTPUT_DIR
cyto_in   = DATA_DIR + CYTO  + INPUT_DIR
enri_in   = DATA_DIR + ENRI  + INPUT_DIR
enri_out  = DATA_DIR + ENRI  + OUTPUT_DIR

# Base du projet
script_dir = workflow.basedir

# ——————————————————————————————————————————————————————————————
# Liste des comparaisons à faire basés sur les séparations dans pre-traitement
COMPARISONS = {
    "PCC_Sains":   ("CONTROL", "PCC"),
    "PGL_Sains":   ("CONTROL", "PGL"),
    "PPGL_Sains":  ("CONTROL", "PPGL"),
    "PCC_PGL":     ("PGL",     "PCC"),
}
COMPS = list(COMPARISONS.keys())

# Règle finale : on génère tout pour chaque comparaison
rule all:
    input:
        # pré-traitement généraux
        expand(pre_out + "{file}", file=[
            "gene_PCC.csv","gene_PGL.csv","gene_CONTROL.csv","gene_PPGL.csv",
            "miRNA_PCC.csv","miRNA_PGL.csv","miRNA_CONTROL.csv","miRNA_PPGL.csv",
        ]),
        # DEG pour chaque comparaison
        expand(deg_out + "{comp}/DEG_results_gene.csv",   comp=COMPS),
        expand(deg_out + "{comp}/volcano_plot_gene.png",  comp=COMPS),
        expand(deg_out + "{comp}/DEG_results_miRNA.csv",  comp=COMPS),
        expand(deg_out + "{comp}/volcano_plot_miRNA.png", comp=COMPS),

        # Conversion outputs
        expand(conv_out + "{file}", file=[
            "update_gene_mapping.csv",
            "dico_id.json",
            "id_updated.json",
            "dico_controle.json",
        ]),

        # Mise à jour des résultats DEG pour chaque comparaison
        expand(trans_in + "{comp}/DEG_updated_gene.csv",   comp=COMPS),

        # Types de transcript pour chaque comparaison
        expand(trans_out + "{comp}/DEG_final_gene.csv",   comp=COMPS),

        # Séparation des types de transcript pour chaque comparaison
        expand(inter_in + "{comp}/DEG_TF.csv",      comp=COMPS),
        expand(inter_in + "{comp}/DEG_GENE.csv",    comp=COMPS),
        expand(inter_in + "{comp}/DEG_lncRNA.csv",  comp=COMPS),    

        # Mise à jour des miRNA
        expand(inter_in + "{comp}/DEG_miRNA.csv", comp=COMPS),

        # Prédiction des interactions pour chaque comparaison
        #miRNA
        expand(inter_out + "{comp}/interaction_miRNA_TF.json",     comp=COMPS),
        expand(inter_out + "{comp}/interaction_miRNA_GENE.json",     comp=COMPS),
        expand(inter_out + "{comp}/interaction_miRNA_lncRNA.json",     comp=COMPS),
        #TF
        expand(inter_out + "{comp}/interaction_TF_TF.json",     comp=COMPS),
        expand(inter_out + "{comp}/interaction_TF_GENE.json",     comp=COMPS),
        expand(inter_out + "{comp}/interaction_TF_lncRNA.json",     comp=COMPS),
        expand(inter_out + "{comp}/interaction_TF_miRNA.json",     comp=COMPS),

        # Fichiers Cytoscape pour chaque comparaison
        expand(cyto_in + "{comp}/nodes.csv", comp=COMPS),
        expand(cyto_in + "{comp}/edges.csv", comp=COMPS),
        expand(enri_in + "{comp}/GO.csv", comp=COMPS),
        expand(enri_in + "{comp}/GO_TF.csv", comp=COMPS),
        expand(enri_in + "{comp}/GO_GENE.csv", comp=COMPS),
        expand(enri_in + "{comp}/GO_miRNA.csv", comp=COMPS),   

        # Enrichissement fonctionnel pour chaque comparaison
        # up/down
        expand(enri_out + "{comp}/GO_all.png", comp=COMPS),
        expand(enri_out + "{comp}/GO_genes.png", comp=COMPS),
        expand(enri_out + "{comp}/GO_tfs.png", comp=COMPS),
        # down
        expand(enri_out + "{comp}/GO_down_all.png", comp=COMPS),
        expand(enri_out + "{comp}/GO_down_genes.png", comp=COMPS),
        expand(enri_out + "{comp}/GO_down_tfs.png", comp=COMPS),
        # up
        expand(enri_out + "{comp}/GO_up_all.png", comp=COMPS),
        expand(enri_out + "{comp}/GO_up_genes.png", comp=COMPS),
        expand(enri_out + "{comp}/GO_up_tfs.png", comp=COMPS) 
# ——————————————————————————————————————————————————————————————
# 1) Pré-traitements
# Pré-traitement des gènes
rule pre_traitement_gene:
    output:
        gene_PCC     = pre_out + "gene_PCC.csv",
        gene_PGL     = pre_out + "gene_PGL.csv",
        gene_CONTROL = pre_out + "gene_CONTROL.csv",
        gene_PPGL    = pre_out + "gene_PPGL.csv",
        merged_map   = conv_in + "merged_mapping_export.csv",
    params:
        pre  = pre_out,
        conv = conv_in,
    shell:
        """
        mkdir -p {params.pre} {params.conv}
        Rscript {script_dir}/pre_traitement_GENE.R \
            {output.gene_PCC} \
            {output.gene_PGL} \
            {output.gene_CONTROL} \
            {output.gene_PPGL} \
            {output.merged_map} \
            {params.pre} \
            {params.conv}
        """

# Pré-traitement des miRNA
rule pre_traitement_miRNA:
    output:
        miRNA_PCC     = pre_out + "miRNA_PCC.csv",
        miRNA_PGL     = pre_out + "miRNA_PGL.csv",
        miRNA_CONTROL = pre_out + "miRNA_CONTROL.csv",
        miRNA_PPGL    = pre_out + "miRNA_PPGL.csv",
    params:
        pre  = pre_out,
    shell:
        """
        mkdir -p {params.pre}
        Rscript {script_dir}/pre_traitement_miRNA.R \
            {output.miRNA_PCC} \
            {output.miRNA_PGL} \
            {output.miRNA_CONTROL} \
            {output.miRNA_PPGL} \
            {params.pre} \
        """

# ——————————————————————————————————————————————————————————————
# 2) DEG (Différentiel d'Expression Génétique)
# DEG — gènes
rule DEG_gene:
    params:
        comp="{comp}"
    input:
        control   = lambda wc: pre_out + f"gene_{COMPARISONS[wc.comp][0]}.csv",
        tumor = lambda wc: pre_out + f"gene_{COMPARISONS[wc.comp][1]}.csv",
    output:
        results = deg_out + "{comp}/DEG_results_gene.csv",
        plot    = deg_out + "{comp}/volcano_plot_gene.png",
    shell:
        """
        mkdir -p {deg_out}{params.comp}
        Rscript {script_dir}/DEG_gene.R \
            {input.control} \
            {input.tumor} \
            {output.results} \
            {output.plot}
        """

# DEG — miRNA
rule DEG_miRNA:
    params:
        comp="{comp}"
    input:
        control   = lambda wc: pre_out + f"miRNA_{COMPARISONS[wc.comp][0]}.csv",
        tumor = lambda wc: pre_out + f"miRNA_{COMPARISONS[wc.comp][1]}.csv",
    output:
        results = deg_out + "{comp}/DEG_results_miRNA.csv",
        plot    = deg_out + "{comp}/volcano_plot_miRNA.png",
    shell:
        """
        mkdir -p {deg_out}{params.comp}
        Rscript {script_dir}/DEG_miRNA.R \
            {input.control} \
            {input.tumor} \
            {output.results} \
            {output.plot}
        """

# ——————————————————————————————————————————————————————————————
# 3) Conversion d'IDs

rule conversion_oldID_to_newID:
    input:
        merged_mapping_export = conv_in + "merged_mapping_export.csv",
    output:
        update_gene_mapping = conv_out + "update_gene_mapping.csv",
        dico_id            = conv_out + "dico_id.json",
        id_updated         = conv_out + "id_updated.json",
        dico_controle      = conv_out + "dico_controle.json",
    params:
        conv = conv_out,
    shell:
        """
        mkdir -p {params.conv}
        python3 conversion_oldID_to_newID.py \
           {input.merged_mapping_export} \
           {output.update_gene_mapping} \
           {output.dico_id} \
           {output.id_updated} \
           {output.dico_controle} \
           {params.conv} \
        """

rule conversion_spot_to_gene:
    params:
        comp = "{comp}",
        trans = trans_in,
    input:
        update_gene_mapping = conv_out + "update_gene_mapping.csv",
        dico_id = conv_out + "dico_id.json",
        results = deg_out + "{comp}/DEG_results_gene.csv",
    output:
        DEG_updated_gene = trans_in + "{comp}/DEG_updated_gene.csv",

    shell:
        """
        mkdir -p {params.trans}{params.comp}
        Rscript {script_dir}/conversion_SpotToGene.R \
            {input.update_gene_mapping} \
            {input.dico_id} \
            {input.results} \
            {output.DEG_updated_gene} \
            {params.trans}{params.comp}
        """

# ——————————————————————————————————————————————————————————————
# 4) Types de transcript et séparation

rule transcript_types:
    params:
        comp = "{comp}",
        trans = trans_out,
    input:
        DEG_updated_gene = trans_in + "{comp}/DEG_updated_gene.csv",
    output: 
        DEG_final = trans_out + "{comp}/DEG_final_gene.csv",
    shell:
        """
        mkdir -p {params.trans}{params.comp}
        python3 {script_dir}/transcript_types.py \
            {input.DEG_updated_gene} \
            {output.DEG_final} \
            {params.trans}{params.comp}
        """

rule transcript_separate:
    params:
        comp = "{comp}",
        inter = inter_in,
    input:
        DEG_final = trans_out + "{comp}/DEG_final_gene.csv",
    output:
        DEG_GENE = inter_in + "{comp}/DEG_GENE.csv",
        DEG_TF = inter_in + "{comp}/DEG_TF.csv",
        DEG_lncRNA = inter_in + "{comp}/DEG_lncRNA.csv",
    shell:
        """
        mkdir -p {params.inter}{params.comp}
        Rscript {script_dir}/transcript_DEG_separateur.R \
            {input.DEG_final} \
            {output.DEG_GENE} \
            {output.DEG_TF} \
            {output.DEG_lncRNA} \
            {params.inter}{params.comp}
        """

# ——————————————————————————————————————————————————————————————
# 5) Prediction des interactions 
rule Interaction_miRNA:
    params:
        comp = "{comp}",
        inter = inter_out,
        user = USER,
        password = PASSWORD,
    input:
        DEG_TF = inter_in + "{comp}/DEG_TF.csv",
        DEG_GENE = inter_in + "{comp}/DEG_GENE.csv",
        DEG_lncRNA = inter_in + "{comp}/DEG_lncRNA.csv",
        DEG_miRNA = deg_out + "{comp}/DEG_results_miRNA.csv",
    output:
        DEG_miRNA_updated = inter_in + "{comp}/DEG_miRNA.csv",
        interaction_miRNA_TF = inter_out + "{comp}/interaction_miRNA_TF.json",
        interaction_miRNA_GENE = inter_out + "{comp}/interaction_miRNA_GENE.json",
        interaction_miRNA_lncRNA = inter_out + "{comp}/interaction_miRNA_lncRNA.json",
    shell:
        """
        mkdir -p {params.inter}{params.comp}
        python3 {script_dir}/interaction_miRNA.py \
            {input.DEG_miRNA} \
            {input.DEG_TF} \
            {input.DEG_GENE} \
            {input.DEG_lncRNA} \
            {output.DEG_miRNA_updated} \
            {output.interaction_miRNA_TF} \
            {output.interaction_miRNA_GENE} \
            {output.interaction_miRNA_lncRNA} \
            {params.inter}{params.comp} \
            {params.user} \
            {params.password}
        """

rule Interaction_TF:
    params:
        comp = "{comp}",
        inter = inter_out,
    input:
        DEG_TF = inter_in + "{comp}/DEG_TF.csv",
        DEG_GENE = inter_in + "{comp}/DEG_GENE.csv",
        DEG_lncRNA = inter_in + "{comp}/DEG_lncRNA.csv",
        DEG_miRNA_updated = inter_in + "{comp}/DEG_miRNA.csv",
    output:
        interaction_TF_TF = inter_out + "{comp}/interaction_TF_TF.json",
        interaction_TF_GENE = inter_out + "{comp}/interaction_TF_GENE.json",
        interaction_TF_lncRNA = inter_out + "{comp}/interaction_TF_lncRNA.json",
        interaction_TF_miRNA = inter_out + "{comp}/interaction_TF_miRNA.json",
    shell:
        """
        mkdir -p {params.inter}{params.comp}
        python3 {script_dir}/interaction_TF.py \
            {input.DEG_TF} \
            {input.DEG_GENE} \
            {input.DEG_lncRNA} \
            {input.DEG_miRNA_updated} \
            {output.interaction_TF_TF} \
            {output.interaction_TF_GENE} \
            {output.interaction_TF_lncRNA} \
            {output.interaction_TF_miRNA} \
            {params.inter}{params.comp}
        """

# ——————————————————————————————————————————————————————————————
# 6) Constructions des fichiers cytoscape 
rule cytoscape: 
    params:
        comp = "{comp}",
        cyto_in = cyto_in,
        enri_in = enri_in,
    input: 
        DEG_miRNA_updated = inter_in + "{comp}/DEG_miRNA.csv",
        DEG_TF = inter_in + "{comp}/DEG_TF.csv",
        DEG_GENE = inter_in + "{comp}/DEG_GENE.csv",
        dico_id = conv_out + "dico_id.json",

        interaction_miRNA_TF = inter_out + "{comp}/interaction_miRNA_TF.json",
        interaction_miRNA_GENE = inter_out + "{comp}/interaction_miRNA_GENE.json",
        interaction_miRNA_lncRNA = inter_out + "{comp}/interaction_miRNA_lncRNA.json",

        interaction_TF_TF = inter_out + "{comp}/interaction_TF_TF.json",
        interaction_TF_GENE = inter_out + "{comp}/interaction_TF_GENE.json",
        interaction_TF_lncRNA = inter_out + "{comp}/interaction_TF_lncRNA.json",
        interaction_TF_miRNA = inter_out + "{comp}/interaction_TF_miRNA.json",

    output: 
        nodes = cyto_in + "{comp}/nodes.csv",
        edges = cyto_in + "{comp}/edges.csv",

        GO = enri_in + "{comp}/GO.csv",
        GO_TF = enri_in + "{comp}/GO_TF.csv",
        GO_GENE = enri_in + "{comp}/GO_GENE.csv",
        GO_miRNA = enri_in + "{comp}/GO_miRNA.csv",    
    shell: 
        """
        mkdir -p {params.cyto_in}{params.comp} {params.enri_in}{params.comp}
        python3 {script_dir}/cytoscape_FileBuilder.py \
            {input.DEG_miRNA_updated} \
            {input.DEG_TF} \
            {input.DEG_GENE} \
            {input.dico_id} \
            {input.interaction_miRNA_TF} \
            {input.interaction_miRNA_GENE} \
            {input.interaction_miRNA_lncRNA} \
            {input.interaction_TF_TF} \
            {input.interaction_TF_GENE} \
            {input.interaction_TF_lncRNA} \
            {input.interaction_TF_miRNA} \
            {output.nodes} \
            {output.edges} \
            {output.GO} \
            {output.GO_TF} \
            {output.GO_GENE} \
            {output.GO_miRNA} \
            {params.cyto_in}{params.comp} \
            {params.enri_in}{params.comp}
        """


# ——————————————————————————————————————————————————————————————
# 7) Enrichissement fonctionnel 
rule enrichissement: 
    params:
        comp = "{comp}",
        enri = enri_out,
        
    input: 
        GO = enri_in + "{comp}/GO.csv",
        GO_GENE = enri_in + "{comp}/GO_GENE.csv",
        GO_TF = enri_in + "{comp}/GO_TF.csv",   
    output: 
        GO_all = enri_out + "{comp}/GO_all.png", 
        GO_genes = enri_out + "{comp}/GO_genes.png", 
        GO_tfs = enri_out + "{comp}/GO_tfs.png", 

        GO_down_all = enri_out + "{comp}/GO_down_all.png", 
        GO_down_genes = enri_out + "{comp}/GO_down_genes.png", 
        GO_down_tfs = enri_out + "{comp}/GO_down_tfs.png", 

        GO_up_all = enri_out + "{comp}/GO_up_all.png", 
        GO_up_genes = enri_out + "{comp}/GO_up_genes.png", 
        GO_up_tfs = enri_out + "{comp}/GO_up_tfs.png", 
    shell: 
        """
        mkdir -p {params.enri}{params.comp}         
        Rscript {script_dir}/Enrichissement.R \
            {input.GO} \
            {input.GO_GENE} \
            {input.GO_TF} \
            {output.GO_all} \
            {output.GO_genes} \
            {output.GO_tfs} \
            {output.GO_down_all} \
            {output.GO_down_genes} \
            {output.GO_down_tfs} \
            {output.GO_up_all} \
            {output.GO_up_genes} \
            {output.GO_up_tfs} \
            {params.enri}{params.comp}
        """
