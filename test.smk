# -*- Snakefile -*-

# Pour pouvoir passer data_dir / input_dir / output_dir via --config
DATA_DIR   = config.get("data_dir",   "Data2/")
INPUT_DIR  = config.get("input_dir",  "Input/")
OUTPUT_DIR = config.get("output_dir", "Output/")

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
# Liste des comparaisons à faire 
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
        treatment = lambda wc: pre_out + f"gene_{COMPARISONS[wc.comp][1]}.csv",
    output:
        results = deg_out + "{comp}/DEG_results_gene.csv",
        plot    = deg_out + "{comp}/volcano_plot_gene.png",
    shell:
        """
        mkdir -p {deg_out}{params.comp}
        Rscript {script_dir}/DEG_gene.R \
            {input.control} \
            {input.treatment} \
            {output.results} \
            {output.plot}
        """

# DEG — miRNA
rule DEG_miRNA:
    params:
        comp="{comp}"
    input:
        control   = lambda wc: pre_out + f"miRNA_{COMPARISONS[wc.comp][0]}.csv",
        treatment = lambda wc: pre_out + f"miRNA_{COMPARISONS[wc.comp][1]}.csv",
    output:
        results = deg_out + "{comp}/DEG_results_miRNA.csv",
        plot    = deg_out + "{comp}/volcano_plot_miRNA.png",
    shell:
        """
        mkdir -p {deg_out}{params.comp}
        Rscript {script_dir}/DEG_miRNA.R \
            {input.control} \
            {input.treatment} \
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
           {params.conv_out} \
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
        mkdir -p {trans_in} {params.comp}
        Rscript {script_dir}/conversion_SpotToGene.R \
            {input.update_gene_mapping} \
            {input.dico_id} \
            {input.results} \
            {output.DEG_updated_gene} \
            {params.trans}
        """
