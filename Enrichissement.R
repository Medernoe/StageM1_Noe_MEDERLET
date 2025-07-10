# =============================================================================
# Pipeline d'analyse d'enrichissement GO et cibles de miRNA
#
# Description :
#   Réalise un enrichissement Gene Ontology (Processus Biologique) avec clusterProfiler
#   pour des listes de gènes, de facteurs de transcription (TF) et pour les cibles validées
#   de microARNs (miRNA) via multiMiR (retiré mais peut etre remplacé par les cibles dans miRabel).
#   Génère des sous-listes régulées à la hausse/à la baisse et affiche tous les graphiques.
#
# Utilisation :
#   1) Installer et charger les packages requis (à faire une seule fois).
#   2) Définir les chemins des fichiers d'entrée et ajuster les paramètres si nécessaire.
#   3) Charger ce script pour exécuter l'ensemble du pipeline.
#
# =============================================================================

# 1) Lecture des arguments
args <- commandArgs(trailingOnly = TRUE)

# =============================================================================

# 2) Chargement des bibliothèques
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})



# =============================================================================
# 3) Fonctions 

# Lance l'enrichissement GO-BP et retourne un objet (résultat + graphique)
run_go_enrichment <- function(entrez_ids, label, org_db = org.Hs.eg.db,
                              ontology = "BP", p_adj = "BH", q_cutoff = 0.05,
                              top_n = 20) {
  
  ego <- enrichGO(gene = entrez_ids, OrgDb = org_db, ont = ontology,
                  pAdjustMethod = p_adj, qvalueCutoff = q_cutoff, readable = TRUE)
  
  plot <- dotplot(ego, showCategory = top_n) +
    ggtitle(paste0("Enrichissement GO-BP – ", label))
  
  list(result = ego, plot = plot)
}



# =============================================================================
# 4) Données d'entrée
# Chaque CSV doit contenir les colonnes 'Symbol' et 'logFC'.

#ouverture des fichiers d'entrée
input_GO      <- args[1]
input_GO_GENE <- args[2]
input_GO_TF   <- args[3]  

data_all    <- read.csv(input_GO,      header = TRUE)
data_genes  <- read.csv(input_GO_GENE, header = TRUE)
data_tfs    <- read.csv(input_GO_TF,   header = TRUE)


id_all   <- data_all$Gene_id
id_genes <- data_genes$Gene_id
id_tf   <- data_tfs$Gene_id
    
# =============================================================================
# 5) Enrichissement GO pour gènes et TFs
results <- list()

res_all   <- run_go_enrichment(id_all,   "Tous gènes + TFs")
res_genes <- run_go_enrichment(id_genes, "Gènes seuls")
res_tfs   <- run_go_enrichment(id_tf,   "TFs seuls")

#stock tous les résultats dans results 
results <- c(results, all = list(res_all), genes = list(res_genes), tfs = list(res_tfs))

# 6) Enrichissements up/down (si logFC)

results$down_all <- run_go_enrichment(data_all$Gene_id[data_all$logFC<0],   "Tous – régulés à la baisse")
results$down_genes <- run_go_enrichment(data_genes$Gene_id[data_genes$logFC<0], "Gènes – régulés à la baisse")
results$down_tfs <- run_go_enrichment(data_tfs$Gene_id[data_tfs$logFC<0], "TFs – régulés à la baisse")

results$up_all   <- run_go_enrichment(data_all$Gene_id[data_all$logFC>0],   "Tous – régulés à la hausse")
results$up_genes   <- run_go_enrichment(data_genes$Gene_id[data_genes$logFC>0], "Gènes – régulés à la hausse")
results$up_tfs   <- run_go_enrichment(data_tfs$Gene_id[data_tfs$logFC>0], "TFs – régulés à la hausse")



# =============================================================================
# 7) Affichage des résultats
print("Analyse d'enrichissement GO terminée. Affichage de tous les graphiques...")
for (nm in names(results)) {
  cat("-- Graphique :", nm, "--\n")
  print(results[[nm]]$plot)
}

# =============================================================================
# 8) Sauvegarde  
output_enri <- args[13]
# Vérifier si le répertoire Output existe
if (!dir.exists(output_enri)) {
  # Créer le répertoire s'il n'existe pas
  dir.create(output_enri, recursive = TRUE)
}

#chemin d'accès des fichiers de sortie
output_GO_all   <- args[4]
output_GO_genes <- args[5]
output_GO_tfs <- args[6]

output_GO_down_all <- args[7]
output_GO_down_genes <- args[8]
output_GO_down_tfs <- args[9]

output_GO_up_all <- args[10]
output_GO_up_genes <- args[11]
output_GO_up_tfs <- args[12]

# Sauvegarde explicite de chaque plot avec le nom de fichier correspondant
ggsave(filename = output_GO_all,   plot = results$all$plot,   width = 14, height = 12, dpi = 500)
ggsave(filename = output_GO_genes, plot = results$genes$plot, width = 14, height = 12, dpi = 500)
ggsave(filename = output_GO_tfs,   plot = results$tfs$plot,   width = 14, height = 12, dpi = 500)

ggsave(filename = output_GO_down_all,   plot = results$down_all$plot,   width = 14, height = 12, dpi = 500)
ggsave(filename = output_GO_down_genes, plot = results$down_genes$plot, width = 14, height = 12, dpi = 500)
ggsave(filename = output_GO_down_tfs,   plot = results$down_tfs$plot,   width = 14, height = 12, dpi = 500)

ggsave(filename = output_GO_up_all,   plot = results$up_all$plot,   width = 14, height = 12, dpi = 500)
ggsave(filename = output_GO_up_genes, plot = results$up_genes$plot, width = 14, height = 12, dpi = 500)
ggsave(filename = output_GO_up_tfs,   plot = results$up_tfs$plot,   width = 14, height = 12, dpi = 500)


message("Fichiers d'enrichissement GO sauvegardés avec succès dans :", output_enri)

