# =============================================================================
# Pipeline d'analyse d'enrichissement GO et cibles de miRNA
#
# Description :
#   Réalise un enrichissement Gene Ontology (Processus Biologique) avec clusterProfiler
#   pour des listes de gènes, de facteurs de transcription (TF) et pour les cibles validées
#   de microARNs (miRNA) via multiMiR.
#   Génère des sous-listes régulées à la hausse/à la baisse et affiche tous les graphiques.
#
# Utilisation :
#   1) Installer et charger les packages requis (à faire une seule fois).
#   2) Définir les chemins des fichiers d'entrée et ajuster les paramètres si nécessaire.
#   3) Charger ce script pour exécuter l'ensemble du pipeline.
#
# =============================================================================

# 1) Installation des packages Bioconductor (à exécuter une seule fois)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "clusterProfiler",   # Enrichissement GO
  "org.Hs.eg.db",      # Base SYMBOL ↔ ENTREZID
  "multiMiR",          # Récupération des cibles de miRNA
  "enrichplot",        # Visualisations d'enrichissement
  "ggplot2"            # Graphiques
))

# 2) Chargement des bibliothèques
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(multiMiR)

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
input_dir   <- "Data/Enrichissement/Input"

data_all    <- read.csv(file.path(input_dir, "GO.csv"),      header = TRUE)
data_genes  <- read.csv(file.path(input_dir, "GO_GENE.csv"), header = TRUE)
data_tfs    <- read.csv(file.path(input_dir, "GO_TF.csv"),   header = TRUE)
data_mir    <- read.csv(file.path(input_dir, "GO_miRNA.csv"), header = TRUE)

id_all   <- data_all$Gene_id
id_genes <- data_genes$Gene_id
id_tf   <- data_tfs$Gene_id
id_mir   <- data_mir$Symbol        

# =============================================================================
# 5) Enrichissement GO pour gènes et TFs
results <- list()

res_all   <- run_go_enrichment(id_all,   "Tous gènes + TFs")
res_genes <- run_go_enrichment(id_genes, "Gènes seuls")
res_tfs   <- run_go_enrichment(id_tf,   "TFs seuls")

#stock tous les résultats dans results 
results <- c(results, all = list(res_all), genes = list(res_genes), tfs = list(res_tfs))

# 6) Enrichissements up/down (si logFC)

results$up_all   <- run_go_enrichment(data_all$Gene_id[data_all$logFC>0],   "Tous – régulés à la hausse")
results$down_all <- run_go_enrichment(data_all$Gene_id[data_all$logFC<0],   "Tous – régulés à la baisse")

results$up_genes   <- run_go_enrichment(data_genes$Gene_id[data_genes$logFC>0], "Gènes – régulés à la hausse")
results$down_genes <- run_go_enrichment(data_genes$Gene_id[data_genes$logFC<0], "Gènes – régulés à la baisse")

results$up_tfs   <- run_go_enrichment(data_tfs$Gene_id[data_tfs$logFC>0], "TFs – régulés à la hausse")
results$down_tfs <- run_go_enrichment(data_tfs$Gene_id[data_tfs$logFC<0], "TFs – régulés à la baisse")


# =============================================================================
# 7) Récupération et enrichissement des cibles miRNA
# Récupérer les interactions validées
mm <- get_multimir(mirna = id_mir, table = "validated", summary = FALSE)
target_ids <- unique(mm@data$target_entrez)

# Enrichissement GO des cibles
res_mir <- run_go_enrichment(target_ids, "Cibles de miRNA")
results <- c(results, mir = list(res_mir))

# =============================================================================
# 8) Affichage des résultats
message("Analyse d'enrichissement GO terminée. Affichage de tous les graphiques...")
for (nm in names(results)) {
  cat("-- Graphique :", nm, "--\n")
  print(results[[nm]]$plot)
}

# =============================================================================
# 9) Sauvegarde  

output_dir <- "Data/Enrichissement/Output"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (nm in names(results)) {
  filename <- file.path(output_dir, paste0("GO_", nm, ".png"))
  ggsave(filename = filename,
         plot     = results[[nm]]$plot,
         width    = 14,
         height   = 12,
         dpi      = 500)
}

