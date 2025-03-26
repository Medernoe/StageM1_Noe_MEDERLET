#!/usr/bin/env Rscript

# Définir le répertoire d'input
setwd("/home/medernoe/Ubuntu_user/Stage_M1_ubuntu/StageM1_git/Data/Pre-traitement/Output")

# Charger les packages nécessaires
library(FactoMineR)
library(factoextra)
library(car)

# Import des données
PCC_19422 <- read.csv("DATA19422_PCC.csv")
colnames(PCC_19422) <- gsub("^GSM([0-9]+)$", "PCC\\1", colnames(PCC_19422))

PCC_51081 <- read.csv("DATA51081.csv")
colnames(PCC_51081) <- gsub("^GSM([0-9]+)$", "2PCC\\1", colnames(PCC_51081))

Control <- read.csv("DATA19422_CONTROL.csv")
colnames(Control) <- gsub("^GSM([0-9]+)$", "Control\\1", colnames(Control))


# Définir le répertoire de sortie des résultats
setwd("/home/medernoe/Ubuntu_user/Stage_M1_ubuntu/StageM1_git/Data")

# Vérifier si le répertoire Output existe
if (!dir.exists("Stats/ACP")) {
  # Créer le répertoire s'il n'existe pas
  dir.create("Stats/ACP", recursive = TRUE)
}


setwd("/home/medernoe/Ubuntu_user/Stage_M1_ubuntu/StageM1_git/Data/Stats/ACP")


# Fusion des données en une seule matrice
temp_merge <- merge(PCC_19422, PCC_51081, by = "X", all = TRUE)
data <- merge(temp_merge, Control, by = "X", all = TRUE)

rownames(data) <- data[, 1]
data <- data[, -1]

# Statistiques descriptives
summary(data)

# PCA
par(cex.axis = 0.5)
res <- PCA(data, scale = TRUE)


# Boxplots des variables quantitatives
X11()  # Ouvrir un nouveau dispositif graphique
par(cex.axis = 1.5)  # Réduire la taille de police des axes
plot1 <- boxplot(data[, 1:76], las = 2, cex.axis = 0.7, main = "Boxplots des variables quantitatives")


# Graphique d'inertie
dev.new()
par(cex.axis = 0.5)
plot2 <- barplot(res$eig[, 2], main = "Graphique d'inerties des axes")


# Nuages des variables 
dev.new()
plot3 <- fviz_pca_var(res, col.var = "cos2",
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE,
                      labelsize = 2)
print(plot3)


# Coloration par groupe
categories <- c(rep("PCC", length(grep("^PCC", colnames(data)))),
                rep("2PCC", length(grep("^2PCC", colnames(data)))),
                rep("Control", length(grep("^Control", colnames(data)))))

colors <- c("PCC" = "blue", "2PCC" = "green", "Control" = "red")


# Visualisations colorées
dev.new()
plot4 <- fviz_pca_var(res,
                      col.var = categories,
                      palette = colors,
                      repel = TRUE,
                      labelsize = 3,
                      axes = c(1, 2))
print(plot4)


dev.new()
plot5 <- fviz_pca_var(res,
                      col.var = categories,
                      palette = colors,
                      repel = TRUE,
                      labelsize = 3,
                      axes = c(1, 3))

print(plot5)


# dev.off() pour fermer 
