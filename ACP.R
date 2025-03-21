library(FactoMineR)
library(factoextra)
library(car)  # Permet de tracer des ellipses autour des groupes


#import des données---------------------------------------------------------------------------------------------------------- 

PCC_19422 = read.csv("DATA19422_PCC.csv")
colnames(PCC_19422) <- gsub("^GSM([0-9]+)$", "PCC\\1", colnames(PCC_19422))

PCC_51081 = read.csv("DATA51081.csv")
colnames(PCC_51081) <- gsub("^GSM([0-9]+)$", "2PCC\\1", colnames(PCC_51081))

Control = read.csv("DATA19422_CONTROL.csv")
colnames(Control) <- gsub("^GSM([0-9]+)$", "Control\\1", colnames(Control))


#fusion des données en une seul matrice -------------------------------------------------------------------------------------

# Fusionner en utilisant une jointure externe
temp_merge <- merge(PCC_19422, PCC_51081, by = "X", all = TRUE)
data <- merge(temp_merge, Control, by = "X", all = TRUE)

rownames(data) <- data[, 1]  
data <- data[, -1]



#Stat descriptive--------------------------------------------------------------------------------------------------------------
summary(data)

par(cex.axis = 1.5)       # Réduire la taille de police des axes
boxplot(data[,1:76],
        las = 2,
        cex.axis = 0.7,
        main = "Boxplots des variables quantitatives")

#répartitions similaire, meme echelle de données, semble proche sauf les PCC des 2 expériences, et PCC 483013

#PCA-----------------------------------------

par(cex.axis = 0.5)  
res=PCA(data,scale=T)

#inertie-----------------------------------------
par(cex.axis = 0.5)    
barplot(res$eig[,2],main="graphique d'inerties des axes")

#la dimension 1 et 2 semble suffisante pour décrire le jeux de données 

plot(res, choix = "var", axes = c(1, 2))


#beau graph des var 1:2 
fviz_pca_var(res, col.var = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Bleu → Orange → Rouge
             repel = TRUE,
             labelsize = 2)


#les gèns sont tous répartie au centre du graph 
#Les données sont difficilement lisible, cependant il ne semblerait pas y avoir de cluster dérangeant dans notre jeux de données 




#coloration par groupe --------------------------------------------------------------------------------------------------------
# Personnalisation de l'affichage des variables sur le graphique PCA
# Colorier les variables en fonction des groupes (2PCC, PCC, Control)
categories <- c(rep("PCC", length(grep("^PCC", colnames(data)))),
                rep("2PCC", length(grep("^2PCC", colnames(data)))),
                rep("Control", length(grep("^Control", colnames(data)))))

# Attribuer des couleurs spécifiques aux groupes
colors <- c("PCC" = "blue", "2PCC" = "green", "Control" = "red")

# Visualisation avec fviz_pca_var
fviz_pca_var(res, 
             col.var = categories,  # Couleur par groupe
             palette = colors,      # Couleurs spécifiées
             repel = TRUE,
             labelsize = 3,         # Ajuster la taille des labels
             axes = c(1, 2))        # Choisir les axes 1 et 2

fviz_pca_var(res, 
             col.var = categories,  # Couleur par groupe
             palette = colors,      # Couleurs spécifiées
             repel = TRUE,
             labelsize = 3,         # Ajuster la taille des labels
             axes = c(1, 3))        # Choisir les axes 1 et 


#en dim 1 et 2 on dirait déja que c'est un peu différents, mais en dim 1 et 3 on dirait encore plus 


#CAH--------------------------------------------------------------------------------------------------------------------------------
# Assurez-vous que le jeu de données est bien normalisé si nécessaire
# Si les données ne sont pas déjà normalisées, vous pouvez les standardiser avant d'effectuer le K-means.
data_scaled <- scale(data)

# Appliquer l'algorithme K-means avec k = 3
set.seed(123)  # Pour assurer la reproductibilité des résultats
kmeans_result <- kmeans(data_scaled, centers = 3, nstart = 25)

# Afficher les résultats
print(kmeans_result)

# Visualiser les clusters sur les 2 premières dimensions (PCA pour la réduction de dimension)
# Si vous souhaitez visualiser les clusters dans un graphique en 2D, utilisez la PCA pour réduire la dimension.
res_pca <- prcomp(data_scaled)
pca_data <- as.data.frame(res_pca$x)

# Ajouter les résultats de K-means à la PCA pour visualiser les clusters
pca_data$cluster <- as.factor(kmeans_result$cluster)

# Tracer les clusters
library(ggplot2)
ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() +
  labs(title = "K-means Clustering (k = 3)", x = "PC1", y = "PC2") +
  scale_color_manual(values = c("red", "blue", "green"))

