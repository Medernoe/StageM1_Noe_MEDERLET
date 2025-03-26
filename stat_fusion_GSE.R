#!/usr/bin/env Rscript

# Définir le répertoire de travail
setwd("/home/medernoe/Ubuntu_user/Stage_M1_ubuntu/StageM1_git/Data/Pre-traitement/Output")

library(car)
library(nortest)


#test stat
#==========================================================================================================================
DATA1_PCC <- read.csv("DATA1_PCC.csv")
DATA1_PGL <- read.csv("DATA1_PGL.csv")
DATA1_CONTROL <- read.csv("DATA1_CONTROL.csv")
DATA2 <- read.csv("DATA2.csv")


# --- ÉTAPE 1 : Préparer les données ---


# vectoriser les données 
DATA2_values <- as.vector(as.matrix(DATA2))
DATA1_PCC_values <- as.vector(as.matrix(DATA1_PCC))




#création d'un dataframe unique 
data_total <- data.frame(
  intensity = c(DATA2_values, DATA1_PCC_values),
  groupe = factor(rep(c("DATA2", "DATA1_PCC"), 
                      times = c(length(DATA2_values), length(DATA1_PCC_values))))
)


head(data_total)
summary(data_total)

mean(data_total$intensity[data_total$groupe == "DATA2"])
mean(data_total$intensity[data_total$groupe == "DATA1_PCC"])

summary(data_total$intensity[data_total$groupe == "DATA2"])
summary(data_total$intensity[data_total$groupe == "DATA1_PCC"])


#création d'un sample pour le teste de normalité 
# Définir une graine (seed) pour rendre l'échantillonnage reproductible
set.seed(12345) 

# Sélectionner aléatoirement 5000 valeurs dans le jeu de données 
DATA1_PCC_sample <- sample(DATA1_PCC_values, size = 5000, replace = FALSE)

DATA2_sample <- sample(DATA2_values, size = 5000, replace = FALSE)




# --- ÉTAPE 2 : Tester la Normalité (Shapiro-Wilk) ---


# Test de normalité pour DATA2
shapiro_DATA2 <- shapiro.test(DATA2_sample) # Shapiro limité à 5000 valeurs
shapiro_DATA2
hist(DATA2_sample)
hist(DATA2_values)

# Test de normalité pour DATA1_PCC
shapiro_DATA1 <- shapiro.test(DATA1_PCC_sample)
shapiro_DATA1 # 2.2e-16 
hist(DATA1_PCC_sample)
hist(DATA1_PCC_values)

# Interprétation :
# Si p-value > 0.05 : pas assez d’évidence pour rejeter la normalité.
# Si p-value < 0.05 : rejet de l’hypothèse de normalité  

#les données semblent normales visuellement cependant les test montre une pvalue = 2.2e-16 donc données pas normales ? 




# --- ÉTAPE 2 - bis  : Tester la Normalité (Anderson-Darling normality test) + QQ-plot ---

#ce test peut etre effectué sur un jeu de données plus grand, pas de limite à 5000 et il est plus sensible aux valeurs extrèmes 

# Exemple :
ad.test(DATA1_PCC_values)
ad.test(DATA2_values)

#Si p-value < 0.05, rejet de l’hypothèse de normalité 
#Si p-value > 0.05, non rejet de l'hypothèse de normalité.

#les données ne semble pas normale non plus meme avec ce test 



# Vérifier visuellement par QQ-plot
qqnorm(DATA1_PCC_sample, main="QQ-plot DATA1_PCC_sample")
qqline(DATA1_PCC_sample, col = "red")

# Vérifier visuellement par QQ-plot
qqnorm(DATA2_sample, main="QQ-plot DATA2_sample")
qqline(DATA2_sample, col = "red")

#Centre (valeurs moyennes) : les points suivent relativement bien la droite théorique (ligne rouge).
#Aux extrémités : nettes déviations des points par rapport à la droite de normalité,indique une distribution avec des valeurs extrêmes plus importantes.
#(queues lourdes)




## Étape 3 : Réalisation du Test non paramétrique (Wilcoxon-Mann-Whitney)


# Exécution du test de Wilcoxon (Mann-Whitney)
wilcox_result <- wilcox.test(DATA2_values, DATA1_PCC_values)
print(wilcox_result)


## Étape 4 : Interprétation des résultats

#Si p-value < 0.05 : 
#Rejet h0 : les distributions (intensités) diffèrent significativement entre DATA2 et DATA1_PCC.

#Si p-value ≥ 0.05 :
#On ne rejette pas h0 . Conclusion : pas de différence significative détectée.


# Interprétation automatique
if(wilcox_result$p.value < 0.05){
  cat("✅ Différence significative entre DATA2 et DATA1_PCC (p-value =", wilcox_result$p.value,")\n")
} else {
  cat("❌ Pas de différence significative détectée (p-value =", wilcox_result$p.value,")\n")
}



## Étape 3-bis : Réalisation du t-Test 

var.test(DATA2_values, DATA1_PCC_values)
#p-value < 2.2e-16 les variance ne sont pas égale 


t.test(DATA2_values, DATA1_PCC_values, var.equal = FALSE)#p-value = 0.182 
#les données ne présente pas de différence significative, on ne rejette pas H0 

ks.test(DATA2_values, DATA1_PCC_values)#p-value < 2.2e-16
#les données présente des différences significatives 
