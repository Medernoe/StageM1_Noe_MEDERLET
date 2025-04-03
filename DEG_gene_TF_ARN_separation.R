#!/usr/bin/env Rscript

# Définir le répertoire de travail
setwd("/home/medernoe/Ubuntu_user/Stage_M1_ubuntu/StageM1_git/Data")


data = read.csv("DEG_complete.csv",row.names = 1 )

table(data$Type)

sum(data$Symbol == "")
#les 146 unknown ne correspondent pas à des gènes 


#====================================================================================================

data = data[data$Symbol != "", ]

DEG_TF = data[data$Type == "TF", ]

DEG_GENE = data[data$Type == "protein_coding", ]

DEG_lncRNA = data[data$Type == 'lncRNA', ]

