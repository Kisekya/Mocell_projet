# PROJET MODELISATION QUALITATIVE 2019 
# Amelie GRUEL - Coralie MULLER
# octobre 2019

setwd("~/M2/MOCELL/Mocell_projet")

# obtenir les données d'intérêt dans la base de données expérimentales
data_experimentale=read.table("TCGA_rsem_norm_Nature2012_RNAseq.tab",header=T)
genes_focus<-read.table("GOI_ids.txt",header=F)
genes_focus<-as.vector(genes_focus[,1])
filter_list<-data_experimentale[,"ID"] %in% genes_focus  # perte d'informations sur fichiers
data_focus<-data_experimentale[filter_list,]  # data_focus contient les données d'intérêt

# modifier data_focus pour récupérer les labels des lignes en première colonne
rownames(data_focus)=data_focus[,1]
data_focus=data_focus[,-1]

# corrélation afin d'obtenir le réseau de co-expression
data_focus=t(data_focus)
correlation=cor(data_focus)
# install.packages("reshape2")
library("reshape2",lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3/")
correlation=melt(correlation)

corr = correlation[correlation["value"]>0.8,]
if(length(corr[,1])>0) corr["type"]="correlation_pair"
anticorr = correlation[correlation["value"]<-0.8,]
if (length(anticorr[,1])!=0) anticorr["type"]="anticorrelation_pair"

if (length(anticorr[,1])>0 && length(corr)>0) {
  correlation=rbind(corr,anticorr)
} else if (length(corr[,1])>0) {
  correlation=corr
} else if (length(anticorr[,1])>0) {
  correlation=anticorr
} else {
  correlation=matrix(, nrow = 0, ncol = 3)
}

correlation["Alphabétique"]<-as.character(correlation[,"Var1"])<as.character(correlation[,"Var2"])
correlation=correlation[correlation[,"Alphabétique"]==TRUE,]
correlation=correlation[,c(1,2,4)]
colnames(correlation)[1:2]=c("id1","id2")
write.table(correlation,"correlation.txt",sep="\t",row.names=F,quote=F)

# donner le type d'interactions à partir des données PPI (Human_STRING_PPI.tab)
data_PPI=read.table("Human_STRING_PPI.tab",header=T)

# filtre sur les PPI à partir de la liste de corrélation
genes_interets=unique(c(as.vector(correlation[,"id1"]),as.vector(correlation[,"id2"])))
filter_list<-data_PPI[,"id1"] %in% genes_interets 
data_PPI<-data_PPI[filter_list,]
filter_list<-data_PPI[,"id2"] %in% genes_interets
data_PPI<-data_PPI[filter_list,]

# fusionner les données PPI et les données de corrélation sur notre liste d'intérêt en un seul dataframe
data_interet=rbind(correlation,data_PPI)
write.table(data_interet,"data_interet.txt",sep="\t",row.names=FALSE,quote=FALSE)

# sélection d'information sur les TF (Human_TF_LABEL.tab)
data_TF=read.table("Human_TF_LABEL.tab",header=T,sep="\t")
filter_list<-data_TF[,"id"] %in% genes_interets
data_TF<-data_TF[filter_list,]
write.table(data_TF,"TF.txt",sep="\t",row.names=F,quote=F)



# filtrer dans PAM50 ceux appartenant dans GOI par True
data_PAM=read.table("PAM50_geneIDs.txt",header=T)
filter_list<-data_PAM[,"id"] %in% genes_focus 
data_PAM["presence"]=filter_list
write.table(data_PAM,"PAM.txt",sep="\t",row.names=F,quote=F)