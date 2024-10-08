rm(list=ls())
rm(list = ls(all.names = TRUE), envir = .rs.CachedDataEnv)
# installer et charger les packages necessaires
library(dplyr)
library(tidyr)
#install.packages("ggplot2")
library(ggplot2)

# récupérer les fichiers avec les données / pollens identifiés et nettoyer

#data cyto
dataCY<- read.csv("ID_2023_2140_CY.csv", sep=",", h=T)
group<-subset(dataCY,select=-X)
dataCY$TOTAL<-rowSums(group)
dataCY$Sample<-dataCY$X
dataCY<-subset(dataCY, select=-c(X,Debris))
dataCY <- pivot_longer(dataCY, cols = -c(Sample, TOTAL), names_to = "Taxon", values_to = "Nombre")
dataCY[is.na(dataCY)] <- 0
numeros <- as.numeric(gsub("CY_", "", dataCY$Sample))
dataCY <- dataCY[order(numeros), ]
dataCY$Sample <- factor(dataCY$Sample, levels = unique(dataCY$Sample))
#write.csv(dataCY, "dataCY_plot.csv", row.names = T)

#data micro
dataMI<- read.csv("ID_2023_2140_MI.csv", sep=";", h=T)
names(dataMI) <- gsub("Total\\.count\\.", "", names(dataMI))
group<-subset(dataMI,select=-Sample)
dataMI$TOTAL<-rowSums(group)
#dataMI<-subset(dataMI, select=-c(no.pollen.found))
dataMI <- pivot_longer(dataMI, cols = -c(Sample, TOTAL), names_to = "Taxon", values_to = "Nombre")
dataMI[is.na(dataMI)] <- 0
numeros1 <- as.numeric(gsub("MI_", "", dataMI$Sample))
dataMI <- dataMI[order(numeros1), ]
dataMI$Sample <- factor(dataMI$Sample, levels = unique(dataMI$Sample))

# faire graphique data cyto
ggplot(dataCY, aes(x = Sample, y = Nombre)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Nombre total de pollens pour la cyto")+
  coord_cartesian(ylim = c(0, 10000))
# faire graphique data micro
ggplot(dataMI, aes(x = Sample, y = Nombre)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Nombre total de pollens pour la micro")+
  coord_cartesian(ylim = c(0, 5000))

# regrouper tous les tableaux en un seul
dataCY$meth<-"cyto"
dataMI$meth<-"micro"
dataCY$Sample <- gsub("CY_", "", dataCY$Sample)
dataMI$Sample <- gsub("MI_", "", dataMI$Sample)
all_data<-full_join(dataCY,dataMI)


#stats pour abondances totales entre cyto et micro de chaque echantillon
# Préparer les données
all_data$Sample <- as.factor(all_data$Sample)
all_data$Taxon <- as.factor(all_data$Taxon)
all_data$meth <- as.factor(all_data$meth)
#write.csv(all_data, "all_data_plotMICY.csv", row.names = T)
#extraire données TOTAL pour chaque meth et chaque Sample
TOTAL <- data.frame(Sample = numeric(), meth = character(), TOTAL = numeric())
for (sample in unique(all_data$Sample)) {
  for (meth_type in unique(all_data$meth)) {
    unique_TOTAL <- unique(all_data$TOTAL[all_data$Sample == sample & all_data$meth == meth_type])
    if (length(TOTAL) > 0) {
      TOTAL <- rbind(TOTAL, data.frame(Sample = sample, meth = meth_type, TOTAL = unique_TOTAL))
    }
  }
}

rownames(TOTAL) <- NULL
print(TOTAL)
ggplot(TOTAL, aes(x = Sample, y = TOTAL, fill = meth)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparaison des abondances totales de pollens 
en fonction de la méthode d'identification
et de l'échantillon")+
  coord_cartesian(ylim = c(0, 125000))

ggplot(TOTAL, aes(x = meth, y = TOTAL, fill = meth)) +
  geom_boxplot() +
  labs(title = "Comparaison des abondances totales de pollens 
en fonction de la méthode d'identification",
       x = "Méthode",
       y = "TOTAL") +
  coord_cartesian(ylim = c(0, 60000))

# stats pour abondances relatives des taxons 
ggplot(all_data[all_data$Sample %in% unique(all_data$Sample)[1:70], ], aes(x = Taxon, y = Nombre, fill=meth)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Sample) +
  labs(title = "Comparaison du nombre de pollens de chaque taxon en fonction méthodes d'ID")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  coord_cartesian(ylim = c(0, 5000))

#t-test pour chaque taxon
donnees_cyto <- subset(all_data, meth == "cyto")
donnees_micro <- subset(all_data, meth == "micro")
ttest_taxon <- data.frame(Taxon = character(), p_value = numeric(), difference = numeric(), stringsAsFactors = FALSE)
taxons <- unique(all_data$Taxon)
for (taxon in taxons) {
  donnees_cyto_taxon <- subset(donnees_cyto, Taxon == taxon)
  donnees_micro_taxon <- subset(donnees_micro, Taxon == taxon)
  diff_mesures <- mean(donnees_micro_taxon$Nombre) - mean(donnees_cyto_taxon$Nombre)
  resultat <- t.test(donnees_cyto_taxon$Nombre, donnees_micro_taxon$Nombre, paired = TRUE)
  ttest_taxon <- rbind(ttest_taxon, data.frame(Taxon = taxon, p_value = resultat$p.value, difference = diff_mesures))
}
print(ttest_taxon)
getwd()
write.csv(ttest_taxon, "ttest_taxon_2023pt2.csv", row.names = T)

#zoom sur certains échantillons
ggplot(all_data[all_data$Sample == "28", ], aes(x = Taxon, y = Nombre, fill=meth)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Sample) +
  labs(title = "Comparaison du nombre de pollens de chaque taxon en fonction méthodes d'ID")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  coord_cartesian(ylim = c(0, 2500))
