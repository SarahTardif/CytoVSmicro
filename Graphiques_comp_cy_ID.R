rm(list=ls())
rm(list = ls(all.names = TRUE), envir = .rs.CachedDataEnv)
# installer et charger les packages necessaires
library(dplyr)
library(tidyr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("httpgd")



# récupérer les fichiers avec les données / pollens identifiés et nettoyer

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
dataMI<-dataMI[dataMI$Sample!=21,]
dataMI<-dataMI[dataMI$Sample!=32,]

#data cyto
#dataCY<- read.csv("ID_Ech2140_balanced.csv", sep=",", h=T)
#dataCY$Sample<-dataCY$X
#dataCY<-subset(dataCY, select=-c(X,Debris))
#group<-subset(dataCY,select=-Sample)
#dataCY$TOTAL<-rowSums(group)


# ajout colonnes dont pas de grains de pollens dans ces taxa
Taxa_wanted <- c("Aesculus_", "Acer_", "Ambrosia_", "Betula_", "Carya_", "Elaeagnus_", "Fagus_", "Fraxinus_","Ginkgo_","Gleditsia_","Juglans_","Ligustrum_","Malus_","Quercus_","Taxus_")
for (col in Taxa_wanted) {
  if (!col %in% colnames(dataCY)) {
    dataCY[[col]] <- 0
  }
}
dataCY <- pivot_longer(dataCY, cols = -c(Sample, TOTAL), names_to = "Taxon", values_to = "Nombre")
#dataCY[is.na(dataCY)] <- 0
#numeros <- as.numeric(gsub("ID_CY ", "", dataCY$Sample))
#dataCY <- dataCY[order(numeros), ]
#dataCY<-dataCY[dataCY$Sample!="ID_CY 21",]
dataCY$Sample <- factor(dataCY$Sample, levels = unique(dataCY$Sample))

#regroupement en fonction genre et famille pour certains groupes (comme dans dataMI)
dataCY$Genus<-sub("_.*", "", dataCY$Taxon)
dataCY <- aggregate(Nombre ~ Genus + Sample+TOTAL, data = dataCY, sum)
spref<-read.csv("species_refcollection.csv")
spref<-subset(spref,select=-species)
dataCY<-merge(dataCY,spref,by="Genus")
dataCY<- unique(dataCY)
dataCY$Genus <- ifelse(dataCY$Family == "Rosaceae", "Rosaceae", dataCY$Genus)
dataCY$Genus <- ifelse(dataCY$Family == "Pinaceae", "Pinaceae", dataCY$Genus)
dataCY$Genus <- ifelse(dataCY$Family == "Cupressaceae", "Cupressaceae", dataCY$Genus)
dataCY$Genus <- ifelse(dataCY$Genus == "Ostrya"|dataCY$Genus =="Corylus", "Corylus.Ostrya", dataCY$Genus)
colnames(dataCY)[colnames(dataCY) == "Genus"] <- "Taxon"
dataCY<-subset(dataCY,select=-Family)
#write.csv(dataCY, "dataCY_plot.csv", row.names = T)

# renettoyage de dataMI et dataCY pour avoir les memes catégories dans taxons
dataCY$Taxon <- ifelse(dataCY$Taxon == "Carpinus"|dataCY$Taxon =="Celtis"|dataCY$Taxon =="Robinia", "Others", dataCY$Taxon)
dataCY$Taxon <- ifelse(dataCY$Taxon == "Syringa", "Olaeceae", dataCY$Taxon)
dataMI$Taxon <- ifelse(dataMI$Taxon == "Larix", "Pinaceae", dataMI$Taxon)
dataMI$Taxon <- ifelse(dataMI$Taxon == "Armoise"|dataMI$Taxon =="Plantago"|dataMI$Taxon =="Autres.herbes"|dataMI$Taxon =="NI"|dataMI$Taxon =="Gramineae"|dataMI$Taxon =="Typha", "Others", dataMI$Taxon)


# faire graphique data cyto
ggplot(dataCY, aes(x = Sample, y = Nombre)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Nombre total de pollens pour la cyto")+
  coord_cartesian(ylim = c(0, 20000))
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
ggplot(all_data, aes(x = meth, y = TOTAL, fill = meth)) +
  geom_boxplot() +
  #labs(title = "Comparaison des abondances totales de pollens 
#en fonction de la méthode d'identification",
 #      x = "Méthode",
  #     y = "TOTAL") +
  #facet_wrap(~ Sample)+
  coord_cartesian(ylim = c(0, 3000))
summary(all_data$TOTAL[all_data$meth == "micro"])
summary(all_data$TOTAL[all_data$meth == "cyto"])

# stats pour abondances relatives des taxons 
# vérifier que les données dataCY sont regroupées dans meme catégories que pour dataMI sinon pas comparable !
ggplot(all_data, aes(x = Taxon, y = Nombre, fill=meth)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Sample) +
  labs(title = "Comparaison du nombre de pollens de chaque taxon en fonction méthodes d'ID")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  coord_cartesian(ylim = c(0, 1000))
#même chose mais en considérant les échantillons comme des réplicas 
ggplot(all_data, aes(x = Taxon, y = Nombre, fill=meth)) +
  geom_boxplot() +
  labs(title = "Distribution du Nombre par Taxon et Méthode",
       x = "Taxon et Méthode",
       y = "Nombre")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  coord_cartesian(ylim=c(0,200))

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
#write.csv(ttest_taxon, "ttest_taxon_2023pt2.csv", row.names = T)

#zoom sur certains échantillons
ggplot(all_data[all_data$Sample == "32", ], aes(x = Taxon, y = Nombre, fill=meth)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Sample) +
  labs(title = "Comparaison du nombre de pollens de chaque taxon en fonction méthodes d'ID")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  coord_cartesian(ylim = c(0, 2500))

### Corrélations entre nb polllen et methode
# Corrélation entre les valeurs uniques de TOTAL et meth
all_data <- all_data %>%
  mutate(
    Sample = as.numeric(Sample)
    Taxon = as.numeric(Taxon)
    TOTAL = as.numeric(TOTAL),
    meth = as.numeric(meth)
  )
unique_data <- all_data %>%
  distinct(TOTAL, meth)  # On garde les paires uniques de TOTAL et meth
correlation_total_meth <- cor(unique_data$TOTAL, unique_data$meth, use = "complete.obs")
# Affichage de la corrélation
cat("Corrélation entre les valeurs uniques de TOTAL et meth :", correlation_total_meth, "\n")

# Calcul de la matrice de corrélation
correlation_matrix <- cor(all_data, use = "complete.obs")

# Affichage de la matrice de corrélation
print("Matrice de corrélation entre les variables numériques :")
print(correlation_matrix)

# Optionnel : Visualisation de la matrice de corrélation avec ggplot2 ou corrplot
library(corrplot)
corrplot(correlation_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 45)
