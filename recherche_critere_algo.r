library(dplyr)
library(here)
library(readr)
library(purrr)
library(fs)
library(tidyr)
library(ggplot2)
library(tibble)

rm(list=ls())





## Crée un vecteur des noms de fichiers, avec tout le chemin d'accès
##vérifier qu'on est dans le bon repertoire
setwd("./Ech_2140_wodebris_ID")
dir_list <- list.files(here(getwd()),
                       pattern = "^ID.*\\.csv$", full.names = TRUE)

## Nomme le vecteur avec seulement le nom de fichier, sans l'extension
names(dir_list) <- path_ext_remove(list.files(here(getwd()),
                       pattern = "^ID.*\\.csv$", full.names = TRUE))
names(dir_list) <-path_ext_remove(basename(dir_list))

files_df <- map_dfr(dir_list, read_csv, .id = "Sample_Name") ## combine tous les fichiers csv en un, ajoute une colonne Sample_name avec le nom de l'échantillon
dataCY$Sample <- gsub("ID_CY ", "", dataCY$Sample)
write.csv(files_df,"../Ech2140_wodebris_ID.csv",row.names=FALSE)





#data
data<-read.csv("ID_Ech2140_balanced_genus.csv",sep=",",h=T)
data$Sample <- gsub("ID_CY ", "", data$Sample)
data<-data[data$Sample!="21",]
data<-data[data$Sample!="32",]
subdata<-data[,-!names(data) %in% "Sample"]
sp_max <- apply(subdata, 1, function(row) names(subdata)[which.max(row)])
prob_max <- apply(subdata, 1, function(row) max(row))
datasp <- data.frame(species = sp_max, prob = prob_max, Sample=data$Sample)

#regarder distribution proba max pollens 
summary(datasp$prob)  # Résumé (min, 1er quartile, médiane, moyenne, 3e quartile, max)

aggregate(prob ~ species, data = datasp, summary)
ggplot(datasp, aes(x = species, y = prob)) +
  geom_boxplot(aes(fill = species), alpha = 0.7, show.legend=FALSE) +
  labs(title = "Distribution de prob par espèce", 
       x = "Espèce", 
       y = "Probabilité") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#enlever tous les pollens dont la proba est inferieure à 0.5
dataCY<- datasp[datasp$prob >= 0.75, ]
dataCY <- dataCY[, -which(names(dataCY) %in% c("prob"))]
dataCY<-table(dataCY$Sample, dataCY$species)
dataCY<-as.data.frame.matrix(dataCY)
#group<-data80CY[,-1]
dataCY$TOTAL<-rowSums(dataCY)
dataCY<-rownames_to_column(dataCY, var = "Sample")

#write.csv(dataCY,"C:/Users/sarah/Downloads/nb_pollen_taxon_sample_genus.csv",row.names=F)

dataCY_long <- reshape(
  dataCY,
  varying = list(names(dataCY)[2:ncol(dataCY)]),
  v.names = "Count",
  timevar = "Genus",
  times = names(dataCY)[2:ncol(dataCY)],
  direction = "long"
)

ggplot(dataCY_long, aes(x = Genus, y = Count)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 500))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




### Trouver meilleur seuil de proba (sans faire arbritairement)
df <- datasp
# Initialisation des seuils de probabilité
seuils <- seq(0, 1, by = 0.01) # seuils de 0 à 1 avec un pas de 0.01
# Calcul de la moyenne et de l'écart-type du nombre de lignes par Sample pour chaque seuil
resultats <- do.call(rbind, lapply(seuils, function(x) {
  df_filtered <- df[df$prob >= x, ] # Filtrer les lignes selon le seuil
  lignes_par_sample <- tapply(df_filtered$species, df_filtered$Sample, length) # Nombre de lignes par Sample
  moyenne <- mean(lignes_par_sample, na.rm = TRUE) # Moyenne des lignes par Sample
  ecart_type <- sd(lignes_par_sample, na.rm = TRUE) # Écart-type des lignes par Sample
  data.frame(seuil = x, moyenne = moyenne, ecart_type = ecart_type)
}))
# Tracé du graphique avec barres d'erreur
ggplot(resultats, aes(x = seuil, y = moyenne)) +
  geom_line(color = "black") +
  geom_errorbar(aes(ymin = moyenne - ecart_type, ymax = moyenne + ecart_type), 
                width = 0.01, color = "grey", alpha = 0.6) +
  geom_hline(yintercept=275)+
  labs(
    title = "recherche seuil de probabilité ",
    x = "probabilité",
    y = "Nombre moyen de pollen total par échantillon"
  ) +
  theme_minimal()

valseuil<-resultats$seuil[which.min(abs(resultats$moyenne - 275))]
valseuil
### prédire par groupe quand peu de dissimilarité entre les especes
data<-read.csv("Ech2140_ID_prob.csv",sep=",",h=T)
data<-data[data$Sample_Name!="21",]
subdata<-data[,-!names(data) %in% "Sample_Name"]
data$prob_max <- colnames(subdata)[apply(subdata, 1, which.max)]
datasansdebris<-data[-which(data$prob_max=="Debris"),]
datasansdebris$prob_max<-NULL
# Transposer la matrice et calculer les distances
matrice <- t(datasansdebris)
distances <- dist(matrice, method = "euclidean")
print(as.matrix(distances))
# Clustering hiérarchique
hclust_resultat <- hclust(distances, method = "ward.D2")
# Visualiser le dendrogramme
plot(hclust_resultat)
# Découper en 2 groupes
groupes_classes <- cutree(hclust_resultat, k = 2)
print(groupes_classes)
# Mapping des prédictions vers les groupes
predictions <- c("cat1", "cat3", "cat5")
groupes_predits <- sapply(predictions, function(classe) {
  groupe <- groupes_classes[classe]
  return(paste0("groupe_", groupe))
})
data.frame(predictions = predictions, groupes = groupes_predits)

### regarder c'est quoi la proportion de debris vs pollens dans echantillons
data<-read.csv("ID_Ech2140_wodebris.csv",sep=",",h=T)
data<-data[data$Sample!="21",]
subdata<-data[,-!names(data) %in% "Sample"]
data$prob_max <- colnames(subdata)[apply(subdata, 1, which.max)]
data_probmax<-data[,c("Sample_Name","prob_max")]
data_probmax$prob_max <- ifelse(data_probmax$prob_max != "Debris", "pollen", data_probmax$prob_max)
data_probmax<-table(data_probmax)
data_probmax<-as.data.frame.matrix(data_probmax)
data_probmax <- rownames_to_column(data_probmax, var = "sample")
data_probmax$diff<- data_probmax$Debris - data_probmax$pollen
write.csv(data_probmax,"C:/Users/sarah/Downloads/data_probmax.csv",row.names=FALSE)
data_probmax_long <- data_probmax %>%
  pivot_longer(cols = c("Debris", "pollen"), 
               names_to = "type", 
               values_to = "value")
ggplot(data_probmax_long,aes(x=type,y=value)) +
  geom_boxplot() +
  labs(title = "pollens vs debris dans echantillons avec prob max") +
  theme_minimal()
summary(data_probmax)

### regarder probabilité max vs proba debris
data<-read.csv("ID_Ech2140_wodebris.csv",sep=",",h=T)
subdata<-data[,-!names(data) %in% "Sample"]
sp_max <- apply(subdata, 1, function(row) names(subdata)[which.max(row)])
prob_max <- apply(subdata, 1, function(row) max(row))
datasp <- data.frame(species = sp_max, prob = prob_max, sample=data$Sample)
datasp$debris_prob <- subdata$Debris
datasansdebris<- datasp[datasp$species != "Debris", ]
datasansdebris$diff<- datasansdebris$prob - datasansdebris$debris_prob
datasansdebris_long <- datasansdebris %>%
  pivot_longer(cols = c("prob", "debris_prob"), 
               names_to = "type", 
               values_to = "value")

ggplot(datasansdebris_long, aes(x = type, y = value, fill = type)) +
  geom_boxplot() +
  labs(title = "Boxplot de prob et debris_prob",
       y = "Valeur", x = "Espèce", fill = "Type") +
  facet_wrap(~ sample) +
  theme_minimal()

ggplot(datasansdebris, aes(x = species, y = diff, fill = species)) +
  geom_boxplot(alpha = 0.7, outlier.color = "black", outlier.size = 2, show.legend=FALSE) +
  theme_minimal() +
  labs(title = "Distribution des différences par espèce",
       x = "Espèce",
       y = "Différence") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

result <- datasansdebris %>%
  filter(!is.na(diff)) %>%  # Supprime les NA
  group_by(species) %>%
  summarise(mean_difference = mean(diff))
