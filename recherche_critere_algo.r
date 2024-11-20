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
setwd("./Ech_2140_ID_prob")
dir_list <- list.files(here(getwd()),
                       pattern = "^ID.*\\.csv$", full.names = TRUE)

## Nomme le vecteur avec seulement le nom de fichier, sans l'extension
names(dir_list) <- path_ext_remove(list.files(here(getwd()),
                       pattern = "^ID.*\\.csv$", full.names = TRUE))
names(dir_list) <-path_ext_remove(basename(dir_list))

files_df <- map_dfr(dir_list, read_csv, .id = "Sample_Name") ## combine tous les fichiers csv en un, ajoute une colonne Sample_name avec le nom de l'échantillon
files_df$Sample_Name <- gsub("ID_CY ", "", files_df$Sample_Name)
write.csv(files_df,"../Ech2140_ID_prob.csv",row.names=FALSE)

## regarder c'est quoi la proportion de debris vs pollens dans echantillons
data<-read.csv("Ech2140_ID_prob.csv",sep=",",h=T)
data<-data[data$Sample_Name!="21",]
subdata<-data[,-!names(data) %in% "Sample_Name"]
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
## regarder probabilité max vs proba debris
data<-read.csv("Ech2140_ID_prob.csv",sep=",",h=T)
subdata<-data[,-!names(data) %in% "Sample_Name"]
sp_max <- apply(subdata, 1, function(row) names(subdata)[which.max(row)])
prob_max <- apply(subdata, 1, function(row) max(row))
datasp <- data.frame(species = sp_max, prob = prob_max, sample=data$Sample_Name)
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


#regarder distribution proba max pollens 
summary(datasansdebris$prob)  # Résumé (min, 1er quartile, médiane, moyenne, 3e quartile, max)

aggregate(prob ~ species, data = datasansdebris, summary)
ggplot(datasansdebris, aes(x = species, y = prob)) +
  geom_boxplot(aes(fill = species), alpha = 0.7, show.legend=FALSE) +
  labs(title = "Distribution de prob par espèce", 
       x = "Espèce", 
       y = "Probabilité") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#enlever tous les pollens dont la proba est inferieure à 0.5
data50<- datasansdebris[datasansdebris$prob >= 0.5, ]
data50 <- data50[, -which(names(data50) %in% c("prob", "debris_prob"))]
data50CY<-table(data50$sample, data50$species)
data50CY<-as.data.frame.matrix(data50CY)
#group<-data50CY[,-1]
data50CY$TOTAL<-rowSums(data50CY)
data50CY<-rownames_to_column(data50CY, var = "Sample")

data50CY_long <- reshape(
  data50CY,
  varying = list(names(data50CY)[2:ncol(data50CY)]),
  v.names = "Count",
  timevar = "Species",
  times = names(data50CY)[2:ncol(data50CY)],
  direction = "long"
)

ggplot(data50CY_long, aes(x = Species, y = Count)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 3000))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
