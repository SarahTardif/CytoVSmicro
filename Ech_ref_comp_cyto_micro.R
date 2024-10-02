rm(list=ls())
remove(list = ls(.rs.CachedDataEnv), envir = .rs.CachedDataEnv)
# installer et charger les packages necessaires
library(dplyr)
library(tidyr)
#install.packages("ggplot2")

library(ggplot2)

# choisir le répertoire de travail
setwd("C:/Users/sarah/OneDrive - UQAM/PhD/0_Papier_Méthodologie/Analyses_comp_methodes")
getwd()

# récupérer les fichiers avec les données / pollens identifiés
data<- read.csv("Ech_ref_cyto_micro.csv", sep=";", h=T)
data$Species<-data$X
data$X<-NULL

# faire graphique comparant les ID micro et cyto

data_plot <- pivot_longer(data, cols = -c(Species,sd_micro), names_to = "Methode", values_to = "Nombre")
data_plot[is.na(data_plot)] <- 0

data_plot <- data_plot %>%
  mutate(ymin = ifelse(Methode == "micro", Nombre - sd_micro, NA),
         ymax = ifelse(Methode == "micro", Nombre + sd_micro, NA))
ggplot(data_plot, aes(x = Species, y = Nombre, fill = Methode)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 0.9), width = 0.25) +
  coord_cartesian(ylim = c(0, 9500)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Calculer les bornes inférieure et supérieure des intervalles
data$lower_bound <- data$micro - data$sd_micro
data$upper_bound <- data$micro + data$sd_micro

# Vérifier si les valeurs cyto se trouvent dans les intervalles
data$overlap <- with(data, cyto >= lower_bound & cyto <= upper_bound)

# Afficher les résultats
print(data)

write.csv(data,"Ech_ref_overlap.csv",row.names = F) #conserver le fichier avec toutes les données propres prêt à être analysés


# test de student
results <- data.frame(especes = character(), p_value = numeric(), diff=numeric(), stringsAsFactors = FALSE)

for (i in 1:nrow(data)) {
  chi2_result <- chisq.test(c(data[i, "micro"], data[i, "cyto"]))
  species<-data[i, "Species"]
  p_value <- chi2_result$p.value
  diffs <- data[i, "micro"] - data[i, "cyto"]
  results <- rbind(results, data.frame(especes = species, p_value = p_value, diff=diffs))
}
print(results)
write.csv(results,"Ech_ref_chi2_cyto_micro.csv",row.names = F) #conserver le fichier avec toutes les données propres prêt à être analysés