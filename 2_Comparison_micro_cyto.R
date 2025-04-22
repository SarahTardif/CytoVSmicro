# import packages and libraries
library(dplyr)
library(tidyr)
library(ggplot2)


# Data 
# run the script "1.find_threshold_values" to obtain dataMI and dataCY

# Graph for dataCY - total number of pollen grains per sample
ggplot(dataCY, aes(x = Sample, y = Nombre)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Nombre total de pollens pour la cyto")+
  coord_cartesian(ylim = c(0, 1000))
# Graph for dataMI - total number of pollen grains per sample
ggplot(dataMI, aes(x = Sample, y = Nombre)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Nombre total de pollens pour la micro")+
  coord_cartesian(ylim = c(0, 5000))

# group dataframes dataMI and dataCY into one
dataCY$meth<-"cyto"
dataMI$meth<-"micro"
dataCY$Sample <- gsub("CY_", "", dataCY$Sample)
dataMI$Sample <- gsub("MI_", "", dataMI$Sample)
all_data<-full_join(dataCY,dataMI)
all_data$Sample <- as.factor(all_data$Sample)
all_data$Taxon <- as.factor(all_data$Taxon)
all_data$meth <- as.factor(all_data$meth)

# graph of total number of pollen grains per sample between for both ID methods (Microscopy and cytometry+ML)
ggplot(all_data, aes(x = meth, y = TOTAL, fill = meth)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 1500))
summary(all_data$TOTAL[all_data$meth == "micro"])
summary(all_data$TOTAL[all_data$meth == "cyto"])

# graph of number of pollen grains per taxon for each sample and for both ID methods (Microscopy and cytometry+ML) 
# check that dataCY data are grouped in the same categories as dataMI, otherwise they are not comparable!
ggplot(all_data, aes(x = Taxon, y = Nombre, fill=meth)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Sample) +
  labs(title = "Comparaison du nombre de pollens de chaque taxon en fonction méthodes d'ID")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  coord_cartesian(ylim = c(0, 1000))

# graph of number of pollen grains per taxon (samples considered as replicates) for both ID methods (Microscopy and cytometry+ML) 
ggplot(all_data, aes(x = Taxon, y = Nombre, fill=meth)) +
  geom_boxplot() +
  labs(title = "",
       x = "Taxon and method",
       y = "Pollen count")+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.4, size =14),
    axis.text.y = element_text(size =14),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    legend.text=element_text(size=14),
    legend.title = element_blank()
    )+
  coord_cartesian(ylim=c(0,350))

# ttest for each taxon
donnees_cyto <- subset(all_data, meth == "cyto")
donnees_micro <- subset(all_data, meth == "micro")
ttest_taxon <- data.frame(Taxon = character(), p_value = numeric(), difference = numeric())
taxons <- unique(all_data$Taxon)
for (taxon in taxons) {
  donnees_cyto_taxon <- subset(donnees_cyto, Taxon == taxon)
  donnees_micro_taxon <- subset(donnees_micro, Taxon == taxon)
  if (nrow(donnees_cyto_taxon) == 0) {
    donnees_cyto_taxon <- data.frame(Nombre = 0)
  }
  if (nrow(donnees_micro_taxon) == 0) {
    donnees_micro_taxon <- data.frame(Nombre = 0)
  }
  diff_mesures <- mean(donnees_micro_taxon$Nombre, na.rm = TRUE) - mean(donnees_cyto_taxon$Nombre, na.rm = TRUE)
  if (length(donnees_cyto_taxon$Nombre) > 1 & length(donnees_micro_taxon$Nombre) > 1) {
    resultat <- t.test(donnees_cyto_taxon$Nombre, donnees_micro_taxon$Nombre, var.equal = TRUE)
    p_value <- resultat$p.value
  } else {
    p_value <- NA  # If not enough data, we set a missing p-value
  }
  ttest_taxon <- rbind(ttest_taxon, data.frame(Taxon = taxon, p_value = p_value, difference = diff_mesures, t=resultat$statistic, confidence_interval= paste(resultat$conf.int[1],resultat$conf.int[2],sep=";"), df=resultat$parameter))
}
print(ttest_taxon)
getwd()
#write.csv(ttest_taxon, "ttest_taxon_genus.csv",row.names = T)

#zoom on a specific sample
ggplot(all_data[all_data$Sample == "22", ], aes(x = Taxon, y = Nombre, fill=meth)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Sample) +
  labs(title = "Comparaison du nombre de pollens de chaque taxon en fonction méthodes d'ID")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  coord_cartesian(ylim = c(0, 150))




### Correlations
# Correlation between total number of pollen grains and the ID method
all_data <- all_data %>%
  mutate(
    Sample = as.numeric(Sample),
    Taxon = as.numeric(Taxon),
    TOTAL = as.numeric(TOTAL),
    meth = as.numeric(meth)
  )
unique_data <- all_data %>%
  distinct(TOTAL, meth) 
correlation_total_meth <- cor(unique_data$TOTAL, unique_data$meth, use = "complete.obs")
cat("Corrélation entre les valeurs uniques de TOTAL et meth :", correlation_total_meth, "\n")
# Correlation matrix
correlation_matrix <- cor(all_data, use = "complete.obs")
print("Matrice de corrélation entre les variables numériques :")
print(correlation_matrix)
library(corrplot)
library(RColorBrewer)
corrplot(correlation_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45, col = rev(brewer.pal(n = 11, name = "RdBu")))

