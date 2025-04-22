# import packages and libraries
library(dplyr)
library(here)
library(readr)
library(purrr)
library(fs)
library(tidyr)
library(ggplot2)
library(tibble)


# import data - pollen ID with microscopy
dataMI<- read.csv("ID_2023_2140_MI.csv", sep=";", h=T)
names(dataMI) <- gsub("Total\\.count\\.", "", names(dataMI))
group<-subset(dataMI,select=-Sample)
dataMI$TOTAL<-rowSums(group) # total number of pollen in each sample
dataMI <- pivot_longer(dataMI, cols = -c(Sample, TOTAL), names_to = "Taxon", values_to = "Nombre")
dataMI[is.na(dataMI)] <- 0
numeros1 <- as.numeric(gsub("MI_", "", dataMI$Sample))
dataMI <- dataMI[order(numeros1), ]
dataMI$Sample <- factor(dataMI$Sample, levels = unique(dataMI$Sample))
dataMI<-dataMI[dataMI$Sample!=21,] # remove sample 21 - Lab manipulations didn't went well
dataMI<-dataMI[dataMI$Sample!=32,] # remove sample 32 - Eleagnus tree just above the sampler 

# import data - pollen ID with Cytometry + RandomForest 
data<-read.csv("ID_Ech2140_balanced_genus.csv",sep=",",h=T)
data$Sample <- gsub("ID_CY ", "", data$Sample)
data<-data[data$Sample!="21",] # remove sample 21 - Lab manipulations didn't went well
data<-data[data$Sample!="32",] # remove sample 32 - Eleagnus tree just above the sampler 
subdata<-data[,-!names(data) %in% "Sample"]
sp_max <- apply(subdata, 1, function(row) names(subdata)[which.max(row)]) # list of species (or genus) associated with the highest classification probability
prob_max <- apply(subdata, 1, function(row) max(row)) # list of the highest classification probability for each pollen grain
datasp <- data.frame(species = sp_max, prob = prob_max, Sample=data$Sample) # dataframe with for each pollen grain (=each row) the species (or genus), the classification probability and the sample 

# distribution of classification probabilities according to species (or genus)
summary(datasp$prob)  
aggregate(prob ~ species, data = datasp, summary)
ggplot(datasp, aes(x = species, y = prob)) +
  geom_boxplot(aes(fill = species), alpha = 0.7, show.legend=FALSE) +
  labs(title = "Distribution de prob par espèce", 
       x = "Espèce", 
       y = "Probabilité") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##### Script to find a probability treshold per taxon ######

## data cleaning to have the same taxa for dataMI and dataCY
datasp$species <- ifelse(datasp$species == "Malus"|datasp$species =="Prunus"|datasp$species =="Pyrus"|datasp$species =="Sorbus"|datasp$species =="Amelanchier", "Rosaceae", datasp$species)
datasp$species <- ifelse(datasp$species == "Juniperus", "Cupressaceae", datasp$species)
datasp$species <- ifelse(datasp$species == "Pinus"|datasp$species =="Picea", "Pinaceae", datasp$species)
datasp$species <- ifelse(datasp$species == "Corylus"|datasp$species =="Ostrya", "Corylus.Ostrya", datasp$species)
datasp$species <- ifelse(datasp$species == "Carpinus"|datasp$species =="Celtis"|datasp$species =="Robinia", "NI.Others", datasp$species)

dataMI$Taxon <- ifelse(dataMI$Taxon == "Larix", "Pinaceae", dataMI$Taxon)
dataMI$Taxon <- ifelse(dataMI$Taxon == "Olaeceae", "Syringa", dataMI$Taxon)
dataMI$Taxon <- ifelse(dataMI$Taxon == "Armoise"|dataMI$Taxon =="Autres.herbes"|dataMI$Taxon =="Gramineae"|dataMI$Taxon =="Plantago"|dataMI$Taxon =="Typha"|dataMI$Taxon == "NI", "NI.Others", dataMI$Taxon)
dataMI <- aggregate(cbind(Nombre)~Sample+Taxon+TOTAL, data = dataMI, FUN = function(x) sum(x, na.rm = TRUE))

meanMI<-list()
valseuil<-list()

for (taxon in unique(datasp$species)){
  df <- datasp[datasp$species==taxon,]
  MItaxon<-dataMI[dataMI$Taxon==taxon,]
  meanMItaxon<- mean(MItaxon$Nombre) # mean number of pollen grains per sample and for the specific taxon
  meanMI[[taxon]] <- meanMItaxon
seuils <- seq(0, 1, by = 0.01) # list of all possible treshold values (= all values between 0 and 1 with a 0.01 step)
# Calculation of the mean and standard deviation of the number of pollen grains per Sample for each possible threshold value
resultats <- do.call(rbind, lapply(seuils, function(x) {
  df_filtered <- df[df$prob >= x, ] 
  lignes_par_sample <- tapply(df_filtered$species, df_filtered$Sample, length) 
  moyenne <- mean(lignes_par_sample, na.rm = TRUE) 
  ecart_type <- sd(lignes_par_sample, na.rm = TRUE) 
  data.frame(seuil = x, moyenne = moyenne, ecart_type = ecart_type)  
}))
# Graph with error bars
plot<-ggplot(resultats, aes(x = seuil, y = moyenne)) +
      geom_line(color = "black") +
      geom_errorbar(aes(ymin = moyenne - ecart_type, ymax = moyenne + ecart_type), 
                width = 0.01, color = "grey", alpha = 0.6) +
      geom_hline(yintercept=meanMItaxon)+
      labs(
        title = paste("recherche seuil de probabilité ",taxon),
        x = "probabilité",
        y = "Nombre moyen de pollen total par échantillon"
      ) +
      theme_minimal()+
      theme(axis.text = element_text(size = 20))
#print(plot)
# dataframe with threshold values for each taxon (= classification probability to be considered so that the difference in pollen count between the two identification methods is minimal)
valseuil[[taxon]]<-resultats$seuil[which.min(abs(resultats$moyenne - meanMItaxon))] 
}
#write.csv(valseuil,"./valseuil.csv",row.names=F)

# avoir lles comptes de pollen de chaque taxon avoir dataCY mais avec seuil de proba différent spécifique à chaque taxon
valseuil<-as.data.frame(valseuil)
valseuil <- pivot_longer(valseuil, cols=everything(),names_to = "taxon", values_to = "probaseuil")
dataCY<-data.frame()
for(taxon in unique(datasp$species)){
  datataxon<-datasp[datasp$species==taxon,]
  proba <- valseuil[valseuil$taxon == taxon, "probaseuil"][[1]]
  if(proba>0){
    datataxon<-datataxon[datataxon$prob >= proba, ]
    dataCY<- rbind(dataCY,datataxon)
  }
  else {
    datataxon<-datataxon[datataxon$prob >= 0, ]
    dataCY<- rbind(dataCY,datataxon)
  }
}
dataCY <- dataCY[, -which(names(dataCY) %in% c("prob"))]
dataCY<-table(dataCY$Sample, dataCY$species)
dataCY<-as.data.frame.matrix(dataCY)
#group<-data80CY[,-1]
dataCY$TOTAL<-rowSums(dataCY)
dataCY<-rownames_to_column(dataCY, var = "Sample")
Taxa_wanted <- c("Fagus")
for (col in Taxa_wanted) {
  if (!col %in% colnames(dataCY)) {
    dataCY[[col]] <- 0
  }
}
dataCY <- pivot_longer(dataCY, cols = -c(Sample, TOTAL), names_to = "Taxon", values_to = "Nombre")
dataCY$Sample <- factor(dataCY$Sample, levels = unique(dataCY$Sample))

