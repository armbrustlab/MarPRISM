library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/g3")

trop <- read_csv("G3_trophicModePredictions_updatedMarferret_marmicroDb2023_noOutliers_fall2023")
meta <- read_csv("G3_surface_tpm_updatedMarferret_marmicroDb2023_sampleTaxa_noOutliers_fixedTPM_fall2023.csv")
prob <- read_csv("G3_trophicModePredictionsProbabilities_updatedMarferret_marmicroDb2023_noOutliers_fall2023")

merged <- cbind(meta, trop, prob)

#Mix, Het, Phot
#1, 0, 2,
merged <- merged %>% mutate(xg_pred = ifelse(xg_pred == 1, "Mix", ifelse(xg_pred == 0, "Het", "Phot")))

colnames(merged)[2]
colnames(merged)[2] <- "taxa"

merged <- merged %>% mutate(xg_pred = str_c(xg_pred, "otrophic"))
merged <- merged %>% mutate(xg_pred = str_replace(xg_pred, "Hetotrophic", "Heterotrophic"))

#downloaded from https://docs.google.com/spreadsheets/d/1lKwvlNjUV8dXHKO5HkMNa2m70QQYydrTf-Fb6PRQ9oE/edit#gid=0
sample <- read_csv("../metaT sample log - metaT.csv")

sample <- sample %>% filter(Cruise == "Gradients3" | Cruise == "Gradients 3")

sample <- sample %>% filter(Type == "polyA")

colnames(merged)[1] <- "sample"

merged <- merged %>% mutate(Sample.ID = str_extract(sample, "UW[0-9]{1,}"))
merged <- merged %>% mutate(Sample.ID = ifelse(is.na(Sample.ID), sample, Sample.ID))
merged %>% filter(is.na(Sample.ID))

sample <- sample %>% mutate(Sample.ID = ifelse(Exp == "Diel", str_c("S",Station,"C",as.character(Cast)), Sample.ID))

sample <- sample %>% select(Cruise:Exp)

merged %>% anti_join(sample, by = c("Sample.ID"))

nrow(merged)
merged <- merged %>% left_join(sample %>% select(-c(Date, Time, Date.Time, Replicate)) %>% distinct(), by = c("Sample.ID"))
nrow(merged)

#exclude nonspecies and nonprotists 
merged <- merged %>% filter(!taxa %in% c("cellular organisms", "Micromonas <green algae>", 
                 "Calanus finmarchicus", "Eucyclops serrulatus", 
                 "Hydra vulgaris", "Oikopleura dioica", "Paracyclopina nana", 
                 "unclassified Acanthoeca", "unclassified Micromonas", 
                 "unclassified Phaeocystis", "unclassified Chrysochromulina", 
                 "Homarus americanus", "Neocalanus flemingeri") & 
                   str_detect(taxa, " "))


totalSample <- merged %>% group_by(taxa, Latitude) %>% summarize(total = n()) 

exclude <- merged %>% group_by(taxa, xg_pred, Latitude) %>% summarize(numSamples = n()) %>% 
  left_join(totalSample, by = c("taxa", "Latitude")) %>% mutate(prop = numSamples/total*100) %>% 
  filter(xg_pred %in% c("Heterotrophic", "Phototrophic")) %>% ungroup() %>% 
  filter(prop > 25) %>%
  group_by(taxa, Latitude) %>% distinct(xg_pred) %>% summarize(n = n()) %>% 
  filter(n > 1) %>% ungroup() %>% distinct(taxa, Latitude)

merged %>% semi_join(exclude, by = c("Latitude", "taxa")) %>% nrow()

merged %>% semi_join(exclude, by = c("Latitude", "taxa")) %>% 
  write_csv("excluded.csv")

merged <- merged %>% anti_join(exclude, by = c("Latitude", "taxa"))

#write cleaned up trophic mode predictions to csv file
merged %>% write_csv("G3_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")
