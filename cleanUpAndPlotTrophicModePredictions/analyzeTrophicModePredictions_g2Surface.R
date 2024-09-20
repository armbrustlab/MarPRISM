library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/g2Surface/")

trop <- read_csv("G2_surface_trophicModePredictions_updatedMarferret_marmicroDb2023_noOutliers_fall2023")
meta <- read_csv("G2_surface_tpm_updatedMarferret_marmicroDb2023_sampleTaxa_noOutliers_fixedTPM_fall2023.csv")
prob <- read_csv("G2_trophicModePredictionsProbabilities_surface_updatedMarferret_marmicroDb2023_noOutliers_fall2023")

merged <- cbind(meta, trop, prob)

#Mix, Het, Phot
#1, 0, 2,
merged <- merged %>% mutate(xg_pred = ifelse(xg_pred == 1, "Mix", ifelse(xg_pred == 0, "Het", "Phot")))

colnames(merged)[2]
colnames(merged)[2] <- "taxa"

merged <- merged %>% mutate(xg_pred = str_c(xg_pred, "otrophic"))
merged <- merged %>% mutate(xg_pred = str_replace(xg_pred, "Hetotrophic", "Heterotrophic"))
merged %>% distinct(xg_pred)

colnames(merged)[1] <- "sample"

merged <- merged %>% mutate(Size = ifelse(str_detect(sample, "2um"), "0.2um", NA))
merged <- merged %>% mutate(Size = ifelse(str_detect(sample, "3um"), "3um", Size))
merged %>% filter(is.na(Size))

merged <- merged %>% mutate(station = str_extract(sample, "S[0-9]{1,}C[0-9]{1,}"))

sample <- read_excel("Gradients2_discrete_samples.xlsx", sheet = 2)

head(sample)
colnames(sample) <- sample[1,]
sample <- sample[-1,]
head(sample)

sample <- sample %>% select(Date:Notes)

sample <- sample %>% as.data.frame()

merged <- merged %>% mutate(Station = parse_number(station))

sample <- sample %>% mutate(Station = parse_number(Station))

merged <- merged %>% mutate(Cast = 1)

sample$Station <- as.numeric(sample$Station)
sample$Cast <- as.numeric(sample$Cast)

sample <- sample %>% distinct(Station, Cast, `Depth (m)`, `Time Start (HST)`, Latitude, Size)

sample$`Time Start (HST)` <- as.numeric(sample$`Time Start (HST)`)

sample <- sample %>% distinct(Station, Cast, `Depth (m)`, `Time Start (HST)`, Latitude, Size)

sample %>% group_by(Station, Cast, Size) %>% summarize(n = n()) %>% filter(n > 1)

merged %>% anti_join(sample, by = c("Station", "Cast", "Size"))

nrow(merged)
merged <- merged %>% left_join(sample, by = c("Station", "Cast", "Size"))
nrow(merged)

merged$Latitude <- as.numeric(merged$Latitude)

#exclude nonspecies and nonprotists
merged <- merged %>% filter(!taxa %in% c("cellular organisms", "Micromonas <green algae>", 
                 "Calanus finmarchicus", "Eucyclops serrulatus", 
                 "Hydra vulgaris", "Oikopleura dioica", "Paracyclopina nana", 
                 "unclassified Acanthoeca", "unclassified Micromonas", 
                 "unclassified Phaeocystis", "unclassified Chrysochromulina") & 
                  str_detect(taxa, " "))

totalSample <- merged %>% group_by(taxa, station) %>% summarize(total = n()) 

exclude <- merged %>% group_by(taxa, xg_pred, station) %>% summarize(numSamples = n()) %>% 
  left_join(totalSample, by = c("taxa", "station")) %>% mutate(prop = numSamples/total*100) %>% 
  filter(xg_pred %in% c("Heterotrophic", "Phototrophic")) %>% ungroup() %>% 
  filter(prop > 25) %>%
  group_by(taxa, station) %>% distinct(xg_pred) %>% summarize(n = n()) %>% 
  filter(n > 1) %>% ungroup() %>% distinct(taxa, station)

merged %>% semi_join(exclude, by = c("station", "taxa")) %>% nrow()

merged %>% semi_join(exclude, by = c("station", "taxa")) %>% write_csv("excluded.csv")

merged <- merged %>% anti_join(exclude, by = c("station", "taxa"))

#write cleaned up trophic mode predictions to csv file
merged %>% write_csv("G2_surface_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")
