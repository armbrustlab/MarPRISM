library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/g1Surface/")

trop <- read_csv("G1_trophicModePredictions_surface_updatedMarferret_marmicroDb2023_noOutliers_fall2023")
meta <- read_csv("G1_surface_tpm_updatedMarferret_marmicroDb2023_sampleTaxa_noOutliers_fixedTPM_fall2023.csv")
prob <- read_csv("G1_trophicModePredictionsProbabilities_surface_updatedMarferret_marmicroDb2023_noOutliers_fall2023")

merged <- cbind(meta, trop, prob)

#Mix, Het, Phot
#1, 0, 2,
merged <- merged %>% mutate(xg_pred = ifelse(xg_pred == 1, "Mix", ifelse(xg_pred == 0, "Het", "Phot")))

colnames(merged)[1]
colnames(merged)[1] <- "sample"

merged <- merged %>% mutate(Size = ifelse(str_detect(sample, "0.2um"), "0.2um", NA))
merged <- merged %>% mutate(Size = ifelse(str_detect(sample, "3um"), "3um", Size))
merged %>% filter(is.na(Size))

sample <- read_csv("ctd_g1_g2_15m_source_Stephen.csv")
sample <- sample %>% select(SAMPLE_ID, DEPTH, LATITUDE, DATETIME)
sample <- sample %>% as.data.frame()
sample <- sample %>% mutate(SAMPLE_ID = str_replace(SAMPLE_ID, "0_2um", "0.2um"))

time <- sample %>% semi_join(merged, by = c("SAMPLE_ID" = "sample"))

time <- time %>% mutate(time = str_extract(DATETIME, "[0-9]{1,}:[0-9]{1,}"))

time %>% distinct(time) %>% arrange(time)

sample <- sample %>% select(-DATETIME)

merged %>% anti_join(sample, by = c("sample" = "SAMPLE_ID"))

nrow(merged)
merged <- merged %>% left_join(sample, by = c("sample" = "SAMPLE_ID"))
nrow(merged)

merged <- merged %>% mutate(xg_pred = str_c(xg_pred, "otrophic"))
merged <- merged %>% mutate(xg_pred = str_replace(xg_pred, "Hetotrophic", "Heterotrophic"))

merged <- merged %>% mutate(station = str_extract(sample, "S[0-9]{1,}C*[0-9]*"))
merged %>% filter(is.na(station))

#exclude non species and non protists
merged <- merged %>% filter(!tax_name %in% c("Calanus finmarchicus", "Hydra vulgaris", "Neocalanus flemingeri", 
                   "Oikopleura dioica", "Paracyclopina nana", 
                   "Micromonas <green algae>", "unclassified Acanthoeca", 
                   "unclassified Chrysochromulina", "unclassified Micromonas", 
                   "unclassified Phaeocystis", "cellular organisms", 
                   "Pleuromamma xiphias") & 
    str_detect(tax_name, " "))

totalSample <- merged %>% group_by(tax_name, station) %>% summarize(total = n()) 

exclude <- merged %>% group_by(tax_name, xg_pred, station) %>% summarize(numSamples = n()) %>% 
  left_join(totalSample, by = c("tax_name", "station")) %>% mutate(prop = numSamples/total*100) %>% 
  filter(xg_pred %in% c("Heterotrophic", "Phototrophic")) %>% ungroup() %>% 
  filter(prop > 25) %>%
  group_by(tax_name, station) %>% distinct(xg_pred) %>% summarize(n = n()) %>% 
  filter(n > 1) %>% ungroup() %>% distinct(tax_name, station)

merged %>% semi_join(exclude, by = c("station", "tax_name")) %>% nrow()

merged %>% semi_join(exclude, by = c("station", "tax_name")) %>% write_csv("excluded.csv")

merged <- merged %>% anti_join(exclude, by = c("station", "tax_name"))

#write cleaned up trophic predictions to csv file
merged %>% write_csv("G1_surface_trophicPredictions_marFerret2023_cleanedUp_noOutliers.csv")
