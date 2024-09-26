library(tidyverse)

setwd("~/Dropbox/grad/research/g1Surface/")

##marprism

trop <- read_csv("G1_trophicModePredictions_surface_updatedMarferret_marmicroDb2023_noOutliers_fall2023")
meta <- read_csv("G1_surface_tpm_updatedMarferret_marmicroDb2023_sampleTaxa_noOutliers_fixedTPM_fall2023.csv")

merged <- cbind(meta, trop)

res <- read_csv("bootstrapPredictions.csv")

merged <- cbind(merged, res)

merged <- merged %>% select(-`...1`)

#Mix, Het, Phot
#1, 0, 2,
merged <- merged %>% mutate(xg_pred = ifelse(xg_pred == 1, "Mix", ifelse(xg_pred == 0, "Het", "Phot")))

merged %>% distinct(xg_pred)
merged <- merged %>% mutate(xg_pred = str_c(xg_pred, "otrophic"))
merged <- merged %>% mutate(xg_pred = str_replace(xg_pred, "Hetotrophic", "Heterotrophic"))
merged %>% distinct(xg_pred)

merged <- merged %>% filter(str_detect(tax_name, " "))
merged <- merged %>% filter(!(tax_name %in% c("Calanus finmarchicus", "Hydra vulgaris", 
                                              "Neocalanus flemingeri", "Oikopleura dioica", "Paracyclopina nana")))
merged <- merged %>% filter(!(tax_name %in% c("cellular organisms", "Micromonas <green algae>")))
merged <- merged %>% filter(!(tax_name %in% c("unclassified Acanthoeca", "unclassified Micromonas", "unclassified Phaeocystis")))
merged <- merged %>% filter(tax_name != "unclassified Chrysochromulina")
merged <- merged%>% filter(!(tax_name %in% c("Calanus finmarchicus", "Eucyclops serrulatus", 
                                             "Hydra vulgaris", "Oikopleura dioica", "Paracyclopina nana")))
merged <- merged %>% filter(tax_name != "Pleuromamma xiphias")

merged %>% distinct(tax_name)

sample <- read_csv("ctd_g1_g2_15m_source_Stephen.csv")

head(sample)
sample <- sample %>% select(SAMPLE_ID, DEPTH, LATITUDE)

sample <- sample %>% as.data.frame()

head(merged)
head(sample)

sample <- sample %>% mutate(SAMPLE_ID = str_replace(SAMPLE_ID, "0_2um", "0.2um"))

merged %>% anti_join(sample, by = c("sample_id" = "SAMPLE_ID"))

nrow(merged)
merged <- merged %>% left_join(sample, by = c("sample_id" = "SAMPLE_ID"))
nrow(merged)

merged <- merged %>% gather(col0:col29, key = "bootNum", value = "boot_xg_pred")

merged <- merged %>% mutate(boot_xg_pred = ifelse(boot_xg_pred == 1, "Mix", ifelse(boot_xg_pred == 0, "Het", "Phot")))

inc <- read_csv("corePfamSummary_boot.csv")

head(merged)
head(inc)

merged <- merged %>% semi_join(inc, by = c("sample_id" = "var", "tax_name", "bootNum"))

merged %>% distinct(boot_xg_pred)
merged <- merged %>% mutate(boot_xg_pred = str_c(boot_xg_pred, "otrophic"))
merged <- merged %>% mutate(boot_xg_pred = str_replace(boot_xg_pred, "Hetotrophic", "Heterotrophic"))
merged %>% distinct(boot_xg_pred)

head(merged)

exclude <- merged %>% filter(boot_xg_pred %in% c("Phototrophic", "Heterotrophic")) %>% group_by(LATITUDE, tax_name, bootNum) %>% 
  distinct(boot_xg_pred) %>% summarize(n = n()) %>% filter(n == 2)
exclude <- exclude %>% ungroup() 
exclude

merged <- merged %>% anti_join(exclude, by = c("LATITUDE", "tax_name", "bootNum"))

#get rid of predictions for transcriptomes with less than 70% of core Pfams expressed

#689/8121 (8.48%)
head(merged)
100*(merged %>% filter(xg_pred != boot_xg_pred) %>% nrow())/(merged %>% nrow())
(merged %>% filter(xg_pred != boot_xg_pred) %>% nrow())
(merged %>% nrow())

merged %>% group_by(xg_pred, boot_xg_pred) %>% summarize(n = n()) %>% 
  mutate(prop = n/(merged %>% nrow())*100)


##ben's model 

trop <- read_csv("G1_trophicModePredictions_surface_updatedMarferret_marmicroDb2023_noOutliers_fall2023_lambert")
meta <- read_csv("G1_surface_tpm_updatedMarferret_marmicroDb2023_sampleTaxa_noOutliers_fixedTPM_fall2023.csv")

merged <- cbind(meta, trop)

pfam <- read_csv("bootstrapPfams_lambert.csv")
res <- read_csv("bootstrapPredictions_lambert.csv")

merged <- cbind(merged, res)

merged <- merged %>% select(-`...1`)

#Mix, Het, Phot
#1, 0, 2,
merged <- merged %>% mutate(xg_pred = ifelse(xg_pred == 1, "Mix", ifelse(xg_pred == 0, "Het", "Phot")))

merged %>% distinct(xg_pred)
merged <- merged %>% mutate(xg_pred = str_c(xg_pred, "otrophic"))
merged <- merged %>% mutate(xg_pred = str_replace(xg_pred, "Hetotrophic", "Heterotrophic"))
merged %>% distinct(xg_pred)

merged <- merged %>% filter(str_detect(tax_name, " "))
merged <- merged %>% filter(!(tax_name %in% c("Calanus finmarchicus", "Hydra vulgaris", 
                                              "Neocalanus flemingeri", "Oikopleura dioica", "Paracyclopina nana")))
merged <- merged %>% filter(!(tax_name %in% c("cellular organisms", "Micromonas <green algae>")))
merged <- merged %>% filter(!(tax_name %in% c("unclassified Acanthoeca", "unclassified Micromonas", "unclassified Phaeocystis")))
merged <- merged %>% filter(tax_name != "unclassified Chrysochromulina")
merged <- merged%>% filter(!(tax_name %in% c("Calanus finmarchicus", "Eucyclops serrulatus", 
                                             "Hydra vulgaris", "Oikopleura dioica", "Paracyclopina nana")))
merged <- merged %>% filter(tax_name != "Pleuromamma xiphias")

merged %>% distinct(tax_name)

sample <- read_csv("ctd_g1_g2_15m_source_Stephen.csv")

head(sample)
sample <- sample %>% select(SAMPLE_ID, DEPTH, LATITUDE)

sample <- sample %>% as.data.frame()

head(merged)
head(sample)

sample <- sample %>% mutate(SAMPLE_ID = str_replace(SAMPLE_ID, "0_2um", "0.2um"))

merged %>% anti_join(sample, by = c("sample_id" = "SAMPLE_ID"))

nrow(merged)
merged <- merged %>% left_join(sample, by = c("sample_id" = "SAMPLE_ID"))
nrow(merged)

pfamsExcluded <- read_csv("bootstrapPfamsExcluded.csv")

head(pfamsExcluded)

pfamsExcluded <- pfamsExcluded %>% select(-`...1`)
pfamsExcluded <- pfamsExcluded %>% gather(col0:col29, key = "boot", value = "excludedPfam")

merged <- merged %>% gather(col0:col29, key = "bootNum", value = "boot_xg_pred")



inc <- read_csv("corePfamSummary_boot_lambert.csv")

head(merged)
head(inc)

merged <- merged %>% semi_join(inc, by = c("sample_id" = "var", "tax_name", "bootNum"))


merged %>% distinct(boot_xg_pred)
merged <- merged %>% mutate(boot_xg_pred = str_c(boot_xg_pred, "otrophic"))
merged <- merged %>% mutate(boot_xg_pred = str_replace(boot_xg_pred, "Hetotrophic", "Heterotrophic"))
merged %>% distinct(boot_xg_pred)

head(merged)

exclude <- merged %>% filter(boot_xg_pred %in% c("Phototrophic", "Heterotrophic")) %>% group_by(LATITUDE, tax_name, bootNum) %>% 
  distinct(boot_xg_pred) %>% summarize(n = n()) %>% filter(n == 2)
exclude <- exclude %>% ungroup() 
exclude

merged <- merged %>% anti_join(exclude, by = c("LATITUDE", "tax_name", "bootNum"))

#get rid of predictions for transcriptomes with less than 70% of core Pfams expressed

#1056/7520 (14.04255)
head(merged)
100*(merged %>% filter(xg_pred != boot_xg_pred) %>% nrow())/(merged %>% nrow())
(merged %>% filter(xg_pred != boot_xg_pred) %>% nrow())
(merged %>% nrow())

merged %>% group_by(xg_pred, boot_xg_pred) %>% summarize(n = n()) %>% 
  mutate(prop = n/(merged %>% nrow())*100)
