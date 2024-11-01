library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/g3Depth/")

merged <- read_csv("g3Depth_speciesWithTrophicPredictions_fall2023_28.csv")

merged_notOfInt <- read_csv("g3Depth_otherThanSpeciesWithTrophicPredictions_fall2023_28.csv")

#downloaded from https://docs.google.com/spreadsheets/d/1lKwvlNjUV8dXHKO5HkMNa2m70QQYydrTf-Fb6PRQ9oE/edit#gid=0
sample <- read_excel("lookup_armbrust_grc_rnaseq_6.xlsx")

sample <- sample %>% filter(str_detect(`Investigator Sample Name`, "D"))
sample <- sample %>% select(2,3,4)

colnames(sample)[1]
colnames(sample)[1] <- "sample"

sample <- sample %>% mutate(stationCast = str_extract(`Add'l ID`, "S[0-9]{1,}C[0-9]{1,}"))
sample <- sample %>% mutate(Depth = str_extract(`Add'l ID`, " [0-9mDCM]{1,}"))
sample <- sample %>% mutate(Depth = str_replace(Depth, " ", ""))
sample <- sample %>% mutate(replicate = str_extract(`Add'l ID`, " [AB]"))
sample <- sample %>% mutate(replicate = str_replace(replicate, " ", ""))

sample <- sample %>% mutate(sample2 = str_c("G3PA.depth.", stationCast, ".", Depth, ".", replicate))

merged %>% anti_join(sample, by = c("sample" = "sample2"))
merged_notOfInt %>% anti_join(sample, by = c("sample" = "sample2"))

nrow(merged)
merged <- merged %>% left_join(sample, by = c("sample" = "sample2"))
nrow(merged)

nrow(merged_notOfInt)
merged_notOfInt <- merged_notOfInt %>% left_join(sample, by = c("sample" = "sample2"))
nrow(merged_notOfInt)

stan <- read_csv("G3Depth_norm_factors.csv", col_names = TRUE)

stan <- stan %>% select(sample, NORM_FACTOR)

stan$sample <- as.character(stan$sample)

merged %>% anti_join(stan, by = c("sample.y" = "sample"))
merged_notOfInt %>% anti_join(stan, by = c("sample.y" = "sample"))

nrow(merged)
merged <- merged %>% left_join(stan, by = c("sample.y" = "sample"))
nrow(merged)

nrow(merged_notOfInt)
merged_notOfInt <- merged_notOfInt %>% left_join(stan, by = c("sample.y" = "sample"))
nrow(merged_notOfInt)

str(merged)
merged <- merged %>% mutate(absoluteCounts = est_counts*NORM_FACTOR)

merged %>% write_csv("speciesWithTrophicModePredictionsAbundance.csv")

str(merged_notOfInt)
merged_notOfInt <- merged_notOfInt %>% mutate(absoluteCounts = est_counts*NORM_FACTOR)

speciesPred <- read_csv("g3Depth_trophicPredictions_cleanedUp_updatedMarferret_marmicroDb2023_noOutliers.csv")

speciesPred %>% anti_join(merged, by = c("taxa" = "tax_name"))

both <- bind_rows(merged %>% mutate(taxaGroup = "speciesWithTrophicMode"), merged_notOfInt %>% mutate(taxaGroup = "Rest"))

both <- both %>% filter(!is.na(tax_name))

both <- both %>% 
  filter(str_detect(tax_name, " ") & !str_detect(tax_name, "unclassified") & !str_detect(tax_name, "uncultured") & 
           !(tax_name %in% c("Chlorella clade", "core chlorophytes", "CS clade", "Elliptochloris clade", "cellular organisms",
                             "Marine_stramenopiles_MASTs group", "Micromonas <green algae>", "Nitzschia <diatoms>", "PX clade", "Rhodomonas <cryptomonads>",
                             "Calanus finmarchicus", "Hydra vulgaris", "cellular organisms", "core chlorophytes",
                             "uncultured eukaryote", "Micromonas <green algae>", "Rhodomonas <cryptomonads>", "PX clade",
                             "CS clade", "Chlorella clade", "Elliptochloris clade", "Nitzschia <diatoms>", "Marine_stramenopiles_MASTs group",
                             "Neocalanus flemingeri", "Oikopleura dioica", "Paracyclopina nana", "Vannella robusta", "Pleuromamma xiphias")))

both <- both %>% mutate(genus = str_extract(tax_name, "[A-Za-z0-9]{1,}"))
both %>% filter(is.na(genus))

exclude <- c("Apostichopus", "Caulerpa", "Acartia", "Adineta", "Amphimedon", "Anisakis", "Aplysia", "Bugula", "Calanus", "Capitella", "Chara", "Chlorokybus", "Chondrus", 
             "Cladosiphon", "Compsopogon", "Crassostrea", "Daphnia", "Debaryomyces", "Dibothriocephalus", "Ectocarpus", "Eucyclops", 
             "Eurytemora", "Gracilariopsis", "Helobdella", "Hippoglossus", "Homarus", "Klebsormidium", "Labidocera", "Lepeophtheirus", 
             "Lingula", "Nemacystus", "Nematostella", "Neopyropia", "Octopus", "Opisthorchis", "Porphyra", "Saccharina", "Strongylocentrotus", 
             "Tigriopus", "Trichuris", "Tursiops", "Ulva", "Undaria", "Vaucheria")

both <- both %>% filter(!(genus %in% exclude))

dinos <- read_csv("../dinoflagellateSpecies.csv")

both %>% anti_join(dinos, by = c("tax_name"))

nrow(both)
both <- both %>% left_join(dinos, by = c("tax_name"))
nrow(both)

both <- both %>% mutate(absoluteCounts_dinoCorr = ifelse(group == "Dinoflagellate", absoluteCounts/6.4, absoluteCounts))

withTrophMode <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts_dinoCorr = sum(absoluteCounts_dinoCorr)) %>% 
  filter(taxaGroup == "speciesWithTrophicMode") %>% select(absoluteCounts_dinoCorr)

rest <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts_dinoCorr = sum(absoluteCounts_dinoCorr)) %>% 
  filter(taxaGroup == "Rest") %>% select(absoluteCounts_dinoCorr)


#55.35 28 species
round(100*withTrophMode$absoluteCounts_dinoCorr/(withTrophMode$absoluteCounts_dinoCorr+rest$absoluteCounts_dinoCorr),2)

