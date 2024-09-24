library(tidyverse)
library(readxl)

setwd("~/Dropbox/grad/research/alohaDiel/")

dat <- read_csv("G1_diel_allSamples_processed_updatedMarferret_marmicroDb2023_noOutliers_speciesWithTrophicPredictions_fall2023_28.csv")

notOfInt <- read_csv("G1_diel_allSamples_processed_updatedMarferret_marmicroDb2023_noOutliers_otherThanSpeciesWithTrophicPredictions_fall2023_28.csv")

stan <- read_csv("diel1_SUM_norm_factors.csv")

dat %>% anti_join(stan, by = c("sample"))
notOfInt %>% anti_join(stan, by = c("sample"))

head(stan)
colnames(stan)[2] <- "convFactor"

nrow(dat)
dat <- dat %>% left_join(stan, by = c("sample"))
nrow(dat)

nrow(notOfInt)
notOfInt <- notOfInt %>% left_join(stan, by = c("sample"))
nrow(notOfInt)

str(dat)
str(notOfInt)

dat <- dat %>% mutate(absoluteCounts = est_counts*convFactor)
notOfInt <- notOfInt %>% mutate(absoluteCounts = est_counts*convFactor)

dat %>% select(tax_name, sample, absoluteCounts) %>% write_csv("speciesWithTrophicModePredictionsAbundance.csv")

both <- bind_rows(dat %>% mutate(taxaGroup = "speciesWithTrophicMode"), notOfInt %>% mutate(taxaGroup = "Rest"))

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

withTrophMode <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts = sum(absoluteCounts)) %>% 
  filter(taxaGroup == "speciesWithTrophicMode") %>% select(absoluteCounts)

rest <- both %>% group_by(taxaGroup) %>% summarize(absoluteCounts = sum(absoluteCounts)) %>% 
  filter(taxaGroup == "Rest") %>% select(absoluteCounts)

#59.24 28 species
round(100*withTrophMode$absoluteCounts/(withTrophMode$absoluteCounts+rest$absoluteCounts),2)







